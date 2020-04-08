version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/ShardUtils.wdl" as SU
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/Finalize.wdl" as FF

workflow PBCCSOnlySingleFlowcell {
    input {
        String gcs_input_dir

        String? sample_name
        Int num_reads_per_split = 100000

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = gcs_input_dir }

    scatter (subread_bam in FindBams.subread_bams) {
        call PB.GetRunInfo { input: subread_bam = subread_bam }

        String SM  = select_first([sample_name, GetRunInfo.run_info["SM"]])
        String PU  = GetRunInfo.run_info["PU"]
        String ID  = PU
        String DIR = SM + "." + ID

        call SU.IndexUnalignedBam { input: input_bam = subread_bam }
        call SU.MakeReadNameManifests { input: input_bri = IndexUnalignedBam.bri }

        scatter (manifest in MakeReadNameManifests.manifest_chunks) {
            call CCS {
                input:
                    input_bam          = subread_bam,
                    input_bri          = IndexUnalignedBam.bri,
                    read_name_manifest = manifest
            }
        }

        call AR.MergeBams as MergeChunks { input: bams = CCS.consensus }

        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report }
    }

    call AR.MergeBams as MergeRuns { input: bams = MergeChunks.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }

    call PB.MergeCCSReports as MergeAllCCSReports { input: reports = MergeCCSReports.report }

    ##########
    # Finalize
    ##########

#    call FF.FinalizeToDir as FinalizeMergedRuns {
#        input:
#            files = [ MergeRuns.merged_bam ],
#            outdir = outdir + "/" + DIR[0] + "/alignments"
#    }
#
#    call FF.FinalizeToDir as FinalizeCCSMetrics {
#        input:
#            files = [ MergeAllCCSReports.report ],
#            #files = [ MergeAllCCSWithClassesReports.report, MergeAllCCSWithClassesClasses.classes ],
#            outdir = outdir + "/" + DIR[0] + "/metrics/ccs_metrics"
#    }
}

task CCS {
    input {
        String input_bam
        File input_bri
        File read_name_manifest

        Int min_passes = 3
        Float min_snr = 2.5
        Int min_length = 10
        Int max_length = 50000
        Float min_rq = 0.99

        Int cpus = 4
        Int mem = 8

        RuntimeAttr? runtime_attr_override
    }

    Array[String] read_names = read_lines(read_name_manifest)
    Int disk_size = 10 + 4*ceil(size(input_bri, "GB"))

    command <<<
        set -o pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        ((samtools view -H ~{input_bam}) && (sed 's/,/\n/g' ~{read_name_manifest} | sort -t'/' -n -k2 | bri get -i ~{input_bri} ~{input_bam})) | samtools view -b > subreads.bam

        ccs --min-passes ~{min_passes} \
            --min-snr ~{min_snr} \
            --min-length ~{min_length} \
            --max-length ~{max_length} \
            --min-rq ~{min_rq} \
            --num-threads ~{cpus} \
            --report-file ccs_report.txt \
            subreads.bam \
            ccs_unmapped.bam
    >>>

    output {
        File consensus = "ccs_unmapped.bam"
        File report = "ccs_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             mem,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  4,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-bri:0.1.19",
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

