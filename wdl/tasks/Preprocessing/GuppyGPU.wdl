version 1.0

##########################################################################################
# A workflow that runs the Guppy basecaller on ONT FAST5 files.
# - The docker tag number will match the version of Guppy that is being run. You can change
#   this value to run a different version of Guppy. Currently supports... [3.5.2, 3.6.0, 4.0.14]
# - All fast5 files within the given GCS dir, gcs_fast5_dir, will be processed
# - Takes a few hours to process 130GB. Best guess is that the processing time scales
#   linearly but untested.
##########################################################################################

import "../Utility/Utils.wdl" as Utils
import "../Utility/ONTUtils.wdl" as ONT
import "../../structs/Structs.wdl"

workflow GuppyGPU {

    meta {
        description: "Run Guppy basecaller on ONT FAST5 files. The docker tag number will match the version of Guppy that is being run. You can change this value to run a different version of Guppy. Currently supports... [3.5.2, 3.6.0, 4.0.14]. All fast5 files within the given GCS dir, gcs_fast5_dir, will be processed. Takes a few hours to process 130GB. Best guess is that the processing time scales linearly but untested."
    }
    parameter_meta {
        gcs_fast5_dir: "GCS path to a directory containing ONT FAST5 files."
        config: "Guppy config file."
        instrument: "Optional. Instrument used for sequencing. Default is 'unknown'."
        flow_cell_id: "Optional. Flow cell ID used for sequencing. Default is 'unknown'."
        protocol_run_id: "Optional. Protocol run ID used for sequencing. Default is 'unknown'."
        sample_name: "Optional. Sample name used for sequencing. Default is 'unknown'."
        num_shards: "Optional. Number of shards to use for parallelization. Default is 1 + ceil(length(read_lines(ListFast5s.manifest))/100)."
        gcs_out_root_dir: "GCS path to a directory where the output will be written."
    }

    input {
        String gcs_fast5_dir

        String config

        String instrument = "unknown"
        String flow_cell_id = "unknown"
        String protocol_run_id = "unknown"
        String sample_name = "unknown"
        Int? num_shards

        String gcs_out_root_dir
    }

    call Utils.Timestamp as TimestampStarted { input: dummy_dependencies = [config] }

    call ListFast5s { input: gcs_fast5_dir = gcs_fast5_dir }

    Int ns = 1 + select_first([num_shards, ceil(length(read_lines(ListFast5s.manifest))/100)])
    call ONT.PartitionManifest as PartitionFast5Manifest { input: manifest = ListFast5s.manifest, N = ns }

    scatter (chunk_index in range(length(PartitionFast5Manifest.manifest_chunks))) {
        call Basecall {
            input:
                fast5_files  = read_lines(PartitionFast5Manifest.manifest_chunks[chunk_index]),
                config       = config,
        }
    }

    call Utils.Timestamp as TimestampStopped { input: dummy_dependencies = Basecall.sequencing_summary }

    call MakeSequencingSummary { input: sequencing_summaries = Basecall.sequencing_summary }

    call MakeFinalSummary {
        input:
            instrument      = instrument,
            flow_cell_id    = flow_cell_id,
            sample_id       = sample_name,
            protocol_run_id = protocol_run_id,
            started         = TimestampStarted.timestamp,
            stopped         = TimestampStopped.timestamp
    }

    call FinalizeBasecalls {
        input:
            pass_bams          = flatten(Basecall.pass_bams),
            sequencing_summary = MakeSequencingSummary.sequencing_summary,
            final_summary      = MakeFinalSummary.final_summary,
            outdir             = gcs_out_root_dir
    }

    output {
        String gcs_dir = FinalizeBasecalls.gcs_dir
    }
}

task ListFast5s {
    input {
        String gcs_fast5_dir

        RuntimeAttr? runtime_attr_override
    }

    String indir = sub(gcs_fast5_dir, "/$", "")

    command <<<
        gsutil ls "~{indir}/**.fast5" > fast5_files.txt
    >>>

    output {
        File manifest = "fast5_files.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task Basecall {
    input {
        Array[File] fast5_files
        String config = "dna_r10.4.1_e8.2_400bps_sup.cfg"

        RuntimeAttr? runtime_attr_override
    }

#    dna_r10.3_450bps_sup.cfg
#    dna_r10.4.1_e8.2_260bps_sup.cfg
#    dna_r10.4.1_e8.2_400bps_sup.cfg
#    dna_r10.4_e8.1_sup.cfg
#    dna_r9.4.1_450bps_modbases_5mc_cg_sup_prom.cfg
#    dna_r9.4.1_450bps_sup.cfg
#    dna_r9.4.1_450bps_sup_prom.cfg
#    dna_r9.4.1_e8.1_sup.cfg
#    dna_r9.5_450bps.cfg

    Int disk_size = 3 * ceil(size(fast5_files, "GB"))

    command <<<
        set -x

        guppy_basecaller \
            -r \
            -i /cromwell_root/ \
            -s guppy_output/ \
            -c ~{config} \
            -x "cuda:all" \
            --bam_out
    >>>

    output {
        Array[File] pass_bams = glob("guppy_output/pass/*.bam")
        File sequencing_summary = "guppy_output/sequencing_summary.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       30,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-guppy:6.5.7"
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
        gpuType:                "nvidia-tesla-p100"
        gpuCount:               1
        nvidiaDriverVersion:    "418.152.00"
        zones:                  ["us-central1-c", "us-central1-f", "us-east1-b", "us-east1-c", "us-west1-a", "us-west1-b"]
        cpuPlatform:            "Intel Haswell"
    }
}

task MakeSequencingSummary {
    input {
        Array[File] sequencing_summaries

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3*ceil(size(sequencing_summaries, "GB"))

    command <<<
        set -euxo pipefail

        head -1 ~{sequencing_summaries[0]} > sequencing_summary.txt

        while read p; do
            awk 'NR > 1 { print }' "$p" >> sequencing_summary.txt
        done <~{write_lines(sequencing_summaries)}
    >>>

    output {
        File sequencing_summary = "sequencing_summary.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task MakeFinalSummary {
    input {
        String instrument
        String sample_id
        String flow_cell_id
        String protocol_run_id
        String started
        String stopped

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1

    command <<<
        set -euxo pipefail

        echo 'instrument=~{instrument}' > final_summary.txt
        echo 'flow_cell_id=~{flow_cell_id}' >> final_summary.txt
        echo 'sample_id=~{sample_id}' >> final_summary.txt
        echo 'protocol_run_id=~{protocol_run_id}' >> final_summary.txt
        echo 'started=~{started}' >> final_summary.txt
        echo 'acquisition_stopped=~{stopped}' >> final_summary.txt
        echo 'processing_stopped=~{stopped}' >> final_summary.txt
        echo 'basecalling_enabled=1' >> final_summary.txt
        echo 'sequencing_summary_file=sequencing_summary.txt' >> final_summary.txt
    >>>

    output {
        File final_summary = "final_summary.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task FinalizeBasecalls {
    input {
        Array[String] pass_bams
        File sequencing_summary
        File final_summary

        String outdir

        RuntimeAttr? runtime_attr_override
    }

    String gcs_output_dir = sub(outdir + "/", "/+$", "")

    command <<<
        set -euxo pipefail

        PASS_BAMS="~{write_lines(pass_bams)}"

        cat $PASS_BAMS | gsutil -m cp -I ~{outdir}

        gsutil cp ~{sequencing_summary} ~{outdir}
        gsutil cp ~{final_summary} ~{outdir}
    >>>

    output {
        String gcs_dir = gcs_output_dir
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
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
