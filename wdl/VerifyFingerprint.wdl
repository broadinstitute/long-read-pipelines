version 1.0

import "tasks/Structs.wdl"

import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF

workflow VerifyFingerprint {
    input {
        File aligned_bam
        File aligned_bai

        File ref_map_file
        String fingerprint_dir

        String gcs_out_root_dir
    }

    parameter_meta {
        aligned_bam:        "GCS path to aligned BAM file"
        aligned_bai:        "GCS path to aligned BAM file index"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"
        fingerprint_dir:    "GCS directory in which fingerprint files are stored"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/VerifyFingerprint"

    call FindFingerprint {
        input:
            aligned_bam     = aligned_bam,
            fingerprint_dir = fingerprint_dir
    }

    call CheckFingerprint {
        input:
            aligned_bam     = aligned_bam,
            aligned_bai     = aligned_bai,
            fingerprint_vcf = FindFingerprint.fingerprint_vcf,
            haplotype_map   = ref_map['haplotype_map']
    }

    call FF.FinalizeToFile as FinalizeFingerprintSummaryMetrics { input: outdir = outdir, file = CheckFingerprint.summary_metrics }
    call FF.FinalizeToFile as FinalizeFingerprintDetailMetrics { input: outdir = outdir, file = CheckFingerprint.detail_metrics }

    output {
        Float lod_expected_sample = CheckFingerprint.metrics_map['LOD_EXPECTED_SAMPLE']
        File fingerprint_metrics = FinalizeFingerprintSummaryMetrics.gcs_path
        File fingerprint_details = FinalizeFingerprintDetailMetrics.gcs_path
    }
}

task FindFingerprint {
    input {
        String aligned_bam
        String fingerprint_dir

        Array[String] filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']

        RuntimeAttr? runtime_attr_override
    }

    String fpdir = sub(fingerprint_dir, "/$", "")
    Int disk_size = 1

    command <<<
        set -x

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        SM=$(samtools view -H ~{aligned_bam} | grep '^@RG' | head -1 | sed 's/\t/\n/g' | grep SM | sed 's/SM://')

        gsutil cp ~{fpdir}/$SM.vcf fingerprint.vcf

        cat fingerprint.vcf | \
            grep -v -e 'placeholder' ~{true='-e' false='' length(filter) > 0} ~{sep=" -e " filter} \
            > fingerprint.fixed.vcf
    >>>

    output {
        File fingerprint_vcf = "fingerprint.fixed.vcf"
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:             2,
        mem_gb:                4,
        disk_gb:               disk_size,
        boot_disk_gb:          10,
        preemptible_tries:     3,
        max_retries:           2,
        docker:                "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                   select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:        select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:           select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:            select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                select_first([runtime_attr.docker, default_attr.docker])
    }
}

task CheckFingerprint {
    input {
        String aligned_bam
        String aligned_bai
        File fingerprint_vcf
        File haplotype_map

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size([fingerprint_vcf, haplotype_map], "GB"))
    String prefix = basename(aligned_bam, ".bam")

    command <<<
        set -x

        gatk CheckFingerprint \
            --INPUT ~{aligned_bam} \
            --GENOTYPES ~{fingerprint_vcf} \
            --HAPLOTYPE_MAP ~{haplotype_map} \
            --OUTPUT ~{prefix}

        grep -v '^#' ~{prefix}.fingerprinting_summary_metrics | grep -A1 READ_GROUP | awk '
            {
                for (i=1; i<=NF; i++)  {
                    a[NR,i] = $i
                }
            }
            NF>p { p = NF }
            END {
                for(j=1; j<=p; j++) {
                    str=a[1,j]
                    for(i=2; i<=NR; i++){
                        str=str" "a[i,j];
                    }
                    print str
                }
            }' | sed 's/ /\t/' > metrics_map.txt

        mv ~{prefix}.fingerprinting_summary_metrics ~{prefix}.fingerprinting_summary_metrics.txt
        mv ~{prefix}.fingerprinting_detail_metrics ~{prefix}.fingerprinting_detail_metrics.txt
    >>>

    output {
        File summary_metrics = "~{prefix}.fingerprinting_summary_metrics.txt"
        File detail_metrics = "~{prefix}.fingerprinting_detail_metrics.txt"
        Map[String, String] metrics_map = read_map("metrics_map.txt")
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:             2,
        mem_gb:                4,
        disk_gb:               disk_size,
        boot_disk_gb:          10,
        preemptible_tries:     3,
        max_retries:           2,
        docker:                "us.gcr.io/broad-gatk/gatk:4.2.0.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                   select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:        select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:           select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:            select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                select_first([runtime_attr.docker, default_attr.docker])
    }
}
