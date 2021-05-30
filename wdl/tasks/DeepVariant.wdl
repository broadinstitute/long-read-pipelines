version 1.0

#######################################################
# This pipeline calls small variants using DeepVariant.
#######################################################

import "Structs.wdl"

task DeepVariant {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fai

        String model_class
        String chr

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: "input BAM from which to call variants"
        bai: "index accompanying the BAM"

        ref_fasta: "reference to which the BAM was aligned to"
        ref_fai:   "index accompanying the reference"

        model_class: "class of model to be applied; currently only 'PACBIO' is accepted"
        chr: "chromsome on which to call variants"
    }

    Int disk_size = ceil(size(bam, "GB")) + 50
    String prefix = basename(bam, ".bam")

    command <<<
        # example from https://github.com/google/deepvariant/blob/r0.8/docs/deepvariant-quick-start.md
        set -euo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        /opt/deepvariant/bin/run_deepvariant \
          --model_type=~{model_class} \
          --ref=~{ref_fasta} \
          --reads=~{bam} \
          --output_vcf=~{prefix}.deepvariant.~{chr}.vcf.gz \
          --output_gvcf=~{prefix}.deepvariant.~{chr}.g.vcf.gz \
          --num_shards=${num_core} \
          --regions=~{chr}
    >>>

    output {
        # save both VCF and gVCF
        File vcf = "~{prefix}.deepvariant.~{chr}.vcf.gz"
        File vcf_tbi = "~{prefix}.deepvariant.~{chr}.vcf.gz.tbi"
        File gvcf = "~{prefix}.deepvariant.~{chr}.g.vcf.gz"
        File gvcf_tbi = "~{prefix}.deepvariant.~{chr}.g.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          12,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "gcr.io/deepvariant-docker/deepvariant:1.1.0-gpu"
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
        cpuPlatform:            "Intel Skylake"
    }
}

task PEPPER {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        String chr
        String preset

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:             "input BAM from which to call variants"
        bai:             "index accompanying the BAM"

        ref_fasta:       "reference to which the BAM was aligned to"
        ref_fasta_fai:   "index accompanying the reference"

        chr:             "chr on which to call variants"
        preset:          "calling preset (CCS, ONT)"
    }

    Int disk_size = ceil(size(bam, "GB")) + 50
    String prefix = basename(bam, ".bam") + ".deepvariant_pepper"
    String mode = if preset == "CCS" then "--ccs" else "--ont"

    command <<<
        # example from https://github.com/kishwarshafin/pepper/blob/r0.4/docs/pipeline_docker/ONT_variant_calling.md
        set -euxo pipefail

        touch ~{bai}

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)
        SM=$(samtools view -H ~{bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g')

        samtools view --no-PG -H ~{bam} | sed 's/\tPU:.\+//g' > header.txt
        samtools reheader -P header.txt ~{bam} > fixed.bam
        samtools index fixed.bam

        samtools view -H fixed.bam

        run_pepper_margin_deepvariant call_variant \
            -b "fixed.bam" \
            -f "~{ref_fasta}" \
            -o ./ \
            -p "~{prefix}" \
            -r "~{chr}" \
            -t ${num_core} \
            -s $SM \
            --gvcf \
            --phased_output \
            ~{mode}

        touch ~{prefix}.phased.vcf.gz ~{prefix}.vcf.gz ~{prefix}.g.vcf.gz
        touch ~{prefix}.phased.vcf.gz.tbi ~{prefix}.vcf.gz.tbi ~{prefix}.g.vcf.gz.tbi
        touch ~{prefix}.visual_report.html ~{prefix}.phaseset.bed
    >>>

    output {
        # save both VCF and gVCF
        File phased_vcf = "~{prefix}.phased.vcf.gz"
        File phased_vcf_tbi = "~{prefix}.phased.vcf.gz.tbi"

        File vcf = "~{prefix}.vcf.gz"
        File vcf_tbi = "~{prefix}.vcf.gz.tbi"
        File gvcf = "~{prefix}.g.vcf.gz"
        File gvcf_tbi = "~{prefix}.g.vcf.gz.tbi"

        File report = "~{prefix}.visual_report.html"
        File phaseset_bed = "~{prefix}.phaseset.bed"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          12,
        mem_gb:             72,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "kishwars/pepper_deepvariant:r0.4.1"
        #docker:             "us.gcr.io/broad-dsp-lrma/lr-dvpepper:r0.4.1"
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
        cpuPlatform:            "Intel Skylake"
    }
}
