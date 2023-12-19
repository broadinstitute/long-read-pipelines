version 1.0

#######################################################
# This pipeline calls small variants using DeepVariant.
#######################################################

import "../../structs/Structs.wdl"


workflow CCSPepper {

    meta {
        description: "Workflow for getting haplotagged BAM, VCF and gVCF from DV-pepper. Note VCF is un-phased."
    }

    parameter_meta {
        bam: "Input BAM file"
        bai: "Input BAM index file"
        ref_fasta: "Reference fasta file"
        ref_fasta_fai: "Reference fasta index file"
        regions:         "genomic regions on which to call variants"
    }

    input {
        File bam
        File bai
        File ref_fasta
        File ref_fasta_fai
        Array[String] regions
    }

    call Pepper as get_hap_tagged_bam {
        input:
            bam = bam,
            bai = bai,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            regions = regions
    }

    call DV as deep_variant {
        input:
            bam = get_hap_tagged_bam.hap_tagged_bam,
            bai = get_hap_tagged_bam.hap_tagged_bai,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            regions = regions
    }

    output {
        File VCF        = deep_variant.VCF
        File VCF_tbi    = deep_variant.VCF_tbi

        File gVCF       = deep_variant.gVCF
        File gVCF_tbi   = deep_variant.gVCF_tbi

        File hap_tagged_bam = get_hap_tagged_bam.hap_tagged_bam
        File hap_tagged_bai = get_hap_tagged_bam.hap_tagged_bai
    }
}

task Pepper {
    parameter_meta {
        bam:             {description: "input BAM from which to call variants",
                          localization_optional: true}
        bai:             {description: "index accompanying the BAM",
                          localization_optional: true}
        ref_fasta:       "reference to which the BAM was aligned"
        ref_fasta_fai:   "index accompanying the reference"
        regions:         "genomic regions on which to call variants"

        runtime_attr_override: "override the default runtime attributes"
    }

    input {
        File bam
        File bai
        File ref_fasta
        File ref_fasta_fai
        Array[String] regions

        RuntimeAttr? runtime_attr_override
    }

    Int bam_sz = ceil(size(bam, "GB"))
    Int disk_size = if bam_sz > 200 then 2*bam_sz else bam_sz + 200

    String output_root = "pepper_output"
    String prefix = basename(bam, ".bam") + ".pepper"

    command <<<
        set -euxo pipefail

        samtools view -h -1 --write-index -o "local.bam##idx##local.bam.bai" -X "~{bai}" "~{bam}" "~{sep='" "' regions}"
        SAMPLE=$(samtools view -H local.bam | sed '/^@RG/!d;s/.*	SM:\([^	]*\).*/\1/' | sed '2,$d')

        mkdir -p "~{output_root}"

        # no gVCF as it Pepper simply doesn't produce gVCF on CCS data
        run_pepper_margin_deepvariant \
            call_variant \
            -b local.bam \
            -f "~{ref_fasta}" \
            -t "$(nproc)" \
            -s "$SAMPLE" \
            -o "~{output_root}" \
            -p "~{prefix}" \
            --phased_output \
            --ccs

        find "~{output_root}/" -print | sed -e 's+[^/]*/+|____+g;s+____|+ |+g' \
            > "~{output_root}/dir_structure.txt"

        if [ -f "~{output_root}/intermediate_files/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam" ]; then
            mv "~{output_root}/intermediate_files/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam" \
               "~{output_root}/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam"
            mv "~{output_root}/intermediate_files/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai" \
               "~{output_root}/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai"
        fi
    >>>

    output {
        File hap_tagged_bam = "~{output_root}/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam"
        File hap_tagged_bai = "~{output_root}/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai"

        # maybe less useful
        File output_dir_structure = "~{output_root}/dir_structure.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "kishwars/pepper_deepvariant:r0.4.1"
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

task DV {
    parameter_meta {
        bam:             {description: "input BAM from which to call variants",
                          localization_optional: true}
        bai:             {description: "index accompanying the BAM",
                          localization_optional: true}
        ref_fasta:       "reference to which the BAM was aligned"
        ref_fasta_fai:   "index accompanying the reference"
        regions:         "genomic regions on which to call variants"

        runtime_attr_override: "override the default runtime attributes"
    }

    input {
        File bam
        File bai
        File ref_fasta
        File ref_fasta_fai
        Array[String] regions

        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, ".bam") + ".deepvariant"
    String output_root = "dv_output"

    Int bam_sz = ceil(size(bam, "GB"))
    Boolean is_big_bam = bam_sz > 100
    Int inflation_factor = if (is_big_bam) then 10 else 5
    Int minimal_disk = 1000
    Int disk_size = if inflation_factor * bam_sz > minimal_disk then inflation_factor * bam_sz else minimal_disk

    command <<<
        set -euxo pipefail

        samtools view -h -1 --write-index -o "local.bam##idx##local.bam.bai" -X "~{bai}" "~{bam}" "~{sep='" "' regions}"

        export MONITOR_MOUNT_POINT=$PWD
        bash vm_local_monitoring_script.sh &> resources.log &
        monitor_pid=$!

        mkdir -p "~{output_root}"
        /opt/deepvariant/bin/run_deepvariant \
            --model_type=PACBIO \
            --ref="~{ref_fasta}" \
            --reads=local.bam \
            --output_vcf="~{output_root}/~{prefix}.vcf.gz" \
            --output_gvcf="~{output_root}/~{prefix}.g.vcf.gz" \
            --num_shards="$(nproc)" \
            --use_hp_information
        dv_result=$?

        if ps -p "$monitor_pid" > /dev/null; then kill "$monitor_pid"; fi

        if [ $dv_result -ne 0 ]; then
            echo "******deepvariant failed with result code $dv_result"
            echo "******monitor log follows"
            cat resources.log
            exit $dv_result
        fi

        find "~{output_root}/" -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g' \
            > "~{output_root}/dir_structure.txt"
    >>>

    output {

        File resouce_monitor_log = "resources.log"

        File output_dir_structure = "~{output_root}/dir_structure.txt"

        File VCF        = "~{output_root}/~{prefix}.vcf.gz"
        File VCF_tbi    = "~{output_root}/~{prefix}.vcf.gz.tbi"

        File gVCF       = "~{output_root}/~{prefix}.g.vcf.gz"
        File gVCF_tbi   = "~{output_root}/~{prefix}.g.vcf.gz.tbi"

        File visual_report_html = "~{output_root}/~{prefix}.visual_report.html"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-deepvariant:1.3.0"
        # docker:             "google/deepvariant:1.2.0-gpu"  # kept here to remind ourselves, occassionally, to review if it's better with GPU
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

task MarginPhase {

    meta {
        description: "Generates phased VCF. Note this runs fast so no need to parallize."
    }

    input {
        File bam
        File bai

        File unphased_vcf
        File? unphased_vcf_tbi

        File ref_fasta
        File ref_fasta_fai

        String? output_bucket

        RuntimeAttr? runtime_attr_override
    }

    Int bam_sz = ceil(size(bam, "GB"))
    Int disk_size = if bam_sz > 200 then 2*bam_sz else bam_sz + 200

    String prefix = basename(bam, ".bam") + ".pepper"
    String output_root = "margin_output"

    command <<<
        set -euxo pipefail

        mkdir -p "~{output_root}/logs"

        # note the -M option was suggested by an author of margin
        # it's unclear which phasedBAM one should use: this, or the one generated from the Pepper step
        margin phase \
            "~{bam}" \
            "~{ref_fasta}" \
            "~{unphased_vcf}" \
            /opt/margin_dir/params/misc/allParams.phase_vcf.json \
            -t "$(nproc)" \
            -M \
            -o "~{output_root}/~{prefix}" \
            2>&1 | tee "~{output_root}/logs/5_margin_phase_vcf.log"

        vcfName="~{output_root}/~{prefix}.phased.vcf.gz"
        bgzip -c "~{output_root}/~{prefix}".phased.vcf > "$vcfName"
        tabix -p vcf "$vcfName"

        tbiName="${vcfName}.tbi"
        bedName="~{output_root}/~{prefix}.phaseset.bed"
        if ~{defined(output_bucket)}; then
            outDir=$(echo "~{output_bucket}" | sed 's+/?$+/+')
            gcloud storage cp "$vcfName" "$tbiName" "$bedName" "$outDir"
        fi
    >>>


    output {
        File phaseset_bed = "$bedName"
        File phasedVCF  = "$vcfName"
        File phasedtbi  = "$tbiName"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          64,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-marginphase:0.1.1"
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
