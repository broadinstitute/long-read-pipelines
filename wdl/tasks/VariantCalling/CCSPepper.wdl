version 1.0

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
        pepper_threads: "Number of threads for Pepper"
        pepper_memory: "Memory for Pepper"
        dv_threads: "Number of threads for DeepVariant"
        dv_memory: "Memory for DeepVariant"
        zones: "select which zone (GCP) to run this task"
    }

    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        Int pepper_threads
        Int pepper_memory

        Int dv_threads
        Int dv_memory

        Array[String] zones = ["us-central1-c", "us-central1-f", "us-central1-a", "us-central1-b"]
    }

    call Pepper as get_hap_tagged_bam {
        input:
            bam = bam,
            bai = bai,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            threads = pepper_threads,
            memory = pepper_memory,
            zones = zones
    }

    call DV as deep_variant {
        input:
            bam = get_hap_tagged_bam.hap_tagged_bam,
            bai = get_hap_tagged_bam.hap_tagged_bai,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            threads = dv_threads,
            memory = dv_memory,
            zones = zones
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
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        Int threads
        Int memory
        Array[String] zones

        RuntimeAttr? runtime_attr_override
    }

	Int disk_size = 100 + 2*ceil(size(bam, "GB")) + 2*ceil(size(bai, "GB")) + 2*ceil(size(ref_fasta, "GB")) + 2*ceil(size(ref_fasta_fai, "GB"))

    String output_root = "/cromwell_root/pepper_output"

    String prefix = basename(bam, ".bam") + ".pepper"

    command <<<
        # avoid the infamous pipefail 141 https://stackoverflow.com/questions/19120263
        set -eux
        SM=$(samtools view -H ~{bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g')

        set -euxo pipefail

        touch ~{bai}
        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        mkdir -p "~{output_root}"

        # no gVCF as it Pepper simply doesn't produce gVCF on CCS data
        run_pepper_margin_deepvariant \
            call_variant \
            -b ~{bam} \
            -f ~{ref_fasta} \
            -t "${num_core}" \
            -s "${SM}" \
            -o "~{output_root}" \
            -p "~{prefix}" \
            --phased_output \
            --ccs

        find "~{output_root}/" -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g' \
            > "~{output_root}/dir_structure.txt"

        if [[ -f "~{output_root}/intermediate_files/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam" ]]; then
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
        cpu_cores:          threads,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "kishwars/pepper_deepvariant:r0.4.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones:                  "${sep=' ' zones}"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task DV {

    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        Int threads
        Int memory
        Array[String] zones

        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, ".bam") + ".deepvariant"
    String output_root = "/cromwell_root/dv_output"

    Int bam_sz = ceil(size(bam, "GB"))
    Boolean is_big_bam = bam_sz > 100
    Int inflation_factor = if (is_big_bam) then 10 else 5
    Int minimal_disk = 1000
	Int disk_size = if inflation_factor * bam_sz > minimal_disk then inflation_factor * bam_sz else minimal_disk

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        mkdir -p "~{output_root}"

        export MONITOR_MOUNT_POINT="/cromwell_root/"
        bash vm_local_monitoring_script.sh &> resources.log &
        job_id=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')

        /opt/deepvariant/bin/run_deepvariant \
            --model_type=PACBIO \
            --ref=~{ref_fasta} \
            --reads=~{bam} \
            --output_vcf="~{output_root}/~{prefix}.vcf.gz" \
            --output_gvcf="~{output_root}/~{prefix}.g.vcf.gz" \
            --num_shards="${num_core}" \
            --use_hp_information || cat resources.log
        if ps -p "${job_id}" > /dev/null; then kill "${job_id}"; fi

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
        cpu_cores:          threads,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
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
        zones:                  "${sep=' ' zones}"
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

        Int memory
        Array[String] zones

        RuntimeAttr? runtime_attr_override
    }

    Int bam_sz = ceil(size(bam, "GB"))
	Int disk_size = if bam_sz > 200 then 2*bam_sz else bam_sz + 200

    Int cores = 64

    String prefix = basename(bam, ".bam") + ".pepper"
    String output_root = "/cromwell_root/margin_output"

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        mkdir -p "~{output_root}" "~{output_root}/logs"
        touch ~{bai}

        # note the -M option was suggested by an author of margin
        # it's unclear which phasedBAM one should use: this, or the one generated from the Pepper step
        margin phase \
            ~{bam} \
            ~{ref_fasta} \
            ~{unphased_vcf} \
            /opt/margin_dir/params/misc/allParams.phase_vcf.json \
            -t "${num_core}" \
            -M \
            -o "~{output_root}/~{prefix}" \
            2>&1 | tee "~{output_root}/logs/5_margin_phase_vcf.log"

        bgzip -c "~{output_root}/~{prefix}".phased.vcf > "~{output_root}/~{prefix}".phased.vcf.gz && \
            tabix -p vcf "~{output_root}/~{prefix}".phased.vcf.gz
    >>>


    output {
        File phaseset_bed = "~{output_root}/~{prefix}.phaseset.bed"
        File phasedVCF  = "~{output_root}/~{prefix}.phased.vcf.gz"
        File phasedtbi  = "~{output_root}/~{prefix}.phased.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cores,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "kishwars/pepper_deepvariant:r0.4.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        zones:                  "${sep=' ' zones}"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
