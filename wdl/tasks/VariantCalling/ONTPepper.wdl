version 1.0

import "../../structs/Structs.wdl"

task Pepper {

    meta {
        description: "A 1-stop shop task offered by Pepper for ONT data.  This pipeline calls small variants using DeepVariant."
    }

    parameter_meta {
        bam: "The input bam file."
        bai: "The input bam index file."
        ref_fasta: "The reference fasta file."
        ref_fasta_fai: "The reference fasta index file."
        threads: "The number of threads to use."
        memory: "The amount of memory to use."
        zones: "select which zone (GCP) to run this task"
        runtime_attr_override: "override the default runtime attributes"
    }

    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        Int threads
        Int memory

        Array[String] zones = ["us-central1-c", "us-central1-f", "us-central1-a", "us-central1-b"]

        RuntimeAttr? runtime_attr_override
    }

    Int bam_sz = ceil(size(bam, "GB"))
    Boolean is_big_bam = bam_sz > 100
    Int inflation_factor = if (is_big_bam) then 10 else 5
    Int minimal_disk = 1000
	Int disk_size = if inflation_factor * bam_sz > minimal_disk then inflation_factor * bam_sz else minimal_disk

    String output_root = "/cromwell_root/pepper_output"

    String prefix = basename(bam, ".bam") + ".deepvariant_pepper"

    command <<<
        # avoid the infamous pipefail 141 https://stackoverflow.com/questions/19120263
        set -eux
        SM=$(samtools view -H ~{bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g')

        set -euxo pipefail

        touch ~{bai}
        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        mkdir -p "~{output_root}"

        run_pepper_margin_deepvariant \
            call_variant \
            -b ~{bam} \
            -f ~{ref_fasta} \
            -t "${num_core}" \
            -s "${SM}" \
            -o "~{output_root}" \
            -p "~{prefix}" \
            --gvcf \
            --phased_output \
            --ont

        df -h .
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
        File VCF        = "~{output_root}/~{prefix}.vcf.gz"
        File VCF_tbi    = "~{output_root}/~{prefix}.vcf.gz.tbi"

        File gVCF       = "~{output_root}/~{prefix}.g.vcf.gz"
        File gVCF_tbi   = "~{output_root}/~{prefix}.g.vcf.gz.tbi"

        File phasedVCF  = "~{output_root}/~{prefix}.phased.vcf.gz"
        File phasedtbi  = "~{output_root}/~{prefix}.phased.vcf.gz.tbi"

        File hap_tagged_bam = "~{output_root}/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam"
        File hap_tagged_bai = "~{output_root}/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai"

        # maybe less useful
        File output_dir_structure = "~{output_root}/dir_structure.txt"
        File phaseset_bed = "~{output_root}/~{prefix}.phaseset.bed"
        File visual_report_html = "~{output_root}/~{prefix}.visual_report.html"
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
        zones:                  "~{sep=' ' zones}"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
