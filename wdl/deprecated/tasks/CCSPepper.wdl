version 1.0

#######################################################
# This pipeline calls small variants using DeepVariant.
#######################################################

import "../../structs/Structs.wdl"

task Pepper {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        Int threads
        Int memory
        String zones
        Boolean use_gpu = false

        RuntimeAttr? runtime_attr_override
    }

    Int bam_sz = ceil(size(bam, "GB"))
	Int disk_size = if bam_sz > 200 then 2*bam_sz else bam_sz + 200

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
        boot_disk_gb:       50,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "kishwars/pepper_deepvariant:r0.8" + if use_gpu then "-gpu" else ""
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
