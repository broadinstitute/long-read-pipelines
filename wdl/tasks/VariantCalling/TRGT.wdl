version 1.0

import "../../structs/Structs.wdl"

task processWithTRGT {
    meta {
        description: "Uses TRGT to size TRs in a PacBio HiFi BAM."
    }

    input {
        File input_bam
        File input_bam_bai
        Boolean is_female

        String outprefix

        File ref_fasta
        File ref_fasta_index

        File repeatCatalog
        String catalog_name  # used in naming, so be linux friendly

        Boolean verbose = false

        RuntimeAttr? runtime_attr_override
    }

    output {
        File trgt_output_vcf     = "~{vcf_out_name}.vcf.gz"
        File trgt_output_vcf_idx = "~{vcf_out_name}.vcf.gz.tbi"
        # File trgt_output_bam   = "~{vcf_out_name}.spanning.bam"  # omitting for now, to save storage space
    }

    String vcf_out_name = outprefix + "_trgt." + catalog_name

    String karyotype = if(is_female) then "XX" else "XY"

    command <<<
    set -euxo pipefail

        nproc=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        time \
        trgt \
            ~{true='--verbose' false=' ' verbose} \
            genotype \
            --genome ~{ref_fasta} \
            --repeats ~{repeatCatalog} \
            --reads ~{input_bam} \
            --threads "${nproc}" \
            --output-prefix ~{vcf_out_name} \
            --karyotype ~{karyotype} \

        find . | sed -e "s/[^-][^\/]*\// |/g" -e "s/|\([^ ]\)/|-\1/"

        # sort & index
        bcftools sort \
            -Ob \
            -o ~{vcf_out_name}.sorted.vcf.gz \
            ~{vcf_out_name}.vcf.gz
        mv  ~{vcf_out_name}.sorted.vcf.gz \
            ~{vcf_out_name}.vcf.gz
        bcftools index \
            -t \
            ~{vcf_out_name}.vcf.gz
    >>>

    #########################
    Int disk_sz = 10 + 2 * ceil(size(input_bam, "GiB"))

    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             12,
        disk_gb:            disk_sz,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-trgt:3.0.0"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
   }
}
