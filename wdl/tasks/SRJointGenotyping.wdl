version 1.0

import "Structs.wdl"

task ImportGVCFs {

    input {
        File sample_name_map

        File interval_list

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String prefix

        Int batch_size

        RuntimeAttr? runtime_attr_override
    }

    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(ref_dict, "GB"))

    Int disk_size = 1 + 4*ref_size

    command <<<
        set -euxo pipefail

        # Make sure that the output directory does not exist:
        [ -e ~{prefix} ] && rm -rf ~{prefix}

        #
        # Notes from WARP Team:
        #
        # We've seen some GenomicsDB performance regressions related to intervals, so we're going to pretend we only have a single interval
        # using the --merge-input-intervals arg
        # There's no data in between since we didn't run HaplotypeCaller over those loci so we're not wasting any compute

        # The memory setting here is very important and must be several GiB lower
        # than the total memory allocated to the VM because this tool uses
        # a significant amount of non-heap memory for native libraries.
        # Also, testing has shown that the multithreaded reader initialization
        # does not scale well beyond 5 threads, so don't increase beyond that.
        gatk --java-options "-Xms8000m -Xmx25000m" \
            GenomicsDBImport \
                --genomicsdb-workspace-path ~{prefix} \
                --batch-size ~{batch_size} \
                -L ~{interval_list} \
                --sample-name-map ~{sample_name_map} \
                --reader-threads 5 \
                --merge-input-intervals \
                --consolidate

        tar -cf ~{prefix}.tar ~{prefix}
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.2.6.1"
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

    output {
        File output_genomicsdb = "~{prefix}.tar"
    }
}

task GenotypeGVCFs {

    input {
        File input_gvcf_data
        File interval_list

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String dbsnp_vcf

        String prefix

        Boolean keep_combined_raw_annotations = false
        RuntimeAttr? runtime_attr_override
    }

    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(ref_dict, "GB"))
    Int db_snp_size = ceil(size(dbsnp_vcf, "GB"))

    Int disk_size = 1 + 4*ceil(size(input_gvcf_data, "GB")) + ref_size + db_snp_size

    parameter_meta {
        input_gvcf_data: { help: "Either a single GVCF file or a GenomicsDB Tar file." }
        interval_list: {
            localization_optional: true
        }
    }

    command <<<
        set -euxo pipefail

        # We must determine if our input variants are in a genomicsdb file or in a VCF.
        # The easiest way is to see if the input is a .tar file:

        is_genomics_db=true
        filename=$(basename -- "~{input_gvcf_data}")
        extension="${filename##*.}"
        if [[ "${extension}" != "tar" ]] ; then
            is_genomics_db=false
        fi

        if $is_genomics_db ; then
            tar -xf ~{input_gvcf_data}
            INPUT_FILE="gendb://$(basename ~{input_gvcf_data} .tar)"
        else
            INPUT_FILE=~{input_gvcf_data}
        fi

        gatk --java-options "-Xms8000m -Xmx25000m" \
            GenotypeGVCFs \
                -R ~{ref_fasta} \
                -O ~{prefix}.vcf \
                -D ~{dbsnp_vcf} \
                -G StandardAnnotation -G AS_StandardAnnotation \
                --only-output-calls-starting-in-intervals \
                -V ${INPUT_FILE} \
                -L ~{interval_list} \
                ~{true='--keep-combined-raw-annotations' false='' keep_combined_raw_annotations} \
                --merge-input-intervals
    >>>
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             26,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.2.6.1"
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

    output {
        File output_vcf = "~{prefix}.vcf"
        File output_vcf_index = "~{prefix}.vcf.tbi"
    }
}