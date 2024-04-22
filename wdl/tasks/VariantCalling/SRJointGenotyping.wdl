version 1.0

import "../../structs/Structs.wdl"
import "../Utility/Utils.wdl" as Utils

task CreateSampleNameMap {

    meta {
        description: "Creates the sample / name-map file of the GVCFs for ingest into ImportGVCFs.  NOTE: Some of this functionality is duplicated from Utils.InferSampleName.  This is intentional - we don't want to localize all these files or shard over potentially thousands of input GVCFs."
    }

    input {
        Array[File] gvcfs
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        gvcfs: {
            help: "Array of single-sample GVCF files.",
            localization_optional: true
        }
    }

    Int disk_size_gb = 20

    String outfile_name = "~{prefix}.sample_name_map.tsv"
    String size_file_gb = "~{prefix}.total_gvcf_file_size.txt"

    # Every so often we should reauthorize so `bcftools` can continue to access our data:
    Int re_auth_interval = 50

    command <<<
        set -euxo pipefail

        # Put our gvcfs into a file we can iterate over:
        gvcf_file_list=~{write_lines(gvcfs)}

        # Initialize a file for the sample names:
        [ -e ~{outfile_name} ] && rm -rf ~{outfile_name}

        # Set our access token:
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        # Create a temporary file to store file sizes in:
        size_file=$(mktemp)

        let i=1
        while read file_path ; do

            # Get our sample list from our file:
            bcftools query -l ${file_path} > sample_names.txt

            # Make sure we only have one sample name:
            [[ $(wc -l sample_names.txt | awk '{print $1}') -ne 1 ]] && echo "Incorrect number of sample names found in GVCF (there can be only one!): ${file_path}" && exit 1

            # Make sure the samplename has an actual name:
            [ $(grep -iq "unnamedsample" sample_names.txt) ] && echo "Sample name found to be unnamedsample in GVCF: ${file_path}" && exit 1

            # Add the sample name and GVCF path to the sample name file:
            echo -e "$(cat sample_names.txt)\t${file_path}" >> ~{outfile_name}

            # Add the file size to the size file:
            gsutil du -sac ${file_path} | tail -n1 | awk '{print $1}' >> ${size_file}

            let i=$i+1
            if [[ $i -gt ~{re_auth_interval} ]] ; then
                # Periodically we should update the token so we don't have problems with long file lists:
                export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
                i=0
            fi
        done < ${gvcf_file_list}

        # Now calculate the final file size in GB:
        # We include an additional GB in case we have a very small dataset:
        awk '{s += $1}END{print int(1+s/(1024*1024*1024))}' ${size_file} > ~{size_file_gb}
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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
        File sample_name_map = outfile_name
        Int total_gvcf_size_gb = read_int("~{size_file_gb}")
    }
}

task ImportGVCFs {

    input {
        File sample_name_map

        File interval_list

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String prefix

        Int batch_size = 50

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
                --genomicsdb-workspace-path ~{prefix}.genomicsDB \
                --batch-size ~{batch_size} \
                -L ~{interval_list} \
                --sample-name-map ~{sample_name_map} \
                --reader-threads 5 \
                --merge-input-intervals \
                --consolidate

        tar -cf ~{prefix}.genomicsDB.tar ~{prefix}.genomicsDB
    >>>

    output {
        File output_genomicsdb = "~{prefix}.genomicsDB.tar"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.5.0.0"
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

task GenotypeGVCFs {

    input {
        File input_gvcf_data
        File? input_gvcf_index  # Required if passing a VCF file.

        File interval_list

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String? dbsnp_vcf

        String prefix

        Boolean keep_combined_raw_annotations = false
        RuntimeAttr? runtime_attr_override
    }

    Int ref_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(ref_dict, "GB"))
    Int db_snp_size = ceil(size(dbsnp_vcf, "GB"))

    Int disk_size = 1 + 4*ceil(size(input_gvcf_data, "GB")) + ref_size + db_snp_size

    String dbsnp_vcf_arg = if defined(dbsnp_vcf) then "-D ~{dbsnp_vcf} " else ""

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
                -O ~{prefix}.vcf.gz \
                ~{dbsnp_vcf_arg} \
                -G StandardAnnotation  \
                --only-output-calls-starting-in-intervals \
                -V ${INPUT_FILE} \
                -L ~{interval_list} \
                ~{true='--keep-combined-raw-annotations' false='' keep_combined_raw_annotations} \
                --merge-input-intervals

        # Removed for now:
        # -G AS_StandardAnnotation
    >>>

    output {
        File output_vcf = "~{prefix}.vcf.gz"
        File output_vcf_index = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             26,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.5.0.0"
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
