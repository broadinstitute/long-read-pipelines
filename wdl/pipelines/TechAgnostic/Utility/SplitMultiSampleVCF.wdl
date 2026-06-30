version 1.0

import "../../../structs/Structs.wdl"

workflow SplitMultiSampleVCF {
    meta {
        description: "Split a multi-sample VCF into individual compressed VCF files, one per sample, with corresponding index files."
        author: "Jonn Smith"
    }

    parameter_meta {
        input_vcf: "Multi-sample VCF file (can be compressed or uncompressed)"
        input_vcf_index: "Index file for the input VCF (required if VCF is compressed)"
        num_samples: "Number of samples in the input VCF (optional; default: 100)"
        sample_names: "Optional list of sample names to extract. If provided, every name must occur in input_vcf; the workflow fails (listing all absent names) when any is missing. Only the listed samples are emitted; otherwise every sample is emitted."
    }

    input {
        File input_vcf
        File? input_vcf_index

        Int num_samples = 100

        Array[String]? sample_names
    }

    call SplitMultiSampleVCFTask {
        input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            num_samples = num_samples,
            sample_names = sample_names
    }

    output {
        Array[File] sample_vcfs = SplitMultiSampleVCFTask.output_vcfs
        Array[File] sample_vcf_indices = SplitMultiSampleVCFTask.output_vcf_indices
    }
}

task SplitMultiSampleVCFTask {
    meta {
        description: "Split a multi-sample VCF into individual compressed VCF files, one per sample, with corresponding index files"
    }

    parameter_meta {
        input_vcf: "Multi-sample VCF file (can be compressed or uncompressed)"
        input_vcf_index: "Index file for the input VCF (required if VCF is compressed)"
        num_samples: "Number of samples in the input VCF (optional; default: 100)"
        sample_names: "Optional list of sample names to extract. If provided, every name must occur in input_vcf; the task fails (listing all absent names) when any is missing. Only the listed samples are emitted; otherwise every sample is emitted."
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File input_vcf
        File? input_vcf_index

        Int num_samples = 100

        Array[String]? sample_names

        RuntimeAttr? runtime_attr_override
    }

    Array[String] requested_samples = select_first([sample_names, []])

    Int disk_size = 10 + num_samples*ceil(size([input_vcf, input_vcf_index], "GB"))

    command <<<
        set -euxo pipefail

        mkdir -p out_dir

        REQUESTED_SAMPLES="~{write_lines(requested_samples)}"

        if [ -s "${REQUESTED_SAMPLES}" ]; then
            # Sample names present in the input VCF.
            bcftools query -l ~{input_vcf} > vcf_samples.txt

            # Collect every requested sample that is absent from the VCF.
            missing=()
            while IFS= read -r sample; do
                [ -z "${sample}" ] && continue
                if ! grep -Fxq -- "${sample}" vcf_samples.txt; then
                    missing+=("${sample}")
                fi
            done < "${REQUESTED_SAMPLES}"

            if [ "${#missing[@]}" -gt 0 ]; then
                echo "ERROR: the following requested sample(s) are not present in the input VCF:"
                for sample in "${missing[@]}"; do
                    echo "  ${sample}"
                done
                exit 1
            fi

            # All requested samples are present: emit only those.
            bcftools +split ~{input_vcf} -S "${REQUESTED_SAMPLES}" -Oz2 -W=tbi -o out_dir
        else
            # No sample subset requested: emit every sample.
            bcftools +split ~{input_vcf} -Oz2 -W=tbi -o out_dir
        fi

    >>>

    output {
        Array[File] output_vcfs = glob("out_dir/*.vcf.gz")
        Array[File] output_vcf_indices = glob("out_dir/*.vcf.gz.tbi")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.3"
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

