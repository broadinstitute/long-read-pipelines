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
        disk_space_multiplier: "Multiplier applied to the input VCF size when sizing the task's disk (optional; default: 4)"
        sample_names: "Optional list of sample names to extract. If provided, every name must occur in input_vcf; the workflow fails (listing all absent names) when any is missing. Only the listed samples are emitted; otherwise every sample is emitted. Mutually exclusive with sample_name_list."
        sample_name_list: "Optional file containing sample names to extract, one per line. Same semantics as sample_names. Mutually exclusive with sample_names."
    }

    input {
        File input_vcf
        File? input_vcf_index

        Int disk_space_multiplier = 4

        Array[String]? sample_names
        File? sample_name_list
    }

    call SplitMultiSampleVCFTask {
        input:
            input_vcf = input_vcf,
            input_vcf_index = input_vcf_index,
            disk_space_multiplier = disk_space_multiplier,
            sample_names = sample_names,
            sample_name_list = sample_name_list
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
        disk_space_multiplier: "Multiplier applied to the input VCF size when sizing the task's disk (optional; default: 4)"
        sample_names: "Optional list of sample names to extract. If provided, every name must occur in input_vcf; the task fails (listing all absent names) when any is missing. Only the listed samples are emitted; otherwise every sample is emitted. Mutually exclusive with sample_name_list."
        sample_name_list: "Optional file containing sample names to extract, one per line. Same semantics as sample_names. Mutually exclusive with sample_names."
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File input_vcf
        File? input_vcf_index

        Int disk_space_multiplier = 4

        Array[String]? sample_names
        File? sample_name_list

        RuntimeAttr? runtime_attr_override
    }

    Array[String] requested_samples = select_first([sample_names, []])
    Boolean have_inline_samples = length(requested_samples) > 0

    Int disk_size = 10 + disk_space_multiplier*ceil(size([input_vcf, input_vcf_index], "GB"))

    command <<<
        set -euxo pipefail

        mkdir -p out_dir

        REQUESTED_SAMPLES="~{write_lines(requested_samples)}"
        SAMPLE_NAME_LIST="~{default='' sample_name_list}"

        # Whether an inline sample_names list was given. Gate on the WDL-known
        # count rather than the size of write_lines() output: write_lines([]) is
        # not guaranteed to be a 0-byte file, so `[ -s ... ]` would spuriously
        # take the subset path for an empty request.
        HAVE_INLINE_SAMPLES=~{true="1" false="0" have_inline_samples}

        # sample_names and sample_name_list are mutually exclusive: reject both.
        if [ "${HAVE_INLINE_SAMPLES}" -eq 1 ] && [ -n "${SAMPLE_NAME_LIST}" ]; then
            echo "ERROR: 'sample_names' and 'sample_name_list' are mutually exclusive; provide at most one." >&2
            exit 1
        fi

        # Resolve the effective list of requested samples (empty => emit all).
        SAMPLES_FILE=""
        if [ "${HAVE_INLINE_SAMPLES}" -eq 1 ]; then
            SAMPLES_FILE="${REQUESTED_SAMPLES}"
        elif [ -n "${SAMPLE_NAME_LIST}" ]; then
            SAMPLES_FILE="${SAMPLE_NAME_LIST}"
        fi

        if [ -n "${SAMPLES_FILE}" ]; then
            # Requested samples absent from the VCF = set difference (requested \ present).
            # comm -23 emits lines unique to the first (requested) sorted list.
            bcftools query -l ~{input_vcf} | sort -u > vcf_samples.sorted
            { grep -v '^[[:space:]]*$' "${SAMPLES_FILE}" || true; } | sort -u > requested_samples.sorted
            comm -23 requested_samples.sorted vcf_samples.sorted > missing_samples.txt

            if [ -s missing_samples.txt ]; then
                echo "ERROR: the following requested sample(s) are not present in the input VCF:" >&2
                sed 's/^/  /' missing_samples.txt >&2
                exit 1
            fi

            # All requested samples are present: emit only those.
            bcftools +split ~{input_vcf} -S "${SAMPLES_FILE}" -Oz2 -W=tbi -o out_dir
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

