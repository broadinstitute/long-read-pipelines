version 1.0

workflow GetVcfSampleName {
    meta {
        description: "Extract the sample name from a single-sample VCF. Fails if the VCF does not contain exactly one sample."
    }
    parameter_meta {
        vcf: "VCF assumed to hold exactly one sample."
    }

    input {
        File vcf
    }

    call GetSampleName { input: vcf = vcf }

    output {
        String sample_name = GetSampleName.sample_name
    }
}

task GetSampleName {
    input {
        File vcf
    }

    Int disk_size = 10 + ceil(size(vcf, "GB"))

    command <<<
        set -euxo pipefail

        # the #CHROM line is the last line of the header and the sample columns are field 10 onward
        bcftools view -h ~{vcf} | tail -n1 | cut -f10- | tr '\t' '\n' | sed '/^$/d' > sample_names.txt

        n=$(wc -l < sample_names.txt)
        if [[ "${n}" -ne 1 ]]; then
            echo "Expected exactly one sample in the VCF, found ${n}: $(paste -sd, sample_names.txt)" >&2
            exit 1
        fi

        cat sample_names.txt
    >>>

    output {
        String sample_name = read_string("sample_names.txt")
    }

    runtime {
        cpu:            1
        memory:         "2 GiB"
        disks:          "local-disk ~{disk_size} HDD"
        preemptible:    2
        maxRetries:     1
        docker:         "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}
