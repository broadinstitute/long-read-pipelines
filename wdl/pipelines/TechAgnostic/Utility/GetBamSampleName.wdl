version 1.0

workflow GetBamSampleName {
    meta {
        description: "Report the SM (sample) tag(s) on a BAM's @RG header lines WITHOUT localizing the BAM: only the header bytes are streamed over GCS."
    }
    parameter_meta {
        bam: "BAM to inspect. Its header is streamed; the object itself is never downloaded."
    }

    input {
        File bam
    }

    call GetSampleName { input: bam = bam }

    output {
        String sample_name      = GetSampleName.sample_name       # unique SM value(s); newline-joined if >1
        Int    num_sample_names = GetSampleName.num_sample_names   # how many distinct SM values were found
    }
}

task GetSampleName {
    parameter_meta {
        bam: {
            localization_optional: true,
            description: "Streamed for its header only (not localized)."
        }
    }

    input {
        File bam
    }

    command <<<
        set -euxo pipefail

        # stream only the header; htslib reads it straight from GCS
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view --no-PG -H ~{bam} > header.txt

        if ! grep -q '^@RG' header.txt; then echo "No @RG line found in header!" >&2; exit 1; fi

        # collect the distinct SM values across all @RG lines
        grep '^@RG' header.txt | tr '\t' '\n' | grep '^SM:' | sed 's/^SM://' | sort -u > sample_names.txt

        wc -l < sample_names.txt > num_sample_names.txt
        cat sample_names.txt  # echo to stdout for the run log
    >>>

    output {
        String sample_name      = read_string("sample_names.txt")
        Int    num_sample_names = read_int("num_sample_names.txt")
    }

    runtime {
        cpu:            1
        memory:         "2 GiB"
        disks:          "local-disk 10 HDD"
        preemptible:    2
        maxRetries:     1
        docker:         "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}
