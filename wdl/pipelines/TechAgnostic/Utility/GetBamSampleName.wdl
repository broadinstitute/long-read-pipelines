version 1.0

workflow GetBamSampleName {
    meta {
        description: "Report the SM (sample) tag(s) on a BAM's @RG header lines WITHOUT localizing the BAM: only the header bytes are streamed over GCS."
    }
    parameter_meta {
        bam: "BAM to inspect. Its header is streamed; the object itself is never downloaded."
        gcs_requester_pays_project: "GCP project to bill when the BAM lives in a requester-pays bucket. Leave unset for normal buckets."
    }

    input {
        File bam
        String? gcs_requester_pays_project
    }

    call GetSampleName { input: bam = bam, gcs_requester_pays_project = gcs_requester_pays_project }

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
        gcs_requester_pays_project: "GCP project to bill for requester-pays buckets. When set, the header is streamed via 'gcloud storage cat --billing-project' (htslib's gs:// backend cannot pass a billing project). When unset, htslib streams the header directly."
    }

    input {
        File bam
        String? gcs_requester_pays_project
    }

    String rp_project = select_first([gcs_requester_pays_project, ""])

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        # stream only the header (samtools closes the pipe after the header, so the BAM is not fully transferred)
        if [[ -n "~{rp_project}" ]]; then
            # requester-pays: gcloud can pass a billing project; htslib's gs:// backend cannot
            set +o pipefail
            gcloud storage cat --billing-project="~{rp_project}" ~{bam} | samtools view --no-PG -H - > header.txt
            set -o pipefail
        else
            samtools view --no-PG -H ~{bam} > header.txt
        fi

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
        docker:         "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}
