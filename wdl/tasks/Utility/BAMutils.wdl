version 1.0

task GetReadGroupInfo {
    meta {
        desciption:
        "Get some read group information given a single-readgroup BAM. If the requested keys are absent, a null value is assigned in the returned entry. If the BAM contains multiple read groups, this will fail."
    }

    parameter_meta {
        uBAM: "The input BAM file."
        keys: "A list of requested fields in the RG line, e.g. ID, SM, LB."
    }

    input {
        String uBAM  # not using file as call-caching brings not much benefit

        Array[String] keys
        String null_value_representation = "None"
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{uBAM} | grep "^@RG" > one_rg_per_line.txt
        num_rgs=$(wc -l one_rg_per_line.txt | awk '{print $1}')
        if [[ "${num_rgs}" -gt 1 ]]; then echo "More than one read groups found!"  && exit 1; fi

        tr '\t' '\n' < one_rg_per_line.txt > rg_fields.txt

        for attribute in ~{sep=' ' keys}; do
            if grep -q "^${attribute}" rg_fields.txt; then
                value=$(grep "^${attribute}" rg_fields.txt | awk -F ':' '{print $2}')
            else
                value="~{null_value_representation}"
            fi
            echo -e "${attribute}\t${value}" >> "result.tsv"
        done
    >>>

    output {
        Map[String, String] read_group_info = read_map("result.tsv")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task ResetSamplename {
    meta {
        desciption: "Reset the SM entry in the input bam's readgroup lines."
    }
    parameter_meta {
        bam: {localization_optional: true}
    }
    input {
        File bam
        File? bai
        String sample_name
    }

    String prefix = basename(bam, ".bam")
    Int disk_size = 50 + 2 * ceil(size(bam, "GiB"))

    String local_bam = "/cromwell_root/~{prefix}.bam"

    command <<<
        set -euxo pipefail

        time gcloud storage cp ~{bam} ~{local_bam}

        samtools view --no-PG -H ~{local_bam} > header.txt
        grep -v "^@SQ" header.txt

        # fix SM in the RG lines
        grep "^@RG" header.txt > rg_lines.txt
        if ! grep -qF "SM:" rg_lines.txt; then
            sed -i "s/$/SM:tbd/" rg_lines.txt
        fi
        awk -v sm="~{sample_name}" -F '\t' 'BEGIN {OFS="\t"} { for (i=1; i<=NF; ++i) { if ($i ~ "SM:") $i="SM:"sm } print}' \
            rg_lines.txt \
            > fixed_rg_lines.txt
        cat fixed_rg_lines.txt

        # paste things back
        grep -v "^@RG" header.txt > otherlines.txt
        cat otherlines.txt fixed_rg_lines.txt > fixed_header.txt

        time samtools reheader fixed_header.txt ~{local_bam} > "~{prefix}.ResetSamplename.bam"
        if ~{defined(bai)}; then
            time samtools index -@1 "~{prefix}.ResetSamplename.bam"
        fi
    >>>

    output {
        File  reheadered_bam = "~{prefix}.ResetSamplename.bam"
        File? reheadered_bai = "~{prefix}.ResetSamplename.bam.bai"
    }

    runtime {
        cpu:            2
        memory:         "8 GiB"
        disks:          "local-disk ~{disk_size} SSD"
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }
}
