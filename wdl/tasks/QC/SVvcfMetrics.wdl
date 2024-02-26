version 1.0

task GatherCounts {
    meta {
        desciption:
        "Gather variant counts of a single-sample SV VCF. SVs are filtered down "
        note:
        "SVs are filtered to only PASS-ing variants and reported SVLEN no shorter than 50bp."
    }
    parameter_meta {
        size_lower_bound:
        "SVs with SVLEN annotation in the INFO column whose absolute value lower than this threshold will be dropped from reports; don't privde negative values"
        filter_pass_expression:
        "the expression to filter PASS-ing variants; call-specific; default is 'PASS'"

        counts_by_type:
        "number of variants for each SVTYPE (as reported by the SV caller)"
    }
    input {
        File  vcf
        File? tbi
        Int size_lower_bound
        String filter_pass_expression = "PASS"
    }
    output {
        Map[String, Int] counts_by_type = read_map("SVcounts.tsv")
    }

    command <<<
    set -euxo pipefail

        sm_cnt=$(zgrep "^#CHROM" ~{vcf} | tr '\t' '\n' | tail +10 | wc -l | awk '{print $1}')
        if [[ ${sm_cnt} -gt 1 ]]; then echo "not a single sample VCF"; exit 1; fi

        # filter to PASS only
        zgrep -v "^#" ~{vcf} \
            | awk -F '\t' -v pass_expr="~{filter_pass_expression}" 'BEGIN{OFS="\t"} {if ($7==pass_expr) print}' \
        > tmp.tsv
        # treat BND separately
        grep -vF 'SVTYPE=BND' tmp.tsv \
            | awk -F '\t' -v sz_limit="~{size_lower_bound}" \
                'BEGIN{OFS="\t"} match($8, /SVLEN=(-)?[0-9]+/) {s=substr($8, RSTART+6, RLENGTH-6); l=s+0; if (l>sz_limit || l<-sz_limit) print}' \
            | grep -Eo "SVTYPE=[A-Za-z]+" \
            | awk -F '=' '{print $2}' \
            | sort | uniq -c \
        > tmp.nonBND.counts.txt
        # now BND
        echo -e "$(grep -cF 'SVTYPE=BND' tmp.tsv)\tBND" > tmp.BND.counts.txt
        # cat and done
        cat tmp.BND.counts.txt tmp.nonBND.counts.txt \
        > "SVcounts.txt"
        # transform to TSV for Map-ing
        cat "SVcounts.txt"
        sed -E 's/^\s+//' "SVcounts.txt" \
            | awk -F ' ' 'BEGIN{OFS="\t"} {print $2, $1}' \
        > "SVcounts.tsv"
    >>>

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk " + (10 + ceil(size(vcf, "GiB"))) + " HDD"
        preemptible:    2
        maxRetries:     1
        docker:         "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
