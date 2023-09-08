version 1.0

import "../../../tasks/Utility/ReadLengths.wdl" as ReLU

workflow SmallVariantsSmellTest {
    meta {
        description: "Basic test on input small variants VCF"
    }

    input {
        File vcf
        File tbi
    }

    call SmallVariantsBasicStats { input: vcf = vcf, tbi = tbi }
    call GetMeanSDMedianMad as snp_GQ_stats   { input: array_of_numbers = SmallVariantsBasicStats.biallelic_snp_GQ_array }
    call GetMeanSDMedianMad as indel_GQ_stats { input: array_of_numbers = SmallVariantsBasicStats.indel_GQ_array }

    # hacky use of readlength stats collection code
    call ReLU.Dyst as SnpHist { input: read_lengths_txt = SmallVariantsBasicStats.biallelic_snp_GQ_array }
    call ReLU.ReverseYield as SnpGQDecile { input: read_lengths_txt = SmallVariantsBasicStats.biallelic_snp_GQ_array }

    call ReLU.Dyst as IndelHist { input: read_lengths_txt = SmallVariantsBasicStats.indel_GQ_array }
    call ReLU.ReverseYield as IndelGQDecile{ input: read_lengths_txt = SmallVariantsBasicStats.indel_GQ_array }

    output {
        # primary output
        Map[String, Float] rename_me_small_variants_stats = SmallVariantsBasicStats.small_variants_stats
        # secondary output (will look ugly, raw, on Terra)
        Map[String, Float] rename_me_small_variants_GQ_stats = {
            "biallelic_snp_mean": snp_GQ_stats.mean,
            "biallelic_snp_stddev": snp_GQ_stats.stddev,
            "biallelic_snp_median":snp_GQ_stats.median,
            "biallelic_snp_mad":snp_GQ_stats.mad,
            "indel_mean": indel_GQ_stats.mean,
            "indel_stddev": indel_GQ_stats.stddev,
            "indel_median":indel_GQ_stats.median,
            "indel_mad":indel_GQ_stats.mad
        }
        Array[Float] rename_me_biallelicSnp_GQ_deciles = SnpGQDecile.reverse_yield
        File rename_me_biallelicSnp_GQ_histogram = SnpHist.histogram
        Array[Float] rename_me_indel_GQ_deciles = IndelGQDecile.reverse_yield
        File rename_me_indel_GQ_histogram = IndelHist.histogram
    }
}

task SmallVariantsBasicStats {
    input {
        File vcf
        File tbi
    }

    String prefix = basename(vcf, ".vcf.gz")

    command <<<
        set -euxo pipefail

        # filter
        bcftools view -f .,PASS "~{vcf}" \
        | bcftools filter -i "QUAL>=40" -O z -o "~{prefix}.QUAL40-PASS.vcf.gz"

        # stats
        stats_txt="~{prefix}.QUAL40-PASS.bcftools_stats.txt"
        bcftools stats "~{prefix}.QUAL40-PASS.vcf.gz" \
        > "${stats_txt}"

        # snp count
        grep -v "^#" "${stats_txt}" \
        | grep -F 'number of SNPs' \
        | awk -F '\t' '{print $NF}' \
        > cnt_snps.txt

        # indel count
        grep -v "^#" "${stats_txt}" \
        | grep -F 'number of indels' \
        | awk -F '\t' '{print $NF}' \
        > cnt_indels.txt

        # Ti/TV ratio
        grep -v "^#" "${stats_txt}" \
        | grep "^TSTV" | head -n1 | awk '{print $NF}' \
        | awk -F '\t' '{print $NF}' \
        > tstv.txt

        # bi-allelic HET HOM ratio
        bcftools view -m2 -M2 -v snps "~{prefix}.QUAL40-PASS.vcf.gz" \
        | grep -v "^#" \
        | awk -F '\t' '{print $NF}' | awk -F ':' '{print $1}' \
        | tr '|' '/' \
        | sort | uniq -c \
        > "GT.summary.txt"
        cat "GT.summary.txt"
        het_cnt_1=$(grep "0/1$" "GT.summary.txt" | head -n1 | awk '{print $1}')
        het_cnt_2=$(grep "1/0$" "GT.summary.txt" | head -n1 | awk '{print $1}')
        hom_cnt=$(grep "1/1$" "GT.summary.txt" | head -n1 | awk '{print $1}')
        het_cnt=$((het_cnt_1 + het_cnt_2))

        rm -f output_stats.tsv
        echo -e "cnt_snps\t$(cat cnt_snps.txt)" >> output_stats.tsv
        echo -e "cnt_indels\t$(cat cnt_indels.txt)" >> output_stats.tsv
        echo -e "titv_ratio\t$(cat tstv.txt)" >> output_stats.tsv
        echo -e "het_snp_cnt\t${het_cnt}" >> output_stats.tsv
        echo -e "hom_snp_cnt\t${hom_cnt}" >> output_stats.tsv

        # GQ array
        bcftools view \
            -m2 -M2 -v snps \
            "~{prefix}.QUAL40-PASS.vcf.gz" \
        | bcftools query -f'[%GQ]\n'  \
        > biallelic.snps.GQ.array.txt
        # indels of an individual might be triallelic site, more often than snp sites
        bcftools view \
            -m2 -M3 -v indels \
            "~{prefix}.QUAL40-PASS.vcf.gz" \
        | bcftools query -f'[%GQ]\n'  \
        > indels.GQ.array.txt

    >>>

    output {
        Map[String, Float] small_variants_stats = read_map("output_stats.tsv")
        File biallelic_snp_GQ_array = "biallelic.snps.GQ.array.txt"
        File indel_GQ_array = "indels.GQ.array.txt"
    }

    runtime {
        disks: "local-disk 20 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.1"
    }
}

task GetMeanSDMedianMad {
    input {
        File array_of_numbers
    }

    output {
        Float mean = read_float("mean.txt")
        Float stddev = read_float("stddev.txt")
        Float median = read_float("median.txt")
        Float mad = read_float("mad.txt")
    }
    command <<<
        set -euxo pipefail

        datamash -H \
            mean 1 median 1 sstdev 1 mad 1 \
        < ~{array_of_numbers} \
        | datamash transpose \
        > result.tsv

        grep "^mean"   result.tsv | awk -F '\t' '{print $2}'> "mean.txt"
        grep "^sstdev" result.tsv | awk -F '\t' '{print $2}'> "stddev.txt"
        grep "^median" result.tsv | awk -F '\t' '{print $2}'> "median.txt"
        grep "^mad"    result.tsv | awk -F '\t' '{print $2}'> "mad.txt"
    >>>

    runtime {
        disks: "local-disk 10 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
}
