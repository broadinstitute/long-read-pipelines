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

    output {
        Array[Pair[String, Map[String, String]]] titrated_qualFilter_smallvar_stats = zip(qual_range, GatherInfoForOneThreshold.res)
    }

    call GetBetterRange { input: start = 10, end = 60, step = 1 }
    Array[Int] qual_range = GetBetterRange.your_range
    scatter (qual in qual_range) {
        call SmallVariantsBasicStats { input: vcf = vcf, tbi = tbi, qual_threshold = qual}

        call GetMeanSDMedianMad as snp_GQ_stats   { input: array_of_numbers = SmallVariantsBasicStats.biallelic_snp_GQ_array }
        call GetMeanSDMedianMad as indel_GQ_stats { input: array_of_numbers = SmallVariantsBasicStats.indel_GQ_array }

        # hacky use of readlength stats collection code
        # call ReLU.Dyst as SnpHist { input: read_lengths_txt = SmallVariantsBasicStats.biallelic_snp_GQ_array }
        call ReLU.ReverseYield as SnpGQDecile { input: read_lengths_txt = SmallVariantsBasicStats.biallelic_snp_GQ_array }

        # call ReLU.Dyst as IndelHist { input: read_lengths_txt = SmallVariantsBasicStats.indel_GQ_array }
        call ReLU.ReverseYield as IndelGQDecile{ input: read_lengths_txt = SmallVariantsBasicStats.indel_GQ_array }

        call PlotIndelSizeDist { input: indel_size_dist_tsv = SmallVariantsBasicStats.indel_sizes, prefix = "~{qual}"}

        call GatherInfoForOneThreshold { input:
            small_variants_stats = SmallVariantsBasicStats.small_variants_stats,
            small_variants_GQ_stats = {
                "GQ_biallelic_snp_mean": snp_GQ_stats.mean,
                "GQ_biallelic_snp_stddev": snp_GQ_stats.stddev,
                "GQ_biallelic_snp_median":snp_GQ_stats.median,
                "GQ_biallelic_snp_mad":snp_GQ_stats.mad,
                "GQ_indel_mean": indel_GQ_stats.mean,
                "GQ_indel_stddev": indel_GQ_stats.stddev,
                "GQ_indel_median":indel_GQ_stats.median,
                "GQ_indel_mad":indel_GQ_stats.mad
                },
            biallelicSnp_GQ_deciles = SnpGQDecile.reverse_yield,
            indel_GQ_deciles = IndelGQDecile.reverse_yield,
            substitution_matrix_tsv = SmallVariantsBasicStats.st_matrix,
            indel_size_dist_plot = PlotIndelSizeDist.plot,
            indel_size_dist_raw = SmallVariantsBasicStats.indel_sizes
        }
    }
}

task GetBetterRange {
    meta {
        desciption: "Becareful what you provide, as I'm not that smart."
    }

    input {
        Int start
        Int end
        Int step
    }

    output {
        Array[String] your_range = read_lines("result.txt")
    }
    #############################################################
    command <<<
        set -euxo pipefail
        seq ~{start} ~{step} ~{end} | tee > "result.txt"
    >>>

    #############################################################
    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
        preemptible: 1
        maxRetries: 1
    }
}

task SmallVariantsBasicStats {
    input {
        File vcf
        File tbi
        Int qual_threshold
    }

    String prefix = basename(vcf, ".vcf.gz")

    String filtered_vcf = "~{prefix}.QUAL~{qual_threshold}-PASS.vcf.gz"

    command <<<
        set -euxo pipefail


        # filter
        bcftools view -f .,PASS "~{vcf}" \
        | bcftools filter \
            -i "QUAL>=~{qual_threshold}" \
            -O z \
        -o "~{filtered_vcf}"

        # stats
        stats_txt="~{prefix}.QUAL~{qual_threshold}-PASS.bcftools_stats.txt"
        bcftools stats "~{filtered_vcf}" \
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
        bcftools view -m2 -M2 -v snps "~{filtered_vcf}" \
        | grep -v "^#" \
        | awk -F '\t' '{print $NF}' | awk -F ':' '{print $1}' \
        | tr '|' '/' \
        | sort | uniq -c \
        | tee \
        > "GT.summary.txt"
        if grep -q "0/1$" "GT.summary.txt"; then
            het_cnt_1=$(grep "0/1$" "GT.summary.txt" | head -n1 | awk '{print $1}')
        else
            het_cnt_1=0
        fi
        if grep -q "1/0$" "GT.summary.txt"; then
            het_cnt_2=$(grep "1/0$" "GT.summary.txt" | head -n1 | awk '{print $1}')
        else
            het_cnt_2=0
        fi
        het_cnt=$((het_cnt_1 + het_cnt_2))
        if grep "1/1$" "GT.summary.txt"; then
            hom_cnt=$(grep "1/1$" "GT.summary.txt" | head -n1 | awk '{print $1}')
        else
            hom_cnt=0
        fi

        rm -f output_stats.tsv
        echo -e "cnt_snps\t$(cat cnt_snps.txt)" >> output_stats.tsv
        echo -e "cnt_indels\t$(cat cnt_indels.txt)" >> output_stats.tsv
        echo -e "titv_ratio\t$(cat tstv.txt)" >> output_stats.tsv
        echo -e "het_snp_cnt\t${het_cnt}" >> output_stats.tsv
        echo -e "hom_snp_cnt\t${hom_cnt}" >> output_stats.tsv
        cat output_stats.tsv

        # GQ array
        bcftools view \
            -m2 -M2 -v snps \
            "~{filtered_vcf}" \
        | bcftools query -f'[%GQ]\n' \
        > biallelic.snps.GQ.array.txt
        touch biallelic.snps.GQ.array.txt
        # indels of an individual might be triallelic site, more often than snp sites
        bcftools view \
            -m2 -M3 -v indels \
            "~{filtered_vcf}" \
        | bcftools query -f'[%GQ]\n'  \
        > indels.GQ.array.txt
        touch indels.GQ.array.txt  # this might be empty

        wc -l ./*.GQ.array.txt

        grep -A13 "^# ST, Substitution types:$" "${stats_txt}" \
            | tail +3 | awk -F '\t' '{print $3"\t"$4}' \
        > substitution.matrix.tsv

        sed -n '/^# IDD, InDel distribution:$/,/^# ST, Substitution types:$/p' "${stats_txt}" \
            | tail +3 | head -n -1 | awk -F '\t' '{print $3"\t"$4}' \
        > indel.size.distribution.tsv
    >>>

    output {
        File post_filter_vcf = "~{filtered_vcf}"
        Map[String, Float] small_variants_stats = read_map("output_stats.tsv")
        File biallelic_snp_GQ_array = "biallelic.snps.GQ.array.txt"
        File indel_GQ_array = "indels.GQ.array.txt"

        File st_matrix = "substitution.matrix.tsv"
        File indel_sizes = "indel.size.distribution.tsv"
    }

    runtime {
        disks: "local-disk 20 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.1"
        preemptible: 1
        maxRetries: 1
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
        preemptible: 1
        maxRetries: 1
    }
}

task PlotIndelSizeDist {
    input {
        File indel_size_dist_tsv
        String prefix
    }
    output {
        File plot = "~{prefix}.InDel.sizes.pdf"
    }

    command <<<
        set -euxo pipefail

        python << CODE
        import pandas as pd
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_style("white")

        df = pd.read_csv("~{indel_size_dist_tsv}", header=None, sep='\t', names=['size', 'count'])

        fig, axs = plt.subplots(figsize=(16,9), dpi=300)
        _ = sns.scatterplot(data=df, x='size', y='count', s=50, ax = axs)

        axs.set_xlim(-20, 20)
        axs.set_ylim(100, 2_000_000)

        axs.set_xlabel('InDel size', fontdict={'size':20, 'weight':'bold'})
        axs.set_ylabel('Count', fontdict={'size':20, 'weight':'bold'})
        axs.tick_params(axis='both', which='major', labelsize=18)
        axs.set_yscale('log')

        plt.savefig("~{prefix}.InDel.sizes.pdf", dpi=300)
        CODE
    >>>
    runtime {
        disks: "local-disk 10 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-papermill-base:2.3.4"
        preemptible: 1
        maxRetries: 1
    }
}

task GatherInfoForOneThreshold {
    meta {
        desciption: ""
    }
    parameter_meta {

    }
    input {
        Map[String, Float] small_variants_stats
        Map[String, Float] small_variants_GQ_stats
        Array[Float] biallelicSnp_GQ_deciles
        Array[Float] indel_GQ_deciles
        String substitution_matrix_tsv
        String indel_size_dist_plot
        String indel_size_dist_raw
    }
    output {
        Map[String, String] res = read_map("result.tsv")
    }

    #############################################################
    command <<<
        set -euxo pipefail

        # first concat input map
        cat \
            ~{write_map(small_variants_stats)} \
            ~{write_map(small_variants_GQ_stats)} \
        > result.tsv

        # then pack GQ deciles
        echo -e "biallelicSnp_GQ_deciles\t~{sep=',' biallelicSnp_GQ_deciles}" >> result.tsv
        echo -e "indel_GQ_deciles\t~{sep=',' indel_GQ_deciles}" >> result.tsv

        # then point to files
        echo -e "indel_size_dist_plot\t~{indel_size_dist_plot}" >> result.tsv
        echo -e "indel_size_dist_raw\t~{indel_size_dist_raw}" >> result.tsv
        echo -e "subst_mat\t~{substitution_matrix_tsv}" >> result.tsv

        cat result.tsv
    >>>

    #############################################################
    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
        preemptible: 1
        maxRetries: 1
    }
}
