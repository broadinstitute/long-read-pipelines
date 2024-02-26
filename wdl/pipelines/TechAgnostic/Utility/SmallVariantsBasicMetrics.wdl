version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/VariantUtils.wdl" as VU

import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/Utility/ReadLengths.wdl" as ReLU

workflow SmallVariantsBasicMetrics {
    meta {
        description:
        "Basic metrics on small variants VCF (not gVCF)"
        note:
        "Since we aren't sure which QUAL filter will be good for you to look at, we run a range of potential QUAL filter values."
    }

    parameter_meta {
        vcf:
        "Note this is not the gVCF."

        titrated_qualFilter_smallvar_stats:
        "Because we run a range of potential QUAL filter threshold values, the result is one entry per value, with the threshold value as the 'key'."
    }

    input {
        File vcf
        File tbi

        String gcs_out_root_dir
    }

    output {
        Array[Pair[Int, Map[String, String]]] titrated_qualFilter_smallvar_stats =
            zip(qual_range, GatherInfoForOneThreshold.res)
    }


    call VU.GetVCFSampleName { input: fingerprint_vcf = vcf }
    String sample_name = GetVCFSampleName.sample_name

    String files_outdir = sub(gcs_out_root_dir, "/$", "") + "/SmallVariantsMetrics/~{sample_name}"

    call GU.GetBetterRange { input: start = 20, end = 50, step = 1 }
    Array[Int] qual_range = GetBetterRange.your_range
    scatter (qual in qual_range) {

        call VU.SmallVariantsBasicStats { input: vcf = vcf, tbi = tbi, qual_threshold = qual}

        call Utils.GetMeanSDMedianMad as snp_GQ_stats { input:
            array_of_numbers = SmallVariantsBasicStats.biallelic_snp_GQ_array
        }
        call Utils.GetMeanSDMedianMad as indel_GQ_stats { input:
            array_of_numbers = SmallVariantsBasicStats.indel_GQ_array
        }

        # hacky use of readlength stats collection code
        # # call ReLU.Dyst as SnpHist { input: read_lengths_txt = SmallVariantsBasicStats.biallelic_snp_GQ_array }
        call ReLU.ReverseYield as SnpGQDecile { input:
            read_lengths_txt = SmallVariantsBasicStats.biallelic_snp_GQ_array
        }

        # # call ReLU.Dyst as IndelHist { input: read_lengths_txt = SmallVariantsBasicStats.indel_GQ_array }
        call ReLU.ReverseYield as IndelGQDecile { input:
            read_lengths_txt = SmallVariantsBasicStats.indel_GQ_array
        }

        call VU.PlotIndelSizeDist { input:
            indel_size_dist_tsv = SmallVariantsBasicStats.indel_sizes,
            prefix = basename(vcf, ".vcf.gz") + ".QUAL~{qual}-PASS"
        }

        call FF.FinalizeToFile as SaveSubstitutionMatrix { input:
            file = SmallVariantsBasicStats.st_matrix,
            outdir = files_outdir
        }
        call FF.FinalizeToFile as SaveInDelSizeArray { input:
            file = SmallVariantsBasicStats.indel_sizes,
            outdir = files_outdir
        }
        call FF.FinalizeToFile as SaveInDelSizePlot { input:
            file = PlotIndelSizeDist.plot,
            outdir = files_outdir
        }

        call GatherInfoForOneThreshold { input:
            small_variants_stats = SmallVariantsBasicStats.small_variants_stats,
            small_variants_GQ_stats = {
                "GQ_biallelic_snp_mean":    snp_GQ_stats.mean,
                "GQ_biallelic_snp_stddev":  snp_GQ_stats.stddev,
                "GQ_biallelic_snp_median":  snp_GQ_stats.median,
                "GQ_biallelic_snp_mad":     snp_GQ_stats.mad,
                "GQ_indel_mean":            indel_GQ_stats.mean,
                "GQ_indel_stddev":          indel_GQ_stats.stddev,
                "GQ_indel_median":          indel_GQ_stats.median,
                "GQ_indel_mad":             indel_GQ_stats.mad
                },
            biallelicSnp_GQ_deciles = SnpGQDecile.reverse_yield,
            indel_GQ_deciles = IndelGQDecile.reverse_yield,
            substitution_matrix_tsv = SaveSubstitutionMatrix.gcs_path,
            indel_size_dist_plot = SaveInDelSizeArray.gcs_path,
            indel_size_dist_raw = SaveInDelSizePlot.gcs_path
        }
    }
}

task GatherInfoForOneThreshold {
    meta {
        desciption:
        "A custom task for formatting outputs"
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
