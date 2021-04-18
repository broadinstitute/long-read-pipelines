version 1.0

##########################################################################################
## A workflow that performs CCS correction and variant calling on PacBio HiFi reads from a
## single flow cell. The workflow shards the subreads into clusters and performs CCS in
## parallel on each cluster.  Error-corrected reads are then variant-called.  A number of
## metrics and figures are produced along the way.
##########################################################################################

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/CallVariantsPBCLR.wdl" as VAR
import "tasks/Figures.wdl" as FIG
import "tasks/Finalize.wdl" as FF

workflow PBCLRWholeGenome {
    input {
        Array[File] aligned_bams

        File ref_map_file

        String participant_name
        Int num_shards = 50

        String gcs_out_root_dir
    }

    parameter_meta {
        bams:               "GCS path to raw subreads or CCS data"
        ref_map_file:       "table indicating reference sequence and auxillary file locations"

        participant_name:   "name of the participant from whom these samples were obtained"
        num_shards:         "[default-valued] number of sharded BAMs to create (tune for performance)"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCLRWholeGenome/~{participant_name}"

    # gather across (potential multiple) input raw BAMs
    call Utils.MergeBams as MergeAllReads { input: bams = aligned_bams, prefix = participant_name }
    call PB.PBIndex as IndexMergedReads { input: bam = MergeAllReads.merged_bam }

    File bam = MergeAllReads.merged_bam
    File bai = MergeAllReads.merged_bai
    File pbi = IndexMergedReads.pbi

    call VAR.CallVariants {
        input:
            bam               = bam,
            bai               = bai,

            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
            ref_dict          = ref_map['dict'],
            tandem_repeat_bed = ref_map['tandem_repeat_bed'],
    }

    ##########
    # Finalize
    ##########

    call FF.FinalizeToFile as FinalizeBam {
        input:
            file = bam,
            outfile = outdir + "/alignments/" + basename(bam)
    }

    call FF.FinalizeToFile as FinalizeBai {
        input:
            file = bai,
            outfile = outdir + "/alignments/" + basename(bai)
    }

    call FF.FinalizeToFile as FinalizePbi {
        input:
            file = pbi,
            outfile = outdir + "/alignments/" + basename(pbi)
    }

    call FF.FinalizeToFile as FinalizePBSV {
        input:
            file = CallVariants.pbsv_vcf,
            outfile = outdir + "/variants/sv/" + basename(CallVariants.pbsv_vcf)
    }

    call FF.FinalizeToFile as FinalizeSniffles {
        input:
            file = CallVariants.sniffles_vcf,
            outfile = outdir + "/variants/sv/" + basename(CallVariants.sniffles_vcf)
    }

    call FF.FinalizeToFile as FinalizeSVIM {
        input:
            file = CallVariants.svim_vcf,
            outfile = outdir + "/variants/sv/" + basename(CallVariants.svim_vcf)
    }

    call FF.FinalizeToFile as FinalizeCuteSV {
        input:
            file = CallVariants.cutesv_vcf,
            outfile = outdir + "/variants/sv/" + basename(CallVariants.cutesv_vcf)
    }

    call FF.FinalizeToFile as FinalizeLongshot {
        input:
            file = CallVariants.longshot_vcf,
            outfile = outdir + "/variants/small/" + basename(CallVariants.longshot_vcf)
    }

    call FF.FinalizeToFile as FinalizeLongshotTbi {
        input:
            file = CallVariants.longshot_tbi,
            outfile = outdir + "/variants/small/" + basename(CallVariants.longshot_tbi)
    }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File merged_bam = FinalizeBam.gcs_path
        File merged_bai = FinalizeBai.gcs_path
        File merged_pbi = FinalizePbi.gcs_path

        File pbsv_vcf = FinalizePBSV.gcs_path
        File sniffles_vcf = FinalizeSniffles.gcs_path
        File svim_vcf = FinalizeSVIM.gcs_path
        File cutesv_vcf = FinalizeCuteSV.gcs_path

        File longshot_vcf = FinalizeLongshot.gcs_path
        File longshot_tbi = FinalizeLongshotTbi.gcs_path
    }
}
