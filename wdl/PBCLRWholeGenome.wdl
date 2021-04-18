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

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCLRWholeGenome/~{SM}"

    # gather across (potential multiple) input raw BAMs
    call Utils.MergeBams as MergeAllReads { input: bams = aligned_bams, prefix = participant_name }

    File bam = MergeAllReads.merged_bam
    File bai = MergeAllReads.merged_bai

    # compute alignment metrics
    call AM.AlignedMetrics as PerSampleMetrics {
        input:
            aligned_bam    = bam,
            aligned_bai    = bai,
            ref_fasta      = ref_map['fasta'],
            ref_dict       = ref_map['dict'],
            gcs_output_dir = outdir + "/metrics/" + participant_name
    }

    call VAR.CallVariants {
        input:
            bam               = bam,
            bai               = bai,

            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
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

    call FF.FinalizeToFile as FinalizePBSV {
        input:
            file = CallSVs.pbsv_vcf,
            outfile = outdir + "/variants/" + basename(CallSVs.pbsv_vcf)
    }

    call FF.FinalizeToFile as FinalizeSniffles {
        input:
            file = CallSVs.sniffles_vcf,
            outfile = outdir + "/variants/" + basename(CallSVs.sniffles_vcf)
    }

    call FF.FinalizeToFile as FinalizeSVIM {
        input:
            file = CallSVs.svim_vcf,
            outfile = outdir + "/variants/" + basename(CallSVs.svim_vcf)
    }

    call FF.FinalizeToFile as FinalizeLongshot {
        input:
            file = CallSmallVariants.longshot_vcf,
            outfile = outdir + "/variants/" + basename(CallSmallVariants.longshot_vcf)
    }

    ##########
    # store the results into designated bucket
    ##########

    call FF.FinalizeToDir as FinalizeSVs {
        input:
            files = [ CallSVs.pbsv_vcf, CallSVs.sniffles_vcf, CallSVs.svim_vcf, CallSVs.cutesv_vcf ],
            outdir = outdir + "/variants"
    }

    call FF.FinalizeToDir as FinalizeSmallVariants {
        input:
            files = [ CallSmallVariants.longshot_vcf, CallSmallVariants.longshot_tbi ],
            outdir = outdir + "/variants"
    }
}
