version 1.0

import "tasks/ONTUtils.wdl" as ONT
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/Figures.wdl" as FIG
import "tasks/CallSVsONT.wdl" as SV
import "tasks/CallSmallVariantsONT.wdl" as SMV
import "tasks/Methylation.wdl" as Meth
import "tasks/Finalize.wdl" as FF

workflow ONTWholeGenome {
    input {
        Array[File] aligned_bams
        File ref_map_file

        String participant_name
        Int num_shards = 50

        String gcs_out_root_dir
    }

    parameter_meta {
        aligned_bams:     "aligned bam files"
        ref_map_file:     "table indicating reference sequence and auxillary file locations"

        participant_name: "name of the participant from whom these samples were obtained"
        num_shards:       "[default-valued] number of shards into which fastq files should be batched"

        gcs_out_root_dir: "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTWholeGenome/" + participant_name

    call Utils.MergeBams as MergeAllReads { input: bams = aligned_bams, prefix = participant_name }

    File bam = MergeAllReads.merged_bam
    File bai = MergeAllReads.merged_bai

    call AM.AlignedMetrics as PerFlowcellRunMetrics {
        input:
            aligned_bam    = bam,
            aligned_bai    = bai,
            ref_fasta      = ref_map['fasta'],
            ref_dict       = ref_map['dict'],
            gcs_output_dir = outdir + "/metrics/" + participant_name
    }

#    call FIG.Figures as PerSampleFigures {
#        input:
#            summary_files  = sequencing_summaries,
#
#            gcs_output_dir = outdir + "/metrics/combined/" + participant_name
#    }

    call SV.CallSVsONT as CallSVs {
        input:
            bam               = bam,
            bai               = bai,

            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
            tandem_repeat_bed = ref_map['tandem_repeat_bed'],
    }

    call SMV.CallSmallVariantsONT as CallSmallVariants {
        input:
            bam               = bam,
            bai               = bai,

            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
            ref_dict          = ref_map['dict'],
    }

    ##########
    # Finalize
    ##########

#    call FF.FinalizeToDir as FinalizeSVs {
#        input:
##            files = [ CallSVs.pbsv_vcf, CallSVs.sniffles_vcf, CallSVs.svim_vcf, CallSVs.cutesv_vcf ],
#            files = [ CallSVs.sniffles_vcf ],
#            outdir = outdir + "/variants"
#    }
#
#    call FF.FinalizeToDir as FinalizeSmallVariants {
#        input:
#            files = [ CallSmallVariants.longshot_vcf, CallSmallVariants.longshot_tbi ],
##                      CallSmallVariants.clair_vcf, CallSmallVariants.clair_tbi ],
#            outdir = outdir + "/variants"
#    }
#
#    call FF.FinalizeToDir as FinalizeMergedRuns {
#        input:
#            files = [ bam, bai ],
#            outdir = outdir + "/alignments"
#    }
}
