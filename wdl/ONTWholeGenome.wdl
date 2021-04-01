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

#    call FF.FinalizeToFile as FinalizePBSV {
#        input:
#            file = CallSVs.pbsv_vcf,
#            outfile = outdir + "/variants/" + basename(CallSVs.pbsv_vcf)
#    }

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

    output {
        File aligned_bam = FinalizeBam.gcs_path
        File aligned_bai = FinalizeBai.gcs_path
        File sniffles_vcf = FinalizeSniffles.gcs_path
        File svim_vcf = FinalizeSVIM.gcs_path
        File longshot_vcf = FinalizeLongshot.gcs_path
    }
}
