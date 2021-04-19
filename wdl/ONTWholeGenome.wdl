version 1.0

######################################################################################
## A workflow that performs single sample variant calling on Oxford Nanopore reads
## from one or more flow cells. The workflow merges multiple samples into a single BAM
## prior to variant calling.
######################################################################################

import "tasks/ONTUtils.wdl" as ONT
import "tasks/Utils.wdl" as Utils
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/CallVariantsONT.wdl" as VAR
import "tasks/Figures.wdl" as FIG
import "tasks/Finalize.wdl" as FF

workflow ONTWholeGenome {
    input {
        Array[File] aligned_bams
        Array[File] aligned_bais

        File ref_map_file
        String participant_name

        String gcs_out_root_dir
    }

    parameter_meta {
        aligned_bams:       "GCS path to aligned BAM files"
        aligned_bais:       "GCS path to aligned BAM file indices"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"
        participant_name:   "name of the participant from whom these samples were obtained"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTWholeGenome/~{participant_name}"

    # gather across (potential multiple) input raw BAMs
    if (length(aligned_bams) > 1) {
        call Utils.MergeBams as MergeAllReads { input: bams = aligned_bams, prefix = participant_name }
    }

    File bam = select_first([MergeAllReads.merged_bam, aligned_bams[0]])
    File bai = select_first([MergeAllReads.merged_bai, aligned_bais[0]])

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
            outfile = outdir + "/alignments/~{participant_name}.bam"
    }

    call FF.FinalizeToFile as FinalizeBai {
        input:
            file = bai,
            outfile = outdir + "/alignments/~{participant_name}.bai"
    }

#    call FF.FinalizeToFile as FinalizePBSV {
#        input:
#            file = CallVariants.pbsv_vcf,
#            outfile = outdir + "/variants/sv/" + basename(CallVariants.pbsv_vcf)
#    }
#
#    call FF.FinalizeToFile as FinalizeSniffles {
#        input:
#            file = CallVariants.sniffles_vcf,
#            outfile = outdir + "/variants/sv/" + basename(CallVariants.sniffles_vcf)
#    }

    call FF.FinalizeToFile as FinalizeDVPEPPERPhasedVcf {
        input:
            file = CallVariants.dvp_phased_vcf,
            outfile = outdir + "/variants/small/" + basename(CallVariants.dvp_phased_vcf)
    }

    call FF.FinalizeToFile as FinalizeDVPEPPERPhasedTbi {
        input:
            file = CallVariants.dvp_phased_tbi,
            outfile = outdir + "/variants/small/" + basename(CallVariants.dvp_phased_tbi)
    }

    call FF.FinalizeToFile as FinalizeDVPEPPERGVcf {
        input:
            file = CallVariants.dvp_g_vcf,
            outfile = outdir + "/variants/small/" + basename(CallVariants.dvp_g_vcf)
    }

    call FF.FinalizeToFile as FinalizeDVPEPPERGTbi {
        input:
            file = CallVariants.dvp_g_tbi,
            outfile = outdir + "/variants/small/" + basename(CallVariants.dvp_g_tbi)
    }

    call FF.FinalizeToFile as FinalizeDVPEPPERVcf {
        input:
            file = CallVariants.dvp_vcf,
            outfile = outdir + "/variants/small/" + basename(CallVariants.dvp_vcf)
    }

    call FF.FinalizeToFile as FinalizeDVPEPPERTbi {
        input:
            file = CallVariants.dvp_tbi,
            outfile = outdir + "/variants/small/" + basename(CallVariants.dvp_tbi)
    }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File merged_bam = FinalizeBam.gcs_path
        File merged_bai = FinalizeBai.gcs_path

#        File pbsv_vcf = FinalizePBSV.gcs_path
#        File sniffles_vcf = FinalizeSniffles.gcs_path

        File dvp_phased_vcf = FinalizeDVPEPPERPhasedVcf.gcs_path
        File dvp_phased_tbi = FinalizeDVPEPPERPhasedTbi.gcs_path
        File dvp_g_vcf = FinalizeDVPEPPERGVcf.gcs_path
        File dvp_g_tbi = FinalizeDVPEPPERGTbi.gcs_path
        File dvp_vcf = FinalizeDVPEPPERVcf.gcs_path
        File dvp_tbi = FinalizeDVPEPPERTbi.gcs_path
    }

}
