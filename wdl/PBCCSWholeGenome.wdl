version 1.0

######################################################################################
## A workflow that performs single sample variant calling on PacBio HiFi reads from
## one or more flow cells. The workflow merges multiple samples into a single BAM
## prior to variant calling.
######################################################################################

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/CallVariantsPBCCS.wdl" as VAR
import "tasks/Finalize.wdl" as FF

workflow PBCCSWholeGenome {
    input {
        Array[File] aligned_bams
        Array[File] aligned_bais

        File ref_map_file
        String participant_name

        Boolean call_variants = true

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

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCCSWholeGenome/~{participant_name}"

    # gather across (potential multiple) input CCS BAMs
    if (length(aligned_bams) > 1) {
        call Utils.MergeBams as MergeAllReads { input: bams = aligned_bams, prefix = participant_name }
    }

    File bam = select_first([MergeAllReads.merged_bam, aligned_bams[0]])
    File bai = select_first([MergeAllReads.merged_bai, aligned_bais[0]])

    call PB.PBIndex as IndexCCSUnalignedReads { input: bam = bam }
    File pbi = IndexCCSUnalignedReads.pbi

    if (call_variants) {
        call VAR.CallVariants {
            input:
                bam               = bam,
                bai               = bai,

                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                tandem_repeat_bed = ref_map['tandem_repeat_bed'],
        }

        String svdir = outdir + "/variants/sv"
        String smalldir = outdir + "/variants/small"

        call FF.FinalizeToFile as FinalizePBSV { input: outdir = svdir, file = CallVariants.pbsv_vcf }
        call FF.FinalizeToFile as FinalizeSniffles { input: outdir = svdir, file = CallVariants.sniffles_vcf }

        call FF.FinalizeToFile as FinalizeDVPEPPERPhasedVcf { input: outdir = smalldir, file = CallVariants.dvp_phased_vcf }
        call FF.FinalizeToFile as FinalizeDVPEPPERPhasedTbi { input: outdir = smalldir, file = CallVariants.dvp_phased_tbi }
        call FF.FinalizeToFile as FinalizeDVPEPPERGVcf { input: outdir = smalldir, file = CallVariants.dvp_g_vcf }
        call FF.FinalizeToFile as FinalizeDVPEPPERGTbi { input: outdir = smalldir, file = CallVariants.dvp_g_tbi }
        call FF.FinalizeToFile as FinalizeDVPEPPERVcf { input: outdir = smalldir, file = CallVariants.dvp_vcf }
        call FF.FinalizeToFile as FinalizeDVPEPPERTbi { input: outdir = smalldir, file = CallVariants.dvp_tbi }
    }

    # Finalize
    String dir = outdir + "/alignments"

    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = dir, file = bam, name = "~{participant_name}.bam" }
    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = dir, file = bai, name = "~{participant_name}.bam.bai" }
    call FF.FinalizeToFile as FinalizeAlignedPbi { input: outdir = dir, file = pbi, name = "~{participant_name}.bam.pbi" }

    output {
        File aligned_bam = FinalizeAlignedBam.gcs_path
        File aligned_bai = FinalizeAlignedBai.gcs_path
        File aligned_pbi = FinalizeAlignedPbi.gcs_path

        File? pbsv_vcf = FinalizePBSV.gcs_path
        File? sniffles_vcf = FinalizeSniffles.gcs_path

        File? dvp_phased_vcf = FinalizeDVPEPPERPhasedVcf.gcs_path
        File? dvp_phased_tbi = FinalizeDVPEPPERPhasedTbi.gcs_path
    }
}
