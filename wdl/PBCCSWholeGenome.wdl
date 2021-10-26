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
        Boolean fast_less_sensitive
    }

    parameter_meta {
        aligned_bams:       "GCS path to aligned BAM files"
        aligned_bais:       "GCS path to aligned BAM file indices"
        ref_map_file:       "table indicating reference sequence and auxillary file locations"
        participant_name:   "name of the participant from whom these samples were obtained"
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
        fast_less_sensitive:"for SV calling specifically, true indicates fast/suboptimal processing, false indicates slower but more sensitive processing"
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
                prefix = participant_name,
                fast_less_sensitive = fast_less_sensitive
        }

        String svdir = outdir + "/variants/sv"
        String smalldir = outdir + "/variants/small"

        call FF.FinalizeToFile as FinalizePBSV { input: outdir = svdir, file = CallVariants.pbsv_vcf }
        call FF.FinalizeToFile as FinalizeSniffles { input: outdir = svdir, file = CallVariants.sniffles_vcf }
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
    }
}
