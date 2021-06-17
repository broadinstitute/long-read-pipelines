version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/StringTie2.wdl"
import "tasks/Finalize.wdl" as FF

workflow QuantifyIsoformsWithStringTie2 {
    input {
        File aligned_bam
        File aligned_bai
        File gtf

        String participant_name
#        String gcs_out_root_dir
    }

#    parameter_meta {
#        aligned_bam:        "GCS path to aligned BAM file"
#        aligned_bai:        "GCS path to aligned BAM file index"
#
#        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
#    }

#    Map[String, String] ref_map = read_map(ref_map_file)

#    String outdir = sub(gcs_out_root_dir, "/$", "") + "/QuantifyIsoformsWithStringTie2/~{participant_name}"

    call StringTie2.Quantify as QuantifyWithoutRetainedIntrons {
        input:
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai,
            gtf = gtf,
            keep_retained_introns = false,
            prefix = participant_name + ".without_retained_introns"
    }

    call StringTie2.Quantify as QuantifyWithRetainedIntrons {
        input:
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai,
            gtf = gtf,
            keep_retained_introns = true,
            prefix = participant_name + ".with_retained_introns"
    }

#    # Finalize
#    String dir = outdir + "/alignments"
#
#    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = dir, file = bam, name = "~{participant_name}.bam" }
#    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = dir, file = bai, name = "~{participant_name}.bam.bai" }
#    call FF.FinalizeToFile as FinalizeAlignedPbi { input: outdir = dir, file = pbi, name = "~{participant_name}.bam.pbi" }
#
#    output {
#        File aligned_bam = FinalizeAlignedBam.gcs_path
#        File aligned_bai = FinalizeAlignedBai.gcs_path
#        File aligned_pbi = FinalizeAlignedPbi.gcs_path
#
#        File? pbsv_vcf = FinalizePBSV.gcs_path
#        File? sniffles_vcf = FinalizeSniffles.gcs_path
#
#        File? dvp_phased_vcf = FinalizeDVPEPPERPhasedVcf.gcs_path
#        File? dvp_phased_tbi = FinalizeDVPEPPERPhasedTbi.gcs_path
#    }
}
