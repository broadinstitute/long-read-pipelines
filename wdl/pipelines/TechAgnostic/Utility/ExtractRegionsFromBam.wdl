version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/SRUtils.wdl" as SRUtils
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ExtractRegionsFromBam {
    meta {
        desciption: "Extract reads from the given bam file which overlap the regions in the given bed file."
    }

    input {
        String gcs_bam_path
        File regions_bed

        String participant_name
        String extraction_comment

        String gcs_out_root_dir
    }

    parameter_meta {
        gcs_bam_path: "GCS URL to bam file from which to extract reads."
        regions_bed: "Bed file containing regions for which to extract reads."
        participant_name:    "Participant (or sample) name for the given bam file."
        extraction_comment: "Comment to add to the end of the output filename."
        gcs_out_root_dir:    "Output folder into which to place the results of this workflow."
    }

    # First clean the extraction comment:
    String clean_comment = sub(extraction_comment, "[!$\t\n\r/<>:\"\\|?*&%#@'`~\[\]\{\}]", "_")

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ExtractRegionsFromBam/~{participant_name}_~{extraction_comment}"

    call Utils.GetReadsInBedFileRegions as GetReadsInBedFileRegions {
        input:
            gcs_bam_path = gcs_bam_path,
            regions_bed = regions_bed,
            prefix = "~{participant_name}_~{extraction_comment}",
    }

    call SRUtils.BamToFq as Bam2Fastq {
        input:
            bam = GetReadsInBedFileRegions.bam,
            prefix = "~{participant_name}_~{extraction_comment}"
    }

    call FF.FinalizeToFile as FinalizeBam { input: outdir = outdir, file = GetReadsInBedFileRegions.bam }
    call FF.FinalizeToFile as FinalizeBai { input: outdir = outdir, file = GetReadsInBedFileRegions.bai }
    call FF.FinalizeToFile as FinalizeFqEnd1 { input: outdir = outdir, file = Bam2Fastq.fq_end1 }
    call FF.FinalizeToFile as FinalizeFqEnd2 { input: outdir = outdir, file = Bam2Fastq.fq_end2 }
    call FF.FinalizeToFile as FinalizeFqUnpaired { input: outdir = outdir, file = Bam2Fastq.fq_unpaired }

    output {
        File bam = FinalizeBam.gcs_path
        File bai = FinalizeBai.gcs_path
        File fq_end1 = FinalizeFqEnd1.gcs_path
        File fq_end2 = FinalizeFqEnd2.gcs_path
        File fq_unpaired = FinalizeFqUnpaired.gcs_path

    }
}
