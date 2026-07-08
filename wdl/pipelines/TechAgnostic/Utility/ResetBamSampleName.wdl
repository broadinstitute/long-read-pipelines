version 1.0

import "../../../tasks/Utility/BAMutils.wdl" as BU
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ResetBamSampleName {
    meta {
        desciption:
        "Reset the SM entry on every @RG line of a BAM's header to a single, consistent sample name. Fixes BAMs whose header erroneously carries multiple SM values. Uses 'samtools reheader' (header-only rewrite; does not touch alignment records)."
    }
    parameter_meta {
        bam:         "BAM whose read-group SM tags should be reset"
        bai:         "Index for the BAM (optional; a new index is produced when supplied)"
        sample_name: "The single sample name to write into every @RG SM field"
        gcs_out_root_dir: "GCS directory under which the reheadered BAM (and index) are finalized"
    }
    input {
        File bam
        File? bai
        String sample_name
        String gcs_out_root_dir
    }

    String workflow_name = "ResetBamSampleName"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_name}/~{sample_name}"

    call BU.ResetSamplename { input: bam = bam, bai = bai, sample_name = sample_name }

    call FF.FinalizeToFile as FinalizeBam { input: outdir = outdir, file = ResetSamplename.reheadered_bam }
    if (defined(bai)) {
        call FF.FinalizeToFile as FinalizeBai { input: outdir = outdir, file = select_first([ResetSamplename.reheadered_bai]) }
    }

    output {
        File  reheadered_bam = FinalizeBam.gcs_path
        File? reheadered_bai = FinalizeBai.gcs_path
    }
}
