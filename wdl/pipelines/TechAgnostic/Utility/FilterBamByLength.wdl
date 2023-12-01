version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/BAMutils.wdl" as BU
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow FilterBamByLength {
    meta {
        desciption:
        "Filter a BAM (mapped or not) by sequence length."
    }
    parameter_meta {
        len_threshold_inclusive: "Reads longer than or equal to this length will be included."
    }
    input {
        File bam
        File? bai
        Int len_threshold_inclusive
        Boolean compute_yield
        Boolean conver_2_fq
        String disk_type
        String gcs_out_root_dir
    }

    String workflow_name = "FilterBamByLength"

    call BU.InferSampleName { input: bam = bam, bai = bai }
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_name}/~{InferSampleName.sample_name}"

    call BU.FilterBamByLen { input: bam = bam, bai = bai, len_threshold_inclusive = len_threshold_inclusive, compute_yield = compute_yield, disk_type = disk_type}
    call FF.FinalizeToFile as FinalizeBam { input: outdir = outdir, file = FilterBamByLen.fBAM }
    if (defined(bai)) {
        call FF.FinalizeToFile as FinalizeBai { input: outdir = outdir, file = select_first([FilterBamByLen.fBAI]) }
    }

    if (conver_2_fq) {
        call BU.BamToFastq { input: bam = FilterBamByLen.fBAM, prefix = basename(bam) + ".length-filter-~{len_threshold_inclusive}"}
        call FF.FinalizeToFile as FinalizeFastq { input: outdir = outdir, file = BamToFastq.reads_fq }
    }

    output {
        File  filtered_bam = FinalizeBam.gcs_path
        File? filtered_bai = FinalizeBai.gcs_path
        File? filtered_fastq = FinalizeFastq.gcs_path
        Float?    total_yield = FilterBamByLen.total_yield
        Float? filtered_yield = FilterBamByLen.filtered_yield
    }
}
