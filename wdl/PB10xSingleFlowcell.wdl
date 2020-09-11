version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.30/wdl/tasks/Utils.wdl" as Utils
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.30/wdl/tasks/PBUtils.wdl" as PB
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.30/wdl/tasks/HiFi.wdl" as HIFI
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.30/wdl/tasks/AlignReads.wdl" as AR
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.30/wdl/tasks/AlignedMetrics.wdl" as AM
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.30/wdl/tasks/AnnotateAdapters.wdl" as AA
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.30/wdl/tasks/Figures.wdl" as FIG
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.30/wdl/tasks/Finalize.wdl" as FF

workflow PB10xSingleFlowcell {
    input {
        String gcs_input_dir
        String? sample_name
        Int fastq_shards = 50

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File ref_flat
        File dbsnp_vcf
        File dbsnp_tbi

        File metrics_locus

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = gcs_input_dir }

    scatter (subread_bam in FindBams.subread_bams) {
        call PB.GetRunInfo { input: subread_bam = subread_bam }

        String SM  = select_first([sample_name, GetRunInfo.run_info["SM"]])
        String PL  = "PACBIO"
        String PU  = GetRunInfo.run_info["PU"]
        String DT  = GetRunInfo.run_info["DT"]
        String ID  = PU
        String DS  = GetRunInfo.run_info["DS"]
        String DIR = SM + "." + ID

        String rg_subreads  = "@RG\\tID:~{ID}.subreads\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
        String rg_consensus = "@RG\\tID:~{ID}.consensus\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        call PB.ShardLongReads { input: unmapped_files = [ subread_bam ], num_reads_per_split = 2000000 }

        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads }
        }
    }
}
