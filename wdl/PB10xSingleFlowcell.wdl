version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/PBUtils.wdl" as PB
import "tasks/HiFi.wdl" as HIFI
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/AnnotateAdapters.wdl" as AA
import "tasks/Figures.wdl" as FIG
import "tasks/Finalize.wdl" as FF

workflow PB10xSingleFlowcell {
    input {
        String gcs_input_dir
        String? sample_name
        Int num_shards = 300

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

        # break one raw BAM into fixed number of shards
        File subread_pbi = sub(subread_bam, ".bam$", ".bam.pbi")
        call Utils.ShardLongReads { input: unaligned_bam = subread_bam, unaligned_pbi = subread_pbi, num_shards = num_shards }

        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads }
        }
    }
}
