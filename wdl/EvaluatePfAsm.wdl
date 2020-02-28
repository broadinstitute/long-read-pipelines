version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/AssemblyMetrics.wdl" as AsmMetrics

workflow EvaluatePfAsm {
    input {
        Array[File] asm
        Array[String] asm_name

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File gff
        File trf

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call AsmMetrics.AssemblyMetrics {
        input:
            asm = asm,
            asm_name = asm_name,

            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_dict = ref_dict,

            gff = gff,
            trf = trf,

            gcs_output_dir = outdir
    }

    ##########
    # Finalize
    ##########

#    call FF.FinalizeToDir as FinalizeMergedRuns {
#        input:
#            files = [ MergeRuns.merged_bam, MergeRuns.merged_bai, MergeRuns.merged_bri ],
#            outdir = outdir + "/" + DIR[0] + "/alignments"
#    }
}

