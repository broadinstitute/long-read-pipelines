version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/2.0-dockstore-test/wdl/tasks/Utils.wdl" as Utils
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/2.0-dockstore-test/wdl/tasks/AssemblyMetrics.wdl" as AsmMetrics

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

    scatter (i in range(length(asm))) {
        call AsmMetrics.AssemblyMetrics {
            input:
                asm = asm[i],
                asm_name = asm_name[i],

                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                ref_dict = ref_dict,

                gff = gff,
                trf = trf,

                gcs_output_dir = outdir + "/" + asm_name[i]
        }
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

