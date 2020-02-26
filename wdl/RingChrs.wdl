version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/2.0-dockstore-test/wdl/tasks/ONTUtils.wdl" as ONT
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/2.0-dockstore-test/wdl/tasks/Utils.wdl" as Utils

workflow RingChrs {
    input {
        Array[Array[String]] samples
        String? sample_name

        File ref_fasta
        File ref_fasta_fai
        File ref_dict
        File ref_excl

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    scatter (flowcells in samples) {
        scatter (flowcell in flowcells) {

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


