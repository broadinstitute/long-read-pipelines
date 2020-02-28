version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/2.0-dockstore-test-2/wdl/tasks/PBUtils.wdl" as PB
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/2.0-dockstore-test-2/wdl/tasks/Utils.wdl" as Utils
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/2.0-dockstore-test-2/wdl/tasks/AlignReads.wdl" as AR
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/2.0-dockstore-test-2/wdl/tasks/Finalize.wdl" as FF

workflow PBCCSOnlySingleFlowcell {
    input {
        String gcs_input_dir

        String? sample_name
        Int num_reads_per_split = 100000

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = gcs_input_dir }

    scatter (subread_bam in FindBams.subread_bams) {
        call PB.GetRunInfo { input: subread_bam = subread_bam }

        String SM  = select_first([sample_name, GetRunInfo.run_info["SM"]])
        String PU  = GetRunInfo.run_info["PU"]
        String ID  = PU
        String DIR = SM + "." + ID

        call Utils.ShardLongReads { input: unmapped_files = [ subread_bam ], num_reads_per_split = num_reads_per_split }

        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads }
            #call PB.CCSWithClasses { input: subreads = subreads }
        }

        call AR.MergeBams as MergeChunks { input: bams = CCS.consensus }

        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report }
        #call PB.MergeCCSReports as MergeCCSWithClassesReports { input: reports = CCSWithClasses.report }
        #call PB.MergeCCSClasses as MergeCCSWithClassesClasses { input: classes = CCSWithClasses.classes }
    }

    call AR.MergeBams as MergeRuns { input: bams = MergeChunks.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }

    call PB.MergeCCSReports as MergeAllCCSReports { input: reports = MergeCCSReports.report }
    #call PB.MergeCCSReports as MergeAllCCSWithClassesReports { input: reports = MergeCCSWithClassesReports.report }
    #call PB.MergeCCSClasses as MergeAllCCSWithClassesClasses { input: classes = MergeCCSWithClassesClasses.classes }

    ##########
    # Finalize
    ##########

    call FF.FinalizeToDir as FinalizeMergedRuns {
        input:
            files = [ MergeRuns.merged_bam ],
            outdir = outdir + "/" + DIR[0] + "/alignments"
    }

    call FF.FinalizeToDir as FinalizeCCSMetrics {
        input:
            files = [ MergeAllCCSReports.report ],
            #files = [ MergeAllCCSWithClassesReports.report, MergeAllCCSWithClassesClasses.classes ],
            outdir = outdir + "/" + DIR[0] + "/metrics/ccs_metrics"
    }
}
