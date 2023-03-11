version 1.0

import "tasks/Finalize.wdl" as FF
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/SampleLevelAlignedMetrics.wdl" as worker

workflow CollectCoverageOverBed {

    meta {
        description: "Computes coverage over a Bed, for one unit (flowcell, sample)"
    }
    input {
        File bam
        File bai
        String id

        File targets_bed
        String targets_summary_name

        Int mapQ_filter

        String gcs_out_root_dir
    }

    parameter_meta {
        id:                     "Used for identifying this particular BAM; impacts the final output file name."
        targets_bed:            "The regions of interest."
        targets_summary_name:   "Short, descriptive name for this collection of targets."
    }

    String dir = sub(gcs_out_root_dir, "/$", "") + "metrics/coverage/"

    ##########
    call AM.MosDepthOverBed {
        input:
            bam = bam,
            bai = bai,
            bed = targets_bed,
            mapQ_filter = mapQ_filter
    }
    call worker.SummarizeDepthOverWholeBed {
        input:
            mosdepth_output = MosDepthOverBed.regions
    }

    call FF.FinalizeToFile as FinalizeSummary {
        input:
            outdir = dir, file = SummarizeDepthOverWholeBed.cov_summary,
            name = "~{id}.coverage_summary_over.~{targets_summary_name}.tsv"
    }

    output {
        File please_name_me_for_the_targets = FinalizeSummary.gcs_path
        Map[String, Float] feature_2_cov = SummarizeDepthOverWholeBed.feature_2_cov
    }
}
