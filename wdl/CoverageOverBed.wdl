version 1.0

import "tasks/Finalize.wdl" as FF
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/SampleLevelAlignedMetrics.wdl" as worker

workflow CoverageOverBed {
    input {
        File bam
        File bai
        String participant_name
        File bed_to_compute_coverage

        String workflow_root_name
        String gcs_out_root_dir
    }

    call AM.MosDepthOverBed {
        input:
            bam = bam,
            bai = bai,
            bed = bed_to_compute_coverage
    }
    call worker.SummarizeDepthOverWholeBed as cov_over_region {
        input:
            mosdepth_output = MosDepthOverBed.regions
    }

    String dir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_root_name}/~{participant_name}/alignments"

    call FF.FinalizeToFile as FinalizeRegionalCoverage { input: outdir = dir, file = cov_over_region.cov_summary }

    output {
        File bed_cov_summary = FinalizeRegionalCoverage.gcs_path
    }
}