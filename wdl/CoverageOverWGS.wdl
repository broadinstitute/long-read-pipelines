version 1.0

import "tasks/AlignedMetrics.wdl"

workflow CoverageOverWGS {
    input {
        File bam
        File bai
    }

    call AlignedMetrics.MosDepthWGS { input: bam = bam, bai = bai}
    output {
        Float wgs_cov = MosDepthWGS.wgs_cov
    }
}
