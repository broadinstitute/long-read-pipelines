version 1.0

import "../../../tasks/QC/AlignedMetrics.wdl" as AM
import "../../../tasks/QC/SexConcordance.wdl" as SC

workflow SexCheckNaive {
    input {
        File bam
        File bai
        String expected_sex_type

        File? mosdepth_summary_txt
    }

    if (!defined(mosdepth_summary_txt)) {
        call AM.MosDepthWGS {input: bam = bam, bai = bai}
    }

    call SC.SummarizeCoverages {input: mosdepth_summary_txt = select_first([mosdepth_summary_txt, MosDepthWGS.summary_txt])}
    call SC.MakeACall {
        input:
            cov_chr1 = SummarizeCoverages.cov_chr1,
            cov_chrX = SummarizeCoverages.cov_chrX,
            cov_chrY = SummarizeCoverages.cov_chrY,
            expected_sex_type = expected_sex_type
    }

    output {
        Map[String, String] inferred_sex_info = MakeACall.inferred_sex_info
    }
}
