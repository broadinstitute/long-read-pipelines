version 1.0

import "../../../tasks/Assembly/SPAdes.wdl" as SPAdes

workflow IllumErrorCorrectBayesHammer {
    input {
        File illumina_fq1
        File illumina_fq2
    }

    call SPAdes.BayesHammer as BayesHammer {
        input:
            illumina_fq1=illumina_fq1,
            illumina_fq2=illumina_fq2
    }

    output {
        File corrected_fq1 = BayesHammer.corrected_fqs[0]
        File corrected_fq2 = BayesHammer.corrected_fqs[1]
    }
}
