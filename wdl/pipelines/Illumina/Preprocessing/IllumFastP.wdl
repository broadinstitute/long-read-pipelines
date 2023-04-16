version 1.0

import "../../../tasks/Preprocessing/FastP.wdl" as FastP

workflow IllumFastP {
    input {
        File illumina_fq1
        File illumina_fq2
    }

    call FastP.FastP as CallFastP {
        input:
            illumina_fq1=illumina_fq1,
            illumina_fq2=illumina_fq2
    }

    output {
        File processed_fq1 = CallFastP.processed_fq1
        File processed_fq2 = CallFastP.processed_fq2
        File unpaired_fq = CallFastP.unpaired_fq
        File fastp_report = CallFastP.fastp_report
        File fastp_json = CallFastP.fastp_json
    }
}