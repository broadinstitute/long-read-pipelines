version 1.0

import "../../../tasks/Preprocessing/TrimGalore.wdl" as TrimGalore

workflow IllumTrimReads {
    input {
        File illumina_fq1
        File illumina_fq2

        Int? min_length = 31
    }

    call TrimGalore.TrimGalore as Trim {
        input:
            reads_fq1=illumina_fq1,
            reads_fq2=illumina_fq2,
            min_length=min_length
    }

    output {
        File trimmed_fq1 = Trim.trimmed_fq1
        File trimmed_fq2 = Trim.trimmed_fq2

        File trimming_report1 = Trim.trimming_report1
        File trimming_report2 = Trim.trimming_report2
    }
}
