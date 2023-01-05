version 1.0

import "tasks/Fusilli.wdl" as Fusilli

workflow IllumAssembleWithFusilliContinue {
    input {
        File illumina_fq1
        File illumina_fq2

        Array[File] references

        File previous_run

        Int k = 47
    }

    call Fusilli.FusilliAssemble as Assemble {
        input:
            illumina_fq1 = illumina_fq1,
            illumina_fq2 = illumina_fq2,
            references = references,
            previous_run = previous_run,
            k = k
    }

    output {
        File contigs = Assemble.contigs
        File fusilli_tarbal = Assemble.fusilli_output_tar
    }
}
