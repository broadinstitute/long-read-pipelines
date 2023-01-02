version 1.0

import "tasks/Fusilli.wdl" as Fusilli

workflow IllumAssembleWithFusilliContinue {
    input {
        File illumina_fq1
        File illumina_fq2

        Array[File] references

        File previous_run

        Int k = 47
        String? gcs_output_file
    }

    call Fusilli.FusilliAssemble as Assemble {
        input:
            illumina_fq1 = illumina_fq1,
            illumina_fq2 = illumina_fq2,
            references = references,
            previous_run = previous_run,
            k = k
    }

    if(defined(gcs_output_file)) {
        call Fusilli.FinalizeFusilliRun {
            input:
                fusilli_output_tar = Assemble.fusilli_output_tar,
                gcs_output_file = select_first([gcs_output_file])
        }
    }

    output {
        File contigs = Assemble.contigs
    }
}
