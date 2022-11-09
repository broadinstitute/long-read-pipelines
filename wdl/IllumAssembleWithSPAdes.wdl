version 1.0

import "tasks/SPAdes.wdl" as SPAdes

workflow IllumAssembleWithSPAdes {
    input {
        File illumina_fq1
        File illumina_fq2
    }

    call SPAdes.SPAdesAssemble as Assembly {
        input:
            illumina_fq1=illumina_fq1,
            illumina_fq2=illumina_fq2
    }

    output {
        File scaffolds = Assembly.scaffolds
        File contigs = Assembly.contigs
        File assembly_graph = Assembly.assembly_graph
    }
}
