version 1.0

import "../../../tasks/Assembly/SPAdes.wdl" as SPAdes

workflow IllumPlasmidSPAdes {
    input {
        File illumina_fq1
        File illumina_fq2
    }

    call SPAdes.plasmidSPAdes as Assembly {
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
