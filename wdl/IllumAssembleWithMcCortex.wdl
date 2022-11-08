version 1.0

import "tasks/McCortex.wdl" as McCortex

workflow IllumAssembleWithMcCortex {
    input {
        String sample_id
        Int k

        File illumina_fq1
        File illumina_fq2

        Array[String] ref_ids
        Array[File] ref_fastas
    }

    call McCortex.McCortexBuild as McCortexBuild {
        input:
            sample_id = sample_id,
            k = k,
            illumina_fq1 = illumina_fq1,
            illumina_fq2 = illumina_fq2,
    }

    scatter(i in range(length(ref_ids))) {
        call McCortex.McCortexLinksForRef as BuildRefLinks {
            input:
                k = k,
                mccortex_graph = McCortexBuild.graph_cleaned,
                ref_id = ref_ids[i],
                ref_fasta = ref_fastas[i]
        }
    }

    call McCortex.MergePairedEndReads as MergeReads {
        input:
            illumina_fq1 = illumina_fq1,
            illumina_fq2 = illumina_fq2
    }

    call McCortex.McCortexLinksForReads as BuildReadLinks {
        input:
            sample_id = sample_id,
            k = k,

            merged_fq = MergeReads.merged,
            illumina_fq1 = MergeReads.split1,
            illumina_fq2 = MergeReads.split2,

            mccortex_graph = McCortexBuild.graph_cleaned
    }

    call McCortex.McCortexAssemble as Assemble {
        input:
            sample_id = sample_id,
            k = k,
            mccortex_graph = McCortexBuild.graph_cleaned,
            ref_links = BuildRefLinks.mccortex_links,
            sample_links = BuildReadLinks.mccortex_links
    }

    output {
        File mccortex_graph = McCortexBuild.graph_cleaned
        File contigs = Assemble.contigs
    }
}
