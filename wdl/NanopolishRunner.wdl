version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.31/wdl/tasks/Nanopolish.wdl" as Nanopolish
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.31/wdl/tasks/Finalize.wdl" as FF

workflow NanopolishRunner {
    input {
        String fast5_dir
        File reads_fasta
        File sequencing_summary
        File draft_assembly_fasta
        File output_dir
        Int parallel_instances
    }

    call Nanopolish.PolishAssembly {
        input:
            fast5_dir = fast5_dir,
            reads_fasta = reads_fasta,
            sequencing_summary =  sequencing_summary,
            draft_assembly_fasta = draft_assembly_fasta,
            parallel_instances = parallel_instances
    }

    call FF.FinalizeToDir {
        input:
            files = [PolishAssembly.polished_assembly],
            outdir = output_dir
    }
}