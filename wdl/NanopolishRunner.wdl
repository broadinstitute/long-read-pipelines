version 1.0

##########################################################################################
# A workflow that runs Nanopolish on a draft ONT assembly.
##########################################################################################

import "tasks/Nanopolish.wdl" as Nanopolish
import "tasks/Finalize.wdl" as FF
import "tasks/Utils.wdl" as Utils

workflow NanopolishRunner {
    input {
        String fast5_dir
        File reads
        File sequencing_summary
        File draft_assembly_fasta
        File output_dir
        Int parallel_instances
    }

    call Utils.ConvertReads {
        input:
            reads = reads,
            output_format = "fasta"
    }

    call Nanopolish.PolishAssembly {
        input:
            fast5_dir = fast5_dir,
            reads_fasta = ConvertReads.converted_reads,
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