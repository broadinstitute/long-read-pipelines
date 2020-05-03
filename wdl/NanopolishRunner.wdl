version 1.0

import "tasks/Nanopolish.wdl" as Nanopolish
import "tasks/Finalize.wdl" as FF

workflow NanopolishRunner {
    input {
        String fast5_gcs_dir
        File combined_read_fasta
        File sequencing_summary
        File draft_assembly_fasta
        File output_dir
    }

    call Nanopolish.PolishAssembly {
        input:
            fast5_gcs_dir = fast5_gcs_dir,
            combined_read_fasta = combined_read_fasta,
            sequencing_summary =  sequencing_summary,
            draft_assembly_fasta = draft_assembly_fasta,
            output_dir = output_dir
    }

    call FF.FinalizeToDir {
        input:
            files = [PolishAssembly.polished_assembly],
            outdir = output_dir
    }
}