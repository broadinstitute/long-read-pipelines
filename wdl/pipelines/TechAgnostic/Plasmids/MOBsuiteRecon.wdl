version 1.0

import "../../../tasks/Plasmids/MOBsuite.wdl" as MOBsuite

workflow MOBsuiteRecon {
    input {
        File assembly_fasta
        File? MOBsuite_db
    }

    call MOBsuite.MOBRecon as Recon {
        input:
            assembly_fasta=assembly_fasta,
            MOBsuite_db=MOBsuite_db
    }

    output {
        File chromosome = Recon.chromosome
        File contig_report = Recon.contig_report
        File mge_report = Recon.mge_report
        File mobtyper_results = Recon.mobtyper_results
        Array[File] plasmids = Recon.plasmids
    }
}