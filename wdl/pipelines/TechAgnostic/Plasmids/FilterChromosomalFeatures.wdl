version 1.0

import "../../../tasks/Plasmids/MOBsuite.wdl" as MOBsuite

workflow FilterChromosomalFeatures {
    input {
        File base_gff3
        File mobsuite_contig_report
    }

    call MOBsuite.GFFFilterChromosomal as Filter {
        input:
            base_gff3=base_gff3,
            mobsuite_contig_report=mobsuite_contig_report
    }

    output {
        File plasmids_gff3 = Filter.plasmids_gff3
    }
}