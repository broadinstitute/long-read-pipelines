version 1.0

import "tasks/Structs.wdl"
import "tasks/BWA.wdl" as BWA
import "tasks/AlignReads.wdl" as AlignReads
import "tasks/WGA.wdl" as WGA

# This workflow aids in building a variant calling truth set, by comparing two (high quality) assemblies
# with tools like mummer and minimap2. This workflow then generates several plots and reports for each variant
# identified which can be used to manually validate each variant.

workflow PrepareVarReview {
    input {
        String sample1
        String sample2

        File ref1
        File ref2

        File? ref1_gff

        Array[File]+ s1_illumina_reads
        File? s1_nanopore_reads
        File? s1_pacbio_reads

        Array[File]+ s2_illumina_reads
        File? s2_nanopore_reads
        File? s2_pacbio_reads
    }

    parameter_meta {
        sample1: "Name of the first sample (e.g. strain name)"
        sample2: "Name of the second sample (e.g. strain name)"

        ref1: "FASTA file containing the assembled genome of sample 1"
        ref2: "FASTA file containing the assembled genome of sample 2"

        ref1_gff: "Gene annotations for the first reference, to be displayed in variant review reports"

        s1_illumina_reads: "Illumina read data from sample 1."
        s1_nanopore_reads: "Oxford nanopore reads from sample 1."
        s1_pacbio_reads: "Pacbio reads from sample 1."

        s2_illumina_reads: "Illumina read data from sample 2."
        s2_nanopore_reads: "Oxford nanopore reads from sample 2."
        s2_pacbio_reads: "Pacbio reads from sample 2."
    }

    # Read alignments to both references
    # ----------------------------------
    #
    
    call BWA.BWAMem2Align as S1_SelfAlignmentIllumina {
        input: ref=ref1, reads=s1_illumina_reads, output_prefix="self_alignment/illumina/~{sample1}"
    }

    call BWA.BWAMem2Align as S2_SelfAlignmentIllumina {
        input: ref=ref2, reads=s2_illumina_reads, output_prefix="self_alignment/illumina/~{sample2}"
    }

    call BWA.BWAMem2Align as S2_to_S1_Illumina {
        input: ref=ref1, reads=s2_illumina_reads, output_prefix="alignments/illumina/~{sample2}_to_~{sample1}"
    }

    call BWA.BWAMem2Align as S1_to_S2_Illumina {
        input: ref=ref2, reads=s1_illumina_reads, output_prefix="alignments/illumina/~{sample1}_to_~{sample2}"
    }

    if(defined(s1_nanopore_reads)) {
        File np_reads_s1 = select_first([s1_nanopore_reads])

        # Self alignment
        call AlignReads.Minimap2 as S1_SelfAlignmentNanopore {
            input: ref_fasta=ref1, reads=[np_reads_s1], RG="''", map_preset="map-ont",
                prefix="self_alignment/oxford_nanopore/~{sample1}"
        }

        # Cross alignment
        call AlignReads.Minimap2 as S1_to_S2_Nanopore {
            input: ref_fasta=ref2, reads=[np_reads_s1], RG="''", map_preset="map-ont",
                prefix="alignments/oxford_nanopore/~{sample1}_to_~{sample2}"
        }
    }

    if(defined(s1_pacbio_reads)) {
        File pb_reads_s1 = select_first([s1_pacbio_reads])

        # Self alignment
        call AlignReads.Minimap2 as S1_SelfAlignmentPacBio {
            input: ref_fasta=ref1, reads=[pb_reads_s1], RG="''", map_preset="map-pb",
            prefix="self_alignment/pacbio/~{sample1}"
        }

        # Cross alignment
        call AlignReads.Minimap2 as S1_to_S2_PacBio {
            input: ref_fasta=ref2, reads=[pb_reads_s1], RG="''", map_preset="map-pb",
                prefix="alignments/pacbio/~{sample1}_to_~{sample2}"
        }
    }

    if(defined(s2_nanopore_reads)) {
        File np_reads_s2 = select_first([s2_nanopore_reads])

        # Self-alignment
        call AlignReads.Minimap2 as S2_SelfAlignmentNanopore {
            input: ref_fasta=ref2, reads=[np_reads_s2], RG="''", map_preset="map-ont",
                prefix="self_alignment/oxford_nanopore/~{sample2}"
        }

        # Cross alignment
        call AlignReads.Minimap2 as S2_to_S1_Nanopore {
            input: ref_fasta=ref1, reads=[np_reads_s2], RG="''", map_preset="map-ont",
                prefix="alignments/oxford_nanopore/~{sample2}_to_~{sample1}"
        }
    }

    if(defined(s2_pacbio_reads)) {
        File pb_reads_s2 = select_first([s2_pacbio_reads])

        # Self-alignment
        call AlignReads.Minimap2 as S2_SelfAlignmentPacBio {
            input: ref_fasta=ref2, reads=[pb_reads_s2], RG="''", map_preset="map-pb",
            prefix="self_alignment/pacbio/~{sample2}"
        }

        # Cross alignment
        call AlignReads.Minimap2 as S2_to_S1_PacBio {
            input: ref_fasta=ref1, reads=[pb_reads_s2], RG="''", map_preset="map-pb",
                prefix="alignments/pacbio/~{sample2}_to_~{sample1}"
        }
    }

    # Whole genome alignments with Nucmer and Minimap2
    # ------------------------------------------------
    # 
    call WGA.Minimap2 as S2_to_S1_MinimapWGA {
        input: ref1=ref1, ref2=ref2, prefix="wga/~{sample2}_to_~{sample1}"
    }

    call WGA.Nucmer as S2_to_S1_NucmerWGA {
        input: ref1=ref1, ref2=ref2, prefix="wga/~{sample2}_to_~{sample1}"
    }
}
