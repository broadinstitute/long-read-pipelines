version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/Hifiasm.wdl" as Hifiasm
import "tasks/PBUtils.wdl" as PB

workflow LocalAssembly {
    input {
        Array[String]+ loci
        String aligned_bam
        File   aligned_bai
        String prefix
#        Boolean add_unaligned_reads = false
        File mat_yak
        File pat_yak

        File ref_map_file
    }

    parameter_meta {
        loci:          "Loci to assemble. At least one is required. Reads from all loci will be merged for assembly. Format: [\"chr1:1000-2000\", \"chr1:5000-10000\"]"
        aligned_bam:   "aligned file"
        aligned_bai:   "index file"
        prefix:        "prefix for output files"

#        add_unaligned_reads: "set to true to include unaligned reads in the assembly (default: false)"

        ref_map_file:  "table indicating reference sequence and auxillary file locations"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    scatter (locus in loci) {
        call Utils.SubsetBam {
            input:
                bam = aligned_bam,
                bai = aligned_bai,
                locus = locus
        }
    }

    if (length(loci) > 1)  {
        call Utils.MergeBams {
            input:
                bams = SubsetBam.subset_bam,
                prefix = "merged"
        }
    }

    File subset_bam = select_first([MergeBams.merged_bam, SubsetBam.subset_bam[0]])

    call Utils.BamToFastq {
        input:
            bam = subset_bam,
            prefix = prefix
    }

    call Hifiasm.Assemble_trio {
        input:
            reads = BamToFastq.reads_fq,
            prefix = prefix,
            mat_yak = mat_yak,
            pat_yak = pat_yak
    }

    call PB.Align as align_hap1 {
        input:
            bam = Assemble_trio.h1_fa,
            ref_fasta = ref_map['fasta'],
            sample_name = "~{prefix}_h1",
            map_preset = "asm5",
            drop_per_base_N_pulse_tags = false,
            prefix = "~{prefix}_h1"
    }

        call PB.Align as align_hap2 {
        input:
            bam = Assemble_trio.h2_fa,
            ref_fasta = ref_map['fasta'],
            sample_name = "~{prefix}_h2",
            map_preset = "asm5",
            drop_per_base_N_pulse_tags = false,
            prefix = "~{prefix}_h2"
    }

    output {
        File h1_fa = Assemble_trio.h1_fa
        File h2_fa = Assemble_trio.h2_fa
        File h1_aligned_bam = align_hap1.aligned_bam
        File h1_aligned_bai = align_hap1.aligned_bai
        File h2_aligned_bam = align_hap2.aligned_bam
        File h2_aligned_bai = align_hap2.aligned_bai
    }
}