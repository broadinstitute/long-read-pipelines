version 1.0

import "tasks/Utils.wdl" as Utils
#import "tasks/Structs.wdl"
#import "tasks/AlignReads.wdl"
import "tasks/Wtdbg.wdl" as Wtdbg
import "tasks/Quast.wdl" as Quast
import "tasks/CallAssemblyVariants.wdl" as AV

workflow LocalAssembly {
    input {
        Array[String]+ loci
        String region_size
        String aligned_bam
        File   aligned_bai
        String prefix
        String preset
#        Boolean? add_unaligned_reads = false
#        Boolean? run_quast = false

        File ref_map_file
    }

    parameter_meta {
        loci:          "Loci to assemble. At least one is required. Reads from all loci will be merged for assembly. Format: [\"chr1:1000-2000\", \"chr1:5000-10000\"]"
        region_size:   "Estimated size of region to assemble. Can use k/m/g suffixes (e.g. 3g for the human genome)"
        aligned_bam:   "aligned file"
        aligned_bai:   "index file"
        prefix:        "prefix for output files"
        preset:        "data preset (\"rs\" for PacBio RSII, \"sq\" for PacBio Sequel, \"ccs\" for PacBio CCS reads and \"ont\" for Oxford Nanopore)"

#        add_unaligned_reads: "set to true to include unaligned reads in the assembly (default: false)"
#        run_quast:           "set to true to run Quast on the assembly (default: false)"

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

    call Wtdbg.Assemble {
        input:
            reads = BamToFastq.reads_fq,
            genome_size = region_size,
            prefix = prefix,
            preset = preset
    }

    Boolean run_quast=false

    if (run_quast) {
        call Quast.Quast {
            input:
                ref = ref_map['fasta'],
                assemblies = [ Assemble.fa ]
        }
    }

    call AV.CallAssemblyVariants as CallAssemblyVariants {
        input:
            asm_fasta = Assemble.fa,
            ref_fasta = ref_map['fasta'],
            participant_name = prefix,
            prefix = prefix + ".wtdbg"
    }

    output {
        File local_bam = subset_bam
        File asm_wtdbg_fa = Assemble.fa
        File variants = CallAssemblyVariants.paftools_vcf
    }
}