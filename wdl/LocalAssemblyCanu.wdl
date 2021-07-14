version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/Canu.wdl" as Canu
import "tasks/Quast.wdl" as Quast
import "tasks/CallAssemblyVariants.wdl" as AV

workflow LocalAssembly {
    input {
        Array[String]+ loci
        String region_size
        String aligned_bam
        File   aligned_bai
        String prefix
        String experiment_type
#        Boolean add_unaligned_reads = false
        Boolean run_quast = false

        File ref_map_file
    }

    parameter_meta {
        loci:          "Loci to assemble. At least one is required. Reads from all loci will be merged for assembly. Format: [\"chr1:1000-2000\", \"chr1:5000-10000\"]"
        region_size:   "Estimated size of region to assemble, in mb"
        aligned_bam:   "aligned file"
        aligned_bai:   "index file"
        prefix:        "prefix for output files"
        experiment_type:  "CLR, CCS, or ONT"

#        add_unaligned_reads: "set to true to include unaligned reads in the assembly (default: false)"
        run_quast:           "set to true to run Quast on the assembly (default: false)"

        ref_map_file:  "table indicating reference sequence and auxillary file locations"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    Map[String, String] map_presets = {
        'CLR':    'pacbio',
        'CCS':    'pacbio-hifi',
        'ONT':    'nanopore'
    }

    String preset = map_presets[experiment_type]

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


    if (experiment_type == 'CLR' || experiment_type == 'ONT') {
        ######## Call 3 stages of Canu #########
        call Canu.Correct as Correct {
            input:
                reads = BamToFastq.reads_fq,
                genome_size = region_size,
                prefix = prefix,
                preset = preset
        }

        call Canu.Trim as Trim {
            input:
                genome_size = region_size,
                corrected_reads = Correct.corrected_reads,
                prefix = prefix,
                preset = preset
        }

        call Canu.Assemble as Assemble {
            input:
                genome_size = region_size,
                trimmed_reads = Trim.trimmed_reads,
                prefix = prefix,
                preset = preset
        }
        ######## Done calling Canu #########
        File assembled_contigs = Assemble.canu_contigs_fasta
    }

    if (experiment_type == 'CCS') {
        ######## Call Canu in one step - this is good for CCS reads that don't need correcting and trimming #########
        call Canu.SingleStep as SingleStep {
            input:
                genome_size = region_size,
                reads = BamToFastq.reads_fq,
                prefix = prefix,
                preset = preset
        }
        File assembled_contigs = SingleStep.canu_contigs_fasta
    }

    if (run_quast) {
        call Quast.Quast {
            input:
                ref = ref_map['fasta'],
                assemblies = [ assembled_contigs ]
        }
    }

    call AV.CallAssemblyVariants as CallAssemblyVariants {
        input:
            asm_fasta = assembled_contigs,
            ref_fasta = ref_map['fasta'],
            participant_name = prefix,
            prefix = prefix + ".canu"
    }

    output {
        File local_bam = subset_bam
        File canu_fa = assembled_contigs
        File variants = CallAssemblyVariants.paftools_vcf
        File paf_aln = CallAssemblyVariants.paf
    }
}