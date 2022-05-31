version 1.0

import "tasks/AlignReads.wdl" as AR
import "tasks/Quast.wdl" as Quast
import "tasks/CallAssemblyVariants.wdl" as  CallAssemblyVariants


workflow mito_paf {

    input {
        File asm_fasta
        File ref_fasta
        String participant_name
        String prefix
    }

    parameter_meta {
        asm_fasta:        "haploid assembly"
        ref_fasta:        "reference to which assembly should be aligned"
        prefix: "prefix for output"
        participant_name: "participant name"
    }



    call CallAssemblyVariants.CallAssemblyVariants as CallAssemblyVariants {input:
                                                                                 asm_fasta = asm_fasta,
                                                                                ref_fasta = ref_fasta,
                                                                                participant_name = participant_name,
                                                                                prefix = prefix
                                                                            }

    output {
        File paf = CallAssemblyVariants.paf
        File paftools_vcf = CallAssemblyVariants.paftools_vcf
    }
}