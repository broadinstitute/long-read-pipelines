import "tasks/AlignReads.wdl" as AR
import "tasks/Quast.wdl" as Quast
import "tasks/CallAssemblyVariants.wdl" as  CallAssemblyVariants


workflow mito_paf{

    meta {
        description:
        "this workflow calls variants in mitochondria genome"
    }
    input{
        File asm_fasta
        File ref_fasta
        File participant_name
        File prefix

    }
    parameter_meta{
        prefix: "prefix for output"

    }

    call CallAssemblyVariants.CallAssemblyVariants as  CallAssemblyVariants {input:
                                                                                 asm_fasta = asm_fasta,
                                                                                ref_fasta = ref_fasta,
                                                                                participant_name = participant_name,
                                                                                prefix = prefix}
}