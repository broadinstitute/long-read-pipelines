version 1.0

import "../../../tasks/VariantCalling/CallAssemblyVariants.wdl" as AV


workflow PaftoolAsmVariants{
    meta{
        description: "a workflow that perform variant calling from genome assemblies"
    }
    input{
        File assembly
        File reference
        String prefix
        String samplename
    }
    call AV.CallAssemblyVariants {
        input:
            asm_fasta = assembly,
            ref_fasta = reference,
            participant_name = samplename,
            prefix = prefix + ".canu"
    }

    output {

        File paf = CallAssemblyVariants.paf
        File paftools_vcf = CallAssemblyVariants.paftools_vcf
    }
}