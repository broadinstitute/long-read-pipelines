version 1.0

import "../../../tasks/VariantCalling/CallAssemblyVariants.wdl" as AV


workflow PaftoolAsmVariants{
    meta{
        description: "a workflow that perform variant calling from genome assemblies"
    }
    input{
        File assembly_hap1
        File assembly_hap2
        File reference
        String prefix
        String samplename
    }
    call AV.CallAssemblyVariants as asmvarhap1 {
        input:
            asm_fasta = assembly_hap1,
            ref_fasta = reference,
            participant_name = samplename + "hap1",
            prefix = prefix
    }

    call AV.CallAssemblyVariants as asmvarhap2 {
        input:
            asm_fasta = assembly_hap2,
            ref_fasta = reference,
            participant_name = samplename + "hap2",
            prefix = prefix
    }

    output {

        File paf_hap1 = asmvarhap1.paf
        File paftools_vcf_hap1 = asmvarhap1.paftools_vcf
        File paf_hap2 = asmvarhap2.paf
        File paftools_vcf_hap2 = asmvarhap2.paftools_vcf
    }
}