version 1.0

import "ONTPepper.wdl" as ONP

workflow DVPPipeline {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File ref_scatter_interval_list_locator
        File ref_scatter_interval_list_ids

        Int pepper_memory
        Boolean sort_gVCF

        String vcf_output_prefix
    }

    parameter_meta {
        ref_scatter_interval_list_locator: "A file holding paths to interval_list files"
        ref_scatter_interval_list_ids:     "A file that gives short IDs to the interval_list files"
    }

    Array[String] interval_list_files = read_lines(ref_scatter_interval_list_locator)
    Array[String] interval_list_ids   = read_lines(ref_scatter_interval_list_ids)
    Array[Pair[String, String]] ided_interval_list_files = zip(interval_list_ids, interval_list_files)

    scatter (pair in ided_interval_list_files) {
        call ONP.SubsetBam {
            input:
                bam = bam,
                bai = bai,
                interval_list_file = pair.right, 
                interval_id = pair.left,
                prefix = basename(bam, ".bam")
        }
        
        call ONP.PEPPER {
            input:
                bam           = SubsetBam.subset_bam,
                bai           = SubsetBam.subset_bai,
                ref_fasta     = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                memory        = pepper_memory
        }
    }

    String local_prefix = vcf_output_prefix + ".deepvariant_pepper"

    call ONP.SortVCFsFast as MergeDeepVariantGVCFs {
        input:
            vcfs     = PEPPER.gVCF,
            tbis     = PEPPER.gVCF_tbi,
            prefix   = local_prefix + ".g",
            sort     = sort_gVCF
    }

    call ONP.SortVCFsFast as MergeDeepVariantPhasedVCFs {
        input:
            vcfs     = PEPPER.phasedVCF,
            tbis     = PEPPER.phasedtbi,
            sort     = true,
            prefix   = local_prefix + ".phased"
    }

    call ONP.SortVCFsFast as MergeDeepVariantVCFs {
        input:
            vcfs     = PEPPER.VCF,
            tbis     = PEPPER.VCF_tbi,
            sort     = true,
            prefix   = local_prefix
    }

    output {
        File phased_vcf = MergeDeepVariantPhasedVCFs.vcf
        File phased_tbi = select_first([MergeDeepVariantPhasedVCFs.tbi, bai]) # hack, must exist

        File g_vcf = MergeDeepVariantGVCFs.vcf
        File? g_tbi = MergeDeepVariantGVCFs.tbi

        File vcf = MergeDeepVariantVCFs.vcf
        File tbi = select_first([MergeDeepVariantVCFs.tbi, bai]) # hack, must exist
    }
}
