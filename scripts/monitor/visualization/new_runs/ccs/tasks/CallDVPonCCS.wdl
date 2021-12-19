version 1.0

import "CCSPepper.wdl" as ONS

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
        call ONS.SubsetBam {
            input:
                bam = bam,
                bai = bai,
                interval_list_file = pair.right, 
                interval_id = pair.left,
                prefix = basename(bam, ".bam")
        }
        
        call ONS.PEPPER {
            input:
                bam           = SubsetBam.subset_bam,
                bai           = SubsetBam.subset_bai,
                ref_fasta     = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                memory        = pepper_memory
        }

        call ONS.DV {
            input:
                bam           = PEPPER.phased_bam,
                bai           = PEPPER.phased_bai,
                ref_fasta     = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                memory        = pepper_memory
        }
    }

    String local_prefix = vcf_output_prefix + ".deepvariant_pepper"

    call ONS.SortVCFs as MergeDeepVariantVCFs {
        input:
            vcfs     = DV.VCF,
            tbis     = DV.VCF_tbi,
            ref_dict = ref_dict,
            prefix   = local_prefix
    }

    call ONS.SortVCFsFast as MergeDeepVariantGVCFs {
        input:
            vcfs     = DV.gVCF,
            tbis     = DV.gVCF_tbi,
            prefix   = local_prefix + ".g", 
            sort     = sort_gVCF
    }

    output {

        File g_vcf = MergeDeepVariantGVCFs.vcf
        File? g_tbi = MergeDeepVariantGVCFs.tbi

        File vcf = MergeDeepVariantVCFs.vcf
        File tbi = MergeDeepVariantVCFs.tbi
    }
}
