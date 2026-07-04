version 1.0

import "../../../tasks/Phasing/Flare.wdl"
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow FlareLocalAncestryInference {
    input {
        File ref_vcf
        File ref_vcf_index
        File ref_panel
        File test_vcf
        File plink_map
        File ref_fasta
        File ref_fasta_fai

        String chromosome
        Boolean is_chr_x = false
        File? sample_sex_map
        File? flare_model

        String gcs_out_root_dir
        String output_prefix

        Int nthreads = 16
        Int mem_gb = 64
    }

    Boolean em = !defined(flare_model)

    String outdir = sub(gcs_out_root_dir, "/+$", "") + "/FlareLocalAncestryInference/~{output_prefix}"

    call Flare.FilterVCFsForFlare as FilterVCFs {
        input:
            joint_vcf = test_vcf,
            ref_vcf = ref_vcf,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            chromosome = chromosome,
            prefix = output_prefix + ".prep",
            is_chr_x = is_chr_x,
            sample_sex_map = sample_sex_map
    }

    call Flare.Flare as F {
        input:
            ref_vcf = FilterVCFs.ref_vcf_out,
            ref_vcf_index = FilterVCFs.ref_vcf_csi,
            ref_panel = ref_panel,
            test_vcf = FilterVCFs.gt_vcf,
            test_vcf_index = FilterVCFs.gt_vcf_csi,
            plink_map = plink_map,
            output_prefix = output_prefix,
            em = em,
            flare_model = flare_model,
            nthreads = nthreads,
            mem_gb = mem_gb
    }

    call FF.FinalizeToFile as FinalizeGlobalAnc {
        input:
            outdir = outdir,
            file = F.global_anc
    }

    call FF.FinalizeToFile as FinalizeAncVCF {
        input:
            outdir = outdir,
            file = F.anc_vcf
    }

    call FF.FinalizeToFile as FinalizeLog {
        input:
            outdir = outdir,
            file = F.log
    }

    if (em) {
        call FF.FinalizeToFile as FinalizeModel {
            input:
                outdir = outdir,
                file = F.model
        }
    }

    output {
        File global_anc = FinalizeGlobalAnc.gcs_path
        File anc_vcf = FinalizeAncVCF.gcs_path
        File log = FinalizeLog.gcs_path
        File? model = FinalizeModel.gcs_path
    }
}
