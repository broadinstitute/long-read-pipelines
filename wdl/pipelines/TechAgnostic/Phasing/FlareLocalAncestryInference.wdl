version 1.0

import "../../../tasks/Phasing/Flare.wdl"

workflow FlareLocalAncestryInference {
    input {
        File ref_vcf
        File ref_vcf_index
        File ref_panel
        File test_vcf
        File plink_map
        File ref_fasta
        File ref_fasta_fai

        File? flare_model

        String output_prefix

        Float maf = 0.01
        Int thin_bp = 20000
        Int nthreads = 16
        Int mem_gb = 64
    }

    Boolean em = !defined(flare_model)

    String prep_prefix = output_prefix + ".prep"

    call Flare.PrepStudyVcfForFlare as PrepStudy {
        input:
            joint_vcf = test_vcf,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            prefix = prep_prefix + ".study",
            maf = maf
    }

    call Flare.PrepRefVcfForFlare as PrepRef {
        input:
            ref_vcf = ref_vcf,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            prefix = prep_prefix + ".ref"
    }

    call Flare.IntersectVCFsForFlare as IntersectSites {
        input:
            gt_vcf = PrepStudy.gt_bcf,
            gt_vcf_index = PrepStudy.gt_bcf_csi,
            ref_vcf = PrepRef.ref_bcf,
            ref_vcf_index = PrepRef.ref_bcf_csi,
            prefix = prep_prefix + ".isec"
    }

    call Flare.ThinVCFsForFlare as ThinSites {
        input:
            gt_vcf = IntersectSites.gt_vcf_out,
            ref_vcf = IntersectSites.ref_vcf_out,
            prefix = prep_prefix + ".thin",
            thin_bp = thin_bp
    }

    call Flare.IntersectVCFsForFlare as FinalizeSites {
        input:
            gt_vcf = ThinSites.gt_vcf_out,
            gt_vcf_index = ThinSites.gt_vcf_csi,
            ref_vcf = ThinSites.ref_vcf_out,
            ref_vcf_index = ThinSites.ref_vcf_csi,
            prefix = prep_prefix + ".flare"
    }

    call Flare.FilterFlareReadySites as FilterSites {
        input:
            gt_vcf = FinalizeSites.gt_vcf_out,
            gt_vcf_index = FinalizeSites.gt_vcf_csi,
            ref_vcf = FinalizeSites.ref_vcf_out,
            ref_vcf_index = FinalizeSites.ref_vcf_csi,
            ref_panel = ref_panel,
            prefix = prep_prefix + ".ready"
    }

    call Flare.Flare as F {
        input:
            ref_vcf = FilterSites.ref_vcf_out,
            ref_vcf_index = FilterSites.ref_vcf_csi,
            ref_panel = ref_panel,
            test_vcf = FilterSites.gt_vcf_out,
            test_vcf_index = FilterSites.gt_vcf_csi,
            plink_map = plink_map,
            output_prefix = output_prefix,
            em = em,
            flare_model = flare_model,
            nthreads = nthreads,
            mem_gb = mem_gb
    }

    if (em) {
        File em_model = F.model
    }

    output {
        File global_anc = F.global_anc
        File anc_vcf = F.anc_vcf
        File log = F.log
        File? model = em_model
    }
}
