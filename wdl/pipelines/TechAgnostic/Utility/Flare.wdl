version 1.0

import "../../../tasks/Utility/FlareUtils.wdl" as FlareUtils
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow Flare {

    meta {
        description: "Run FLARE local ancestry inference for a single chromosome on SHAPEIT4-ligated study BCFs."
    }

    parameter_meta {
        gt_bcf: "Single-chromosome study BCF (e.g. aou_lr_phase2_v1.chr20.shapeit4.ligated.bcf)"
        ref_vcf: "gnomAD HGDP+1KG reference VCF for this chromosome"
        ref_vcf_tbi: "Tabix index for reference VCF"
        ref_panel: "FLARE ref-panel file (sample to population map)"
        genetic_map: "Beagle PLINK genetic map for this chromosome (chr-prefixed)"
        ref_fasta: "GRCh38 reference FASTA"
        chromosome: "Chromosome name (e.g. chr20 or chrX)"
        is_chr_x: "Whether this run is for chrX"
        sample_sex_map: "Optional sample sex map for chrX runs"
        flare_model: "Optional chr1 .model file; if unset, em=true"
        gcs_out_root_dir: "GCS directory for output files"
        prefix: "Output prefix"
        nthreads: "Number of FLARE threads"
        mem_gb: "Java heap size in GB for FLARE"
    }

    input {
        File gt_bcf
        File ref_vcf
        File ref_vcf_tbi
        File ref_panel
        File genetic_map
        File ref_fasta

        String chromosome
        Boolean is_chr_x = false
        File? sample_sex_map
        File? flare_model

        String gcs_out_root_dir
        String prefix

        Int nthreads = 16
        Int mem_gb = 128
    }

    Boolean em = !defined(flare_model)

    String outdir = sub(gcs_out_root_dir, "/+$", "") + "/Flare/~{prefix}"

    call FlareUtils.FilterVCFsForFlare as FilterVCFs {
        input:
            gt_bcf = gt_bcf,
            ref_vcf = ref_vcf,
            ref_fasta = ref_fasta,
            chromosome = chromosome,
            prefix = prefix + ".prep",
            is_chr_x = is_chr_x,
            sample_sex_map = sample_sex_map
    }

    call FlareUtils.RunFlare as RunFlare {
        input:
            gt_vcf = FilterVCFs.gt_vcf,
            ref_vcf = FilterVCFs.ref_vcf_out,
            ref_panel = ref_panel,
            genetic_map = genetic_map,
            prefix = prefix,
            em = em,
            flare_model = flare_model,
            nthreads = nthreads,
            mem_gb = mem_gb
    }

    call FF.FinalizeToFile as FinalizeGlobalAnc {
        input:
            outdir = outdir,
            file = RunFlare.global_anc
    }

    call FF.FinalizeToFile as FinalizeAncVCF {
        input:
            outdir = outdir,
            file = RunFlare.anc_vcf
    }

    call FF.FinalizeToFile as FinalizeLog {
        input:
            outdir = outdir,
            file = RunFlare.log
    }

    if (em) {
        call FF.FinalizeToFile as FinalizeModel {
            input:
                outdir = outdir,
                file = RunFlare.model
        }
    }

    output {
        File global_anc = FinalizeGlobalAnc.gcs_path
        File anc_vcf = FinalizeAncVCF.gcs_path
        File log = FinalizeLog.gcs_path
        File? model = FinalizeModel.gcs_path
    }
}
