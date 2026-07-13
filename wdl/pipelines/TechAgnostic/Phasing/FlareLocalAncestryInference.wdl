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

        Boolean run_flare = true
        Boolean enable_thinning = true
        Boolean training_mode = false

        String output_prefix

        Float maf = 0.01
        Int thin_bp = 20000
        Int n_sample_shards = 10
        Int seed = -99999
        Int nthreads = 16
        Int mem_gb = 64
    }

    Boolean em = training_mode
    Int effective_shards = if training_mode then 1 else n_sample_shards

    String prep_prefix = output_prefix + ".prep"

    if (enable_thinning) {
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
    }

    if (!enable_thinning) {
        call Flare.PrepIntegratedVcfForFlare as PrepIntegrated {
            input:
                joint_vcf = test_vcf,
                ref_vcf = ref_vcf,
                ref_panel = ref_panel,
                prefix = prep_prefix + ".ready"
        }
    }

    File ready_gt_vcf = select_first([FilterSites.gt_vcf_out, PrepIntegrated.gt_vcf_out])
    File ready_gt_vcf_index = select_first([FilterSites.gt_vcf_csi, PrepIntegrated.gt_vcf_csi])
    File ready_ref_vcf = select_first([FilterSites.ref_vcf_out, PrepIntegrated.ref_vcf_out])
    File ready_ref_vcf_index = select_first([FilterSites.ref_vcf_csi, PrepIntegrated.ref_vcf_csi])

    if (run_flare) {
        if (effective_shards > 1) {
            call Flare.SplitStudySamplesForFlare as SplitSamples {
                input:
                    gt_vcf = ready_gt_vcf,
                    n_shards = effective_shards,
                    prefix = prep_prefix + ".samples"
            }

            scatter (shard_idx in range(length(SplitSamples.sample_lists))) {
                call Flare.Flare as F_shard {
                    input:
                        ref_vcf = ready_ref_vcf,
                        ref_vcf_index = ready_ref_vcf_index,
                        ref_panel = ref_panel,
                        test_vcf = ready_gt_vcf,
                        test_vcf_index = ready_gt_vcf_index,
                        plink_map = plink_map,
                        output_prefix = output_prefix + ".shard" + shard_idx,
                        em = em,
                        flare_model = flare_model,
                        gt_samples = SplitSamples.sample_lists[shard_idx],
                        seed = seed,
                        nthreads = nthreads,
                        mem_gb = mem_gb
                }
            }

            call Flare.MergeFlareShardOutputs as MergeShards {
                input:
                    anc_vcfs = F_shard.anc_vcf,
                    global_anc_files = F_shard.global_anc,
                    output_prefix = output_prefix
            }
        }

        if (effective_shards == 1) {
            call Flare.Flare as F {
                input:
                    ref_vcf = ready_ref_vcf,
                    ref_vcf_index = ready_ref_vcf_index,
                    ref_panel = ref_panel,
                    test_vcf = ready_gt_vcf,
                    test_vcf_index = ready_gt_vcf_index,
                    plink_map = plink_map,
                    output_prefix = output_prefix,
                    em = em,
                    flare_model = flare_model,
                    seed = seed,
                    nthreads = nthreads,
                    mem_gb = mem_gb
                }

            if (training_mode) {
                File em_model = F.model
            }
        }
    }

    output {
        File? global_anc = select_first([MergeShards.global_anc, F.global_anc])
        File? anc_vcf = select_first([MergeShards.anc_vcf, F.anc_vcf])
        File? anc_vcf_csi = MergeShards.anc_vcf_csi
        File? log = F.log
        Array[File]? shard_logs = F_shard.log
        File? model = em_model
    }
}
