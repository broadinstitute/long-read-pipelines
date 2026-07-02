version 1.0

import "ONTPepper.wdl"

import "../../tasks/Utility/VariantUtils.wdl"
import "../../tasks/Utility/Utils.wdl"
import "../../tasks/Alignment/WhatsHap.wdl"

workflow Run {
    meta {
        desciption:
        "Runs Clair3 on the input (sharded) BAM."
    }
    parameter_meta {
        phase_and_tag: "having this turned off means phased VCF and haplotagged BAM will not be output."
        how_to_shard_wg_for_calling: "An array of the BAM's shard; each element is assumed to be a tuple of (ID for the shard, (BAM of the shard, BAI of the shard))"
        prefix: "Prefix for output files"
        model_for_pepper_margin_dv: "refer to https://github.com/kishwarshafin/pepper for appropriate values"
    }

    input {
        Array[Pair[String, Pair[File, File]]] how_to_shard_wg_for_calling
        Array[Int] memup_shards
        Boolean all_memup
        String prefix
        String model_for_pepper_margin_dv

        Map[String, String] ref_map

        # read-haplotaging desired or not
        Boolean phase_and_tag

        # optimization
        Int dv_threads
        Int dv_memory
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    output {
        File legacy_ont_dvp_g_vcf = MergePEPPERGVCFs.vcf
        File legacy_ont_dvp_g_tbi = MergePEPPERGVCFs.tbi
        File? legacy_ont_dvp_phased_vcf = MergePEPPERPhasedVCFs.vcf
        File? legacy_ont_dvp_phased_tbi = MergePEPPERPhasedVCFs.tbi
        File? legacy_ont_dvp_haplotagged_bam = MergePEPPERHapTaggedBam.merged_bam
        File? legacy_ont_dvp_haplotagged_bai = MergePEPPERHapTaggedBam.merged_bai
        File? legacy_ont_dvp_phased_vcf_stats_tsv = ONTPhaseStatsLegacy.stats_tsv
        File? legacy_ont_dvp_phased_vcf_stats_gtf = ONTPhaseStatsLegacy.stats_gtf
    }

    Int N = length(how_to_shard_wg_for_calling)
    scatter (i in range(N)) {
        Pair[String, Pair[File, File]] triplet = how_to_shard_wg_for_calling[i]
        if (triplet.left != "alts") {
            if (!all_memup) {
                scatter (j in memup_shards) { Int indicator = if (j == i) then 1 else 0 }
                call SumInts as SumIndicators { input: integers = indicator }
                Boolean memup_for_shard = SumIndicators.total > 0
            }
            Boolean memup = all_memup || select_first([memup_for_shard, false])

            call ONTPepper.Pepper {
                input:
                    bam           = triplet.right.left,
                    bai           = triplet.right.right,
                    ref_fasta     = ref_map['fasta'],
                    ref_fasta_fai = ref_map['fai'],
                    model         = model_for_pepper_margin_dv,
                    phase_and_tag = phase_and_tag,
                    threads       = dv_threads,
                    memory        = if (!memup) then dv_memory else 6*dv_threads,
                    zones         = zones
            }
        }
    }

    String pepper_prefix = prefix + ".PEPPER-Margin-DeepVariant"

    call VariantUtils.MergeAndSortVCFs as MergePEPPERGVCFs { input:
        vcfs     = select_all(Pepper.gVCF),
        prefix   = pepper_prefix + ".g",
        ref_fasta_fai = ref_map['fai'],
        # Right-size: this merge peaked at ~2.2 GiB RAM / ~3 cores on the task's
        # default 48 GB / 8-core VM (bcftools sort spills to --temp-dir, so RAM
        # stays low). Trim CPU+RAM only; disk is left at the task default since
        # the sort's temp spill needs the headroom.
        runtime_attr_override = object { cpu_cores: 4, mem_gb: 8 }
    }

    if (phase_and_tag) {
        call VariantUtils.MergeAndSortVCFs as MergePEPPERPhasedVCFs { input:
            vcfs     = select_all(Pepper.phasedVCF),
            prefix   = pepper_prefix + ".phased",
            ref_fasta_fai = ref_map['fai']
        }

        call Utils.MergeBams as MergePEPPERHapTaggedBam { input:
            bams     = select_all(Pepper.hap_tagged_bam),
            prefix   = prefix +  ".MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged"
        }

        call WhatsHap.Stats as ONTPhaseStatsLegacy { input:
            phased_vcf=MergePEPPERPhasedVCFs.vcf, phased_tbi=MergePEPPERPhasedVCFs.tbi
        }
    }
}

task SumInts {
    input {
        Array[Int] integers
    }
    output {
        Int total = read_int("result.txt")
    }

    # WDL helper to write one Int per line to a file
    File int_file = write_lines(integers)

    command <<<
        # Use awk to sum the file line-by-line (memory efficient!)
        awk '{s+=$1} END {print s}' ~{int_file} > result.txt
    >>>
    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
