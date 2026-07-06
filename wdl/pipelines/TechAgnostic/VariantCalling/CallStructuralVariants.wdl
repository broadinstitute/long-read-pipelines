version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../structs/ReferenceMetadata.wdl"
import "../../../structs/Structs.wdl"

import "../../../tasks/VariantCalling/PBSV.wdl"
import "../../../tasks/VariantCalling/Sniffles2.wdl" as Sniffles2

struct SVCallingConfig {
    Int min_sv_len
    Boolean pbsv_discover_per_chr

    RuntimeAttr? pbsv_discover_runtime_attr_override
    RuntimeAttr? pbsv_call_runtime_attr_override
    RuntimeAttr? sniffles_runtime_attr_override

    String? gcp_zones
}

workflow Work {
    meta {
        description: "Call structual variants using reads-based methods (i.e. not for assembly-contig-based methods)."
    }
    parameter_meta {
        is_hifi: "Indicate if the input is HiFi data"
        is_ont: "If the input data is ONT"
        per_chr_bam_bai_and_id: "Must be provided when pbsv_discover_per_chr is true."
        pbsv_discover_per_chr: "To run the discover stage of PBSV in per-chromosome style or not. If true, then the WGS bam must be sharded accordingly beforehand."
        minsvlen: "Minimum SV length in bp; only affects Sniffles 2 calls."
    }

    input {
        String gcs_out_dir

        # sample info
        File bam
        File bai
        String prefix

        Boolean is_hifi
        Boolean is_ont

        Boolean pbsv_discover_per_chr
        Array[Pair[String, Pair[File, File]]]? per_chr_bam_bai_and_id

        # reference info
        File ref_bundle_json_file

        # sv-specific args
        Int minsvlen = 20
        RuntimeAttr? pbsv_discover_runtime_attr_override
        RuntimeAttr? pbsv_call_runtime_attr_override
        RuntimeAttr? sniffles_runtime_attr_override

        # optimization
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }

    output {
        File sniffles_vcf = FinalizeSnifflesVcf.gcs_path
        File sniffles_tbi = FinalizeSnifflesTbi.gcs_path
        File sniffles_snf = FinalizeSnifflesSnf.gcs_path

        File pbsv_vcf = FinalizePBSVvcf.gcs_path
        File pbsv_tbi = FinalizePBSVtbi.gcs_path
    }

    if (pbsv_discover_per_chr) {
        if (!defined(per_chr_bam_bai_and_id)) {
            call Utils.StopWorkflow { input: reason = "When calling PBSV to work on chromosomes separately, must also provide a list of BAMs sharded by chromosomes"}
        }
    }

    HumanReferenceBundle ref_bundle = read_json(ref_bundle_json_file)
    ##########################################################
    # Sniffles-2
    ##########################################################
    call Utils.InferSampleName { input: bam = bam, bai = bai }
    call Sniffles2.SampleSV as Sniffles2SV {
        input:
            bam    = bam,
            bai    = bai,
            minsvlen = minsvlen,
            sample_id = InferSampleName.sample_name,
            prefix = prefix,
            tandem_repeat_bed = ref_bundle.tandem_repeat_bed,
            runtime_attr_override = sniffles_runtime_attr_override
    }

    ##########################################################
    # PBSV
    ##########################################################
    if (pbsv_discover_per_chr) {

        scatter (triplet in select_first([per_chr_bam_bai_and_id])) {
            String contig = triplet.left
            File shard_bam = triplet.right.left
            File shard_bai = triplet.right.right

            call PBSV.Discover as pbsv_discover_chr {
                input:
                    bam = shard_bam,
                    bai = shard_bai,
                    is_hifi = is_hifi,
                    is_ont = is_ont,
                    chr = contig,
                    prefix = prefix,
                    ref_fasta = ref_bundle.fasta,
                    ref_fasta_fai = ref_bundle.fai,
                    tandem_repeat_bed = ref_bundle.tandem_repeat_bed,
                    zones = zones,
                    runtime_attr_override = pbsv_discover_runtime_attr_override
            }
        }

        call PBSV.Call as pbsv_wg_call {
            input:
                svsigs = pbsv_discover_chr.svsig,
                ref_fasta = ref_bundle.fasta,
                ref_fasta_fai = ref_bundle.fai,
                is_hifi = is_hifi,
                is_ont = is_ont,
                prefix = prefix,
                zones = zones,
                runtime_attr_override = pbsv_call_runtime_attr_override
        }
    }

    if (!pbsv_discover_per_chr) {
        call PBSV.RunPBSV as PBSVslow {
            input:
                bam = bam,
                bai = bai,
                prefix = prefix,
                is_hifi = is_hifi,
                is_ont = is_ont,
                ref_fasta = ref_bundle.fasta,
                ref_fasta_fai = ref_bundle.fai,
                tandem_repeat_bed = ref_bundle.tandem_repeat_bed,
                zones = zones,
                discover_runtime_attr_override = pbsv_discover_runtime_attr_override,
                call_runtime_attr_override = pbsv_call_runtime_attr_override
        }
    }

    ##########################################################
    # Save files
    ##########################################################
    String svdir = sub(gcs_out_dir, "/$", "")

    call FF.FinalizeToFile as FinalizePBSVvcf { input: outdir = svdir, file = select_first([pbsv_wg_call.vcf, PBSVslow.vcf]) }
    call FF.FinalizeToFile as FinalizePBSVtbi { input: outdir = svdir, file = select_first([pbsv_wg_call.tbi, PBSVslow.tbi]) }

    call FF.FinalizeToFile as FinalizeSnifflesVcf { input: outdir = svdir, file = Sniffles2SV.vcf }
    call FF.FinalizeToFile as FinalizeSnifflesTbi { input: outdir = svdir, file = Sniffles2SV.tbi }
    call FF.FinalizeToFile as FinalizeSnifflesSnf { input: outdir = svdir, file = Sniffles2SV.snf }
}
