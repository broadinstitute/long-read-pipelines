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

    # toggle individual SV callers; both default to true when omitted
    Boolean? run_sniffles
    Boolean? run_pbsv

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
        run_sniffles: "Whether to run Sniffles2. Set false (with run_pbsv true) to call PBSV only."
        run_pbsv: "Whether to run PBSV. Set false (with run_sniffles true) to call Sniffles2 only."
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

        # which callers to run; either (but not neither) may be turned off
        Boolean run_sniffles = true
        Boolean run_pbsv = true

        # reference info
        File ref_bundle_json_file

        # sv-specific args
        Int minsvlen = 20
        RuntimeAttr? pbsv_discover_runtime_attr_override
        RuntimeAttr? pbsv_call_runtime_attr_override
        RuntimeAttr? sniffles_runtime_attr_override

        # per-shard memory bump for pbsv discover; shard indices match the order of per_chr_bam_bai_and_id
        Array[Int] pbsv_discover_memup_shards = []
        Int? pbsv_discover_memup_gb

        # optimization
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }

    output {
        File? sniffles_vcf = FinalizeSnifflesVcf.gcs_path
        File? sniffles_tbi = FinalizeSnifflesTbi.gcs_path
        File? sniffles_snf = FinalizeSnifflesSnf.gcs_path

        File? pbsv_vcf = FinalizePBSVvcf.gcs_path
        File? pbsv_tbi = FinalizePBSVtbi.gcs_path
    }

    if (!run_sniffles && !run_pbsv) {
        call Utils.StopWorkflow as NoCallerRequested { input: reason = "Both Sniffles and PBSV were turned off; nothing to do for SV calling."}
    }

    HumanReferenceBundle ref_bundle = read_json(ref_bundle_json_file)
    String svdir = sub(gcs_out_dir, "/$", "")

    ##########################################################
    # Sniffles-2
    ##########################################################
    if (run_sniffles) {
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

        call FF.FinalizeToFile as FinalizeSnifflesVcf { input: outdir = svdir, file = Sniffles2SV.vcf }
        call FF.FinalizeToFile as FinalizeSnifflesTbi { input: outdir = svdir, file = Sniffles2SV.tbi }
        call FF.FinalizeToFile as FinalizeSnifflesSnf { input: outdir = svdir, file = Sniffles2SV.snf }
    }

    ##########################################################
    # PBSV
    ##########################################################
    if (run_pbsv) {
        if (pbsv_discover_per_chr && !defined(per_chr_bam_bai_and_id)) {
            call Utils.StopWorkflow { input: reason = "When calling PBSV to work on chromosomes separately, must also provide a list of BAMs sharded by chromosomes"}
        }

        if (pbsv_discover_per_chr) {

            Array[Pair[String, Pair[File, File]]] pbsv_discover_shards = select_first([per_chr_bam_bai_and_id])

            scatter (i in range(length(pbsv_discover_shards))) {
                Pair[String, Pair[File, File]] triplet = pbsv_discover_shards[i]
                String contig = triplet.left
                File shard_bam = triplet.right.left
                File shard_bai = triplet.right.right

                scatter (j in pbsv_discover_memup_shards) { Int pbsv_discover_memup_indicator = if (j == i) then 1 else 0 }
                call SumInts as SumPbsvDiscoverMemupIndicators { input: integers = pbsv_discover_memup_indicator }
                Boolean pbsv_discover_memup_for_shard = SumPbsvDiscoverMemupIndicators.total > 0

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
                        runtime_attr_override = pbsv_discover_runtime_attr_override,
                        memup = pbsv_discover_memup_for_shard,
                        memup_gb = pbsv_discover_memup_gb
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

        call FF.FinalizeToFile as FinalizePBSVvcf { input: outdir = svdir, file = select_first([pbsv_wg_call.vcf, PBSVslow.vcf]) }
        call FF.FinalizeToFile as FinalizePBSVtbi { input: outdir = svdir, file = select_first([pbsv_wg_call.tbi, PBSVslow.tbi]) }
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
