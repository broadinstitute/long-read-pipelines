version 1.0

import "../Utility/Utils.wdl"

import "PBSV.wdl"
import "Sniffles2.wdl" as Sniffles2

workflow Work {
    meta {
        description: "Call structual variants using reads-based methods (i.e. not for assembly-contig-based methods)."
    }
    parameter_meta {
        is_hifi: "Indicate if the input is HiFi data"
        is_ont: "If the input data is ONT"
        per_chr_bam_bai_and_id: "Must be provided when pbsv_discover_per_chr is true."
        pbsv_discover_per_chr: "To run the discover stage of PBSV in per-chromosome style or not. If true, then the WGS bam must be sharded accordingly beforehand."
    }

    input {
        # sample info
        File bam
        File bai
        String prefix

        Boolean is_hifi
        Boolean is_ont

        Boolean pbsv_discover_per_chr
        Array[Pair[String, Pair[File, File]]]? per_chr_bam_bai_and_id

        # reference info
        Map[String, String] ref_map

        # sv-specific args
        Int minsvlen = 50

        # optimization
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }

    output {
        File sniffles_vcf = Sniffles2SV.vcf
        File sniffles_tbi = Sniffles2SV.tbi
        File sniffles_snf = Sniffles2SV.snf

        File pbsv_vcf = select_first([pbsv_wg_call.vcf, PBSVslow.vcf])
        File pbsv_tbi = select_first([pbsv_wg_call.tbi, PBSVslow.tbi])
    }

    if (pbsv_discover_per_chr) {
        if (!defined(per_chr_bam_bai_and_id)) {
            call Utils.StopWorkflow { input: reason = "When calling PBSV to work on chromosomes separately, must also provide a list of BAMs sharded by chromosomes"}
        }
    }

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
            tandem_repeat_bed = ref_map['tandem_repeat_bed']
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
                    ref_fasta = ref_map['fasta'],
                    ref_fasta_fai = ref_map['fai'],
                    tandem_repeat_bed = ref_map['tandem_repeat_bed'],
                    zones = zones
            }
        }

        call PBSV.Call as pbsv_wg_call {
            input:
                svsigs = pbsv_discover_chr.svsig,
                ref_fasta = ref_map['fasta'],
                ref_fasta_fai = ref_map['fai'],
                is_hifi = is_hifi,
                is_ont = is_ont,
                prefix = prefix,
                zones = zones
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
                ref_fasta = ref_map['fasta'],
                ref_fasta_fai = ref_map['fai'],
                tandem_repeat_bed = ref_map['tandem_repeat_bed'],
                zones = zones
        }
    }
}
