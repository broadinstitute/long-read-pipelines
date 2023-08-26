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
        tandem_repeat_bed: "BED file containing TRF finder for better SV calls (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"
        pbsv_discover_per_chr: "To run the discover stage of PBSV in per-chromosome style or not. If true, then the WGS bam must be sharded accordingly beforehand."
        per_chr_bam_bai_and_id: "Must be provided when pbsv_discover_per_chr is true."
    }
    input {
        # sample info
        File bam
        File bai
        String prefix

        Boolean is_hifi
        Boolean is_ont

        Array[Pair[String, Pair[File, File]]]? per_chr_bam_bai_and_id

        # reference info
        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        # sv-specific args
        Int minsvlen = 50
        File? tandem_repeat_bed

        Boolean pbsv_discover_per_chr

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

    call Utils.InferSampleName { input: bam = bam, bai = bai }
    call Sniffles2.SampleSV as Sniffles2SV {
        input:
            bam    = bam,
            bai    = bai,
            minsvlen = minsvlen,
            sample_id = InferSampleName.sample_name,
            prefix = prefix,
            tandem_repeat_bed = tandem_repeat_bed
    }

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
                    ref_fasta = ref_fasta,
                    ref_fasta_fai = ref_fasta_fai,
                    tandem_repeat_bed = tandem_repeat_bed,
                    chr = contig,
                    prefix = prefix,
                    zones = zones
            }
        }

        call PBSV.Call as pbsv_wg_call {
            input:
                svsigs = pbsv_discover_chr.svsig,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                is_hifi = is_hifi,
                is_ont = is_ont,
                prefix = prefix + ".pbsv",
                zones = zones
        }
    }

    if (!pbsv_discover_per_chr) {
        call PBSV.RunPBSV as PBSVslow {
            input:
                bam = bam,
                bai = bai,
                ref_fasta = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                prefix = prefix,
                tandem_repeat_bed = tandem_repeat_bed,
                is_hifi = is_hifi,
                is_ont = is_ont,
                zones = zones
        }
    }
}
