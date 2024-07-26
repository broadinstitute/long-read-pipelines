version 1.0

import "../../../tasks/VariantCalling/DeepVariant.wdl" as DV
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow PBCCSDVOnInterval {

    meta {
        description: "A workflow that performs single sample variant calling on a given interval in PacBio HiFi reads from one or more flow cells. The workflow merges multiple SMRT cells into a single BAM prior to variant calling."
    }
    parameter_meta {
        aligned_bams:       "GCS path to aligned BAM files"
        aligned_bais:       "GCS path to aligned BAM file indices"

        participant_name:   "name of the participant from whom these samples were obtained"
        region:             "interval to call variants on"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        File aligned_bam
        File aligned_bai

        String participant_name
        String region

        File ref_map_file
        String gcs_out_root_dir

        Int dvp_threads = 32
        Int dvp_memory = 128
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCCSDVOnInterval/~{participant_name}"

    String smalldir = outdir + "/variants/small"

    call DV.DeepVariant {
        input:
            bam = aligned_bam,
            bai = aligned_bai,

            ref_fasta = ref_map['fasta'],
            ref_fasta_fai = ref_map['fai'],

            region = region,

            dv_threads = dvp_threads,
            dv_memory = dvp_memory
    }

    call FF.FinalizeToFile as FinalizeDVVcf { input: outdir = smalldir, file = DeepVariant.VCF }
    call FF.FinalizeToFile as FinalizeDVTbi { input: outdir = smalldir, file = DeepVariant.VCF_tbi }
    call FF.FinalizeToFile as FinalizeDVGVcf { input: outdir = smalldir, file = DeepVariant.gVCF }
    call FF.FinalizeToFile as FinalizeDVGTbi { input: outdir = smalldir, file = DeepVariant.gVCF_tbi }

    output {
        File dvp_interval_vcf = FinalizeDVVcf.gcs_path
        File dvp_interval_tbi = FinalizeDVTbi.gcs_path
        File dvp_interval_g_vcf = FinalizeDVGVcf.gcs_path
        File dvp_interval_g_tbi = FinalizeDVGTbi.gcs_path
    }
}
