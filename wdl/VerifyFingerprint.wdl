version 1.0

import "tasks/Structs.wdl"
import "tasks/Finalize.wdl" as FF

import "tasks/utils/qc/Fingerprinting.wdl" as FPUtils
import "tasks/VariantUtils.wdl"

workflow VerifyFingerprint {

    meta {
        description: "A workflow to check correctness of metadata on a flowcell, by genotyping it's BAM generated with its metadata, against a 'truth' genotyped VCF."
    }

    input {
        File aligned_bam
        File aligned_bai
        String expt_type

        Int? artificial_baseQ_for_CLR = 10

        File fingerprint_vcf

        File ref_map_file

        String gcs_out_root_dir
    }

    parameter_meta {
        aligned_bam:        "GCS path to aligned BAM file, supposed to be of the same sample as from the fingerprinting VCF"
        expt_type:          "There will be special treatment for 'CLR' data (minimum base quality for bases used when computing a fingerprint)"
        artificial_baseQ_for_CLR: "An artificial value for CLR reads used for fingerprint verification (CLR reads come with all 0 base qual)"

        fingerprint_vcf:    "Single sample Fingerprint VCF file from local database"


        ref_map_file:       "table indicating reference sequence and auxillary file locations"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/VerifyFingerprint"

    call VariantUtils.GetVCFSampleName {
        input:
            fingerprint_vcf = fingerprint_vcf
    }

    call FPUtils.FilterGenotypesVCF {
        input:
            fingerprint_vcf = fingerprint_vcf
    }

    call FPUtils.ExtractGenotypingSites {
        input:
            fingerprint_vcf = FilterGenotypesVCF.read_to_use_vcf
    }

    call FPUtils.ExtractRelevantGenotypingReads {
        input:
            aligned_bam     = aligned_bam,
            aligned_bai     = aligned_bai,
            genotyping_sites_bed = ExtractGenotypingSites.sites,
    }

    if (expt_type!='CLR') {
        call FPUtils.CheckFingerprint {
            input:
                aligned_bam     = ExtractRelevantGenotypingReads.relevant_reads,
                aligned_bai     = ExtractRelevantGenotypingReads.relevant_reads_bai,
                fingerprint_vcf = FilterGenotypesVCF.read_to_use_vcf,
                vcf_sample_name = GetVCFSampleName.sample_name,
                haplotype_map   = ref_map['haplotype_map']
        }
    }

    if (expt_type=='CLR') {
        call FPUtils.ResetCLRBaseQual {
            input:
                bam = ExtractRelevantGenotypingReads.relevant_reads,
                bai = ExtractRelevantGenotypingReads.relevant_reads_bai,
                arbitrary_bq = select_first([artificial_baseQ_for_CLR])
        }
        call FPUtils.CheckCLRFingerprint {
            input:
                aligned_bam     = ResetCLRBaseQual.barbequed_bam,
                aligned_bai     = ResetCLRBaseQual.barbequed_bai,
                min_base_q      = artificial_baseQ_for_CLR,
                fingerprint_vcf = FilterGenotypesVCF.read_to_use_vcf,
                vcf_sample_name = GetVCFSampleName.sample_name,
                haplotype_map   = ref_map['haplotype_map']
        }
    }

    File summary_metrics = select_first([CheckFingerprint.summary_metrics, CheckCLRFingerprint.summary_metrics])
    File detail_metrics  = select_first([CheckFingerprint.detail_metrics,  CheckCLRFingerprint.detail_metrics])

    call FF.FinalizeToFile as FinalizeFingerprintSummaryMetrics { input: outdir = outdir, file = summary_metrics }
    call FF.FinalizeToFile as FinalizeFingerprintDetailMetrics  { input: outdir = outdir, file = detail_metrics }

    Map[String, String] metrics_map = select_first([CheckFingerprint.metrics_map, CheckCLRFingerprint.metrics_map])

    output {
        Float lod_expected_sample = metrics_map['LOD_EXPECTED_SAMPLE']

        File fingerprint_metrics = FinalizeFingerprintSummaryMetrics.gcs_path
        File fingerprint_details = FinalizeFingerprintDetailMetrics.gcs_path
    }
}
