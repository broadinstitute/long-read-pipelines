version 1.0

import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/QC/Fingerprinting.wdl" as FPUtils
import "../../../tasks/Utility/VariantUtils.wdl"

workflow VerifyFingerprint {

    meta {
        description: "A workflow to check correctness of metadata on a flowcell, by genotyping it's BAM generated with its metadata, against a 'truth' genotyped VCF."
    }
    parameter_meta {
        aligned_bam:        "GCS path to aligned BAM file, supposed to be of the same sample as from the fingerprinting VCF"
        expt_type:          "There will be special treatment for 'CLR' data (minimum base quality for bases used when computing a fingerprint)"
        artificial_baseQ_for_CLR: "An artificial value for CLR reads used for fingerprint verification (CLR reads come with all 0 base qual)"

        fingerprint_store:  "GS path to where all known fingerprinting GT'ed VCFS are stored"
        smid:               "SM- prefixed ID"
        collaborator_smid:  "Collaborator sample ID"
        collaborator_participant_id:    "Collaborator participant ID"

        use_this_fp_vcf:    "Optional gt VCF, if provided, used for fingerprint verification (fingerprint_store, smid, collaborator_smid, collaborator_participant_id will all be ignored)"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        File aligned_bam
        File aligned_bai
        String expt_type

        Int? artificial_baseQ_for_CLR = 10

        String fingerprint_store
        String smid
        String collaborator_smid
        String collaborator_participant_id

        File? use_this_fp_vcf

        File ref_map_file

        String gcs_out_root_dir
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/VerifyFingerprint"

    if (!defined(use_this_fp_vcf)) {
        call FPUtils.ListGenotypedVCFs { input: fingerprint_store = fingerprint_store }
        String separator = "__"
        String intended_file= smid + separator + collaborator_smid + separator + collaborator_participant_id + ".vcf.gz"

        call FPUtils.PickGenotypeVCF {
            input:
                fingerprinting_vcf_gs_paths = ListGenotypedVCFs.vcf_gs_paths,
                vcf_name = intended_file
        }
    }

    File gt_vcf = if (defined(use_this_fp_vcf)) then select_first([use_this_fp_vcf]) else select_first([PickGenotypeVCF.vcfs])[0]
    call FPUtils.ReheaderFullGRCh38VCFtoNoAlt as reheader {input: full_GRCh38_vcf = gt_vcf}

    call VariantUtils.GetVCFSampleName {
        input:
            fingerprint_vcf = reheader.reheadered_vcf
    }

    call FPUtils.FilterGenotypesVCF {
        input:
            fingerprint_vcf = reheader.reheadered_vcf
    }

    call FPUtils.ExtractGenotypingSites {
        input:
            fingerprint_vcf = FilterGenotypesVCF.ready_to_use_vcf
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
                fingerprint_vcf = FilterGenotypesVCF.ready_to_use_vcf,
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
                min_base_q      = select_first([artificial_baseQ_for_CLR]),
                fingerprint_vcf = FilterGenotypesVCF.ready_to_use_vcf,
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
