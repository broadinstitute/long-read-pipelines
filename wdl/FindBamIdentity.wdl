version 1.0

import "tasks/Utils.wdl"
import "tasks/VariantUtils.wdl"

import "tasks/utils/qc/Fingerprinting.wdl" as FPUtils

workflow FindBamIdentity {

    meta {
        description: "A workflow to identify a flowcell's the true identity, by genotyping it's BAM, against an array of 'truth' genotyped VCF."
    }

    input {
        File aligned_bam
        File aligned_bai
        String expt_type

        Int artificial_baseQ_for_CLR = 10

        String fingerprint_store
        Array[String]? suspect_sample_names

        File ref_map_file
    }

    parameter_meta {
        aligned_bam:        "GCS path to aligned BAM file of the flowcell"
        expt_type:          "There will be special treatment for 'CLR' data (minimum base quality for bases used when computing a fingerprint)"
        artificial_baseQ_for_CLR: "An artificial value for CLR reads used for fingerprint verification (CLR reads come with all 0 base qual)"

        fingerprint_store:      "GS path to where all known fingerprinting GT'ed VCFS are stored"
        suspect_sample_names: "A list of sample names whose FP vcfs to pick (VCFs at FP store are named with these names)"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    Boolean is_clr_bam = expt_type=='CLR'

    call FPUtils.ListGenotypedVCFs { input: fingerprint_store = fingerprint_store }
    call FPUtils.SelectPossibleVCFs { input: fingerprinting_vcf_gs_paths = ListGenotypedVCFs.vcf_gs_paths, suspect_sample_names = suspect_sample_names }

    scatter (vcf in SelectPossibleVCFs.vcfs) {
        call FPUtils.ReheaderFullGRCh38VCFtoNoAlt {input: full_GRCh38_vcf = vcf}

        call FPUtils.FilterGenotypesVCF {
            input:
                fingerprint_vcf = ReheaderFullGRCh38VCFtoNoAlt.reheadered_vcf
        }
        call FPUtils.ExtractGenotypingSites {
            input:
                fingerprint_vcf = FilterGenotypesVCF.ready_to_use_vcf
        }
    }

    call FPUtils.MergeGenotypingSites {input: all_sites = ExtractGenotypingSites.sites}

    call FPUtils.ExtractRelevantGenotypingReads {
        input:
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai,
            genotyping_sites_bed = MergeGenotypingSites.merged_sites
    }

    if (is_clr_bam) {
        call FPUtils.ResetCLRBaseQual {
            input:
                bam = ExtractRelevantGenotypingReads.relevant_reads,
                bai = ExtractRelevantGenotypingReads.relevant_reads_bai,
                arbitrary_bq = artificial_baseQ_for_CLR
        }
    }

    scatter (vcf in FilterGenotypesVCF.ready_to_use_vcf) {

        call VariantUtils.GetVCFSampleName {
            input:
                fingerprint_vcf = vcf
        }

        if ( ! is_clr_bam ) {
            call FPUtils.CheckFingerprint {
                input:
                    aligned_bam     = ExtractRelevantGenotypingReads.relevant_reads,
                    aligned_bai     = ExtractRelevantGenotypingReads.relevant_reads_bai,
                    fingerprint_vcf = vcf,
                    vcf_sample_name = GetVCFSampleName.sample_name,
                    haplotype_map   = ref_map['haplotype_map']
            }
            Float non_clr_lod = CheckFingerprint.metrics_map['LOD_EXPECTED_SAMPLE']
        }

        if (is_clr_bam) {
            call FPUtils.CheckCLRFingerprint {
                input:
                    aligned_bam     = select_first([ResetCLRBaseQual.barbequed_bam]),
                    aligned_bai     = select_first([ResetCLRBaseQual.barbequed_bai]),
                    min_base_q      = artificial_baseQ_for_CLR,
                    fingerprint_vcf = vcf,
                    vcf_sample_name = GetVCFSampleName.sample_name,
                    haplotype_map   = ref_map['haplotype_map']
            }
            Float clr_lod = CheckCLRFingerprint.metrics_map['LOD_EXPECTED_SAMPLE']
        }
    }

    call FindMaxLOD {
        input:
            lod_scores = if is_clr_bam then select_all(clr_lod) else select_all(non_clr_lod)
    }

    if (FindMaxLOD.max_lod < 6) {call Utils.StopWorkflow {input: reason = "No LOD score on the suspected identities are definite." }}

    String matching_vcf = SelectPossibleVCFs.vcfs[FindMaxLOD.idx]

    call GetAllIdentityInfo { input: vcf = matching_vcf }

    output {
        Float lod = FindMaxLOD.max_lod

        String true_id_at_fingerprint_store = GetAllIdentityInfo.resolved_id_at_fingerprint_store
        # String true_smid = GetAllIdentityInfo.resolved_identities[0]
        # String true_collab_sample_id = GetAllIdentityInfo.resolved_identities[1]
        # String true_collab_partic_id = GetAllIdentityInfo.resolved_identities[2]
    }
}

task FindMaxLOD {
    input {
        Array[String] lod_scores
    }

    Int n = length(lod_scores) - 1

    command <<<

        set -eux

        seq 0 ~{n} > indices.txt
        echo "~{sep='\n' lod_scores}" > scores.txt
        paste -d' ' indices.txt scores.txt > to.sort.txt

        sort -k2 -n to.sort.txt > sorted.txt
        cat sorted.txt

        tail -n 1 sorted.txt | awk '{print $1}' > "idx.txt"
        tail -n 1 sorted.txt | awk '{print $2}' > "max_lod.txt"
    >>>

    output {
        Int idx = read_int("idx.txt")
        Float max_lod = read_float("max_lod.txt")
    }

    ###################
    runtime {
        cpu: 2
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task GetAllIdentityInfo {
    meta {
        description: "Get collaborator participant, sample id, and SMID from a fingerprint VCF that follows the naming convension of smid__collabSmId_collabPartId.vcf"
    }
    input {
        String vcf
    }

    String vcf_name = basename(basename(vcf, ".gz"), ".vcf")

    command <<<
        set -eux
        echo ~{vcf_name} | awk -F '.' '{print $1}' > id_at_fingerprint_store.txt
    >>>
    output {
        String resolved_id_at_fingerprint_store = read_string("id_at_fingerprint_store.txt")
    }

    ###################
    runtime {
        cpu: 2
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        docker:"gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
