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

        Array[File] suspect_fingerprint_vcfs

        File ref_map_file
    }

    parameter_meta {
        aligned_bam:        "GCS path to aligned BAM file of the flowcell"
        expt_type:          "There will be special treatment for 'CLR' data (minimum base quality for bases used when computing a fingerprint)"
        artificial_baseQ_for_CLR: "An artificial value for CLR reads used for fingerprint verification (CLR reads come with all 0 base qual)"

        suspect_fingerprint_vcfs:    "An array of single sample fingerprinting VCF files, one of which is suspected to hold the true identity of the BAM"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    Boolean is_clr_bam = expt_type=='CLR'

    scatter (vcf in suspect_fingerprint_vcfs) {
        call FPUtils.ExtractGenotypingSites {
            input:
                fingerprint_vcf = vcf
        }
    }

    call FPUtils.MergeGenotypingSites {input: all_sites = ExtractGenotypingSites.sites}    

    # not as efficient as it can be given genotyping sites will be further filtered down, but this involves less code.
    call FPUtils.ExtractRelevantGenotypingReads {
        input:
            aligned_bam = aligned_bam, 
            aligned_bai = aligned_bai,
            genotyping_sites_bed = MergeGenotypingSites.merged_sites
    }
    
    scatter (vcf in suspect_fingerprint_vcfs) {
        call FPUtils.FilterGenotypesVCF {
            input:
                fingerprint_vcf = vcf
        }

        call VariantUtils.GetVCFSampleName {
            input:
                fingerprint_vcf = vcf
        }

        if ( ! is_clr_bam ) {
            call FPUtils.CheckFingerprint {
                input:
                    aligned_bam     = ExtractRelevantGenotypingReads.relevant_reads,
                    aligned_bai     = ExtractRelevantGenotypingReads.relevant_reads_bai,
                    fingerprint_vcf = FilterGenotypesVCF.read_to_use_vcf,
                    vcf_sample_name = GetVCFSampleName.sample_name,
                    haplotype_map   = ref_map['haplotype_map']
            }
            Float non_clr_lod = CheckFingerprint.metrics_map['LOD_EXPECTED_SAMPLE']
        }

        if (is_clr_bam) {
            call FPUtils.ResetCLRBaseQual {
                input:
                    bam = ExtractRelevantGenotypingReads.relevant_reads,
                    bai = ExtractRelevantGenotypingReads.relevant_reads_bai,
                    arbitrary_bq = artificial_baseQ_for_CLR
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
            Float clr_lod = CheckCLRFingerprint.metrics_map['LOD_EXPECTED_SAMPLE']
        }
    }

    call FindMaxLOD {
        input:
            lod_scores = if is_clr_bam then select_all(clr_lod) else select_all(non_clr_lod)
    }

    if (FindMaxLOD.max_lod < 6) {call Utils.StopWorkflow {input: reason = "No LOD score on the suspected identities are definite." }}

    Array[String] TargetSampleNames = GetVCFSampleName.sample_name
    String sample_name = TargetSampleNames[FindMaxLOD.idx]

    output {
        String true_identity = sample_name
        Float lod = FindMaxLOD.max_lod
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
        docker:"ubuntu:20.04"
    }
}
