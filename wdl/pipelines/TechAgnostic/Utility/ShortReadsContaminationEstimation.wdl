version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/BAMutils.wdl" as BU
import "../../../tasks/QC/Contamination.wdl"

# This workflow estimates contamination levels for short-read BAM files
struct VBID2_config {
    File genotyping_sites

    Boolean is_hgdp_sites
    Boolean is_100k_sites

    Boolean disable_baq

    Int max_retries

    String? tech
}

workflow ShortReadsContaminationEstimation {
    meta {
        description:
        "Estimate the cross-individual contamination level of a GRCh38 BAM file for short-read data."
    }

    input {
        File bam
        File bai
        String tech  # Short-read technology, e.g., NovaSeq
        File ref_map_file

        File gt_sites_bed
        Boolean is_hgdp_sites
        Boolean is_100k_sites

        Boolean disable_baq

        String disk_type

        Int max_retries
    }

    parameter_meta {
        gt_sites_bed:     "Bed file holding the genotyping sites for short-read data."
        is_hgdp_sites:    "Provided BED is HGDP genotyping sites."
        is_100k_sites:    "Provided BED is 100k genotyping sites, suitable for short-read workflows."
        disable_baq:      "If turned on, BAQ computation will be disabled (faster operation). Recommended for short-read data."

        tech: "Technology used for generating the data; accepted values: [NovaSeq, HiSeq]"

        max_retries: "Number of retries for transient errors (e.g., disk or network issues)."
    }

    # Thresholds adjusted for short-read technologies
    Map[String, Int] bam_threshold_per_tech = {'NovaSeq': 50, 'HiSeq': 100}  # Adjusted for short-read data
    Int bam_file_threshold = bam_threshold_per_tech[tech]

    if (bam_file_threshold > ceil(size(bam, "MiB"))) {
        Float extreme_low_cov_val = 50.0  # Adjusted low-coverage value for short reads
    }

    if (bam_file_threshold <= ceil(size(bam, "MiB"))) {
        # Prepare reference map and proceed with pileup
        Map[String, String] ref_map = read_map(ref_map_file)

        Int scaleup_factor = 20
        call BU.BamToRelevantPileup as Pileup {
            input:
                bam = bam,
                bai = bai,
                bed = gt_sites_bed,
                ref_fasta = ref_map['fasta'],
                disable_baq = disable_baq,
                disk_type = disk_type,
                max_retries = max_retries
        }

        call Contamination.VerifyBamID {
            input: 
                pileup = Pileup.pileups, 
                ref_fasta = ref_map['fasta'], 
                is_hgdp_sites = is_hgdp_sites, 
                is_100k_sites = is_100k_sites
        }
    }

    output {
        Float contamination_est = select_first([VerifyBamID.contamination_est, extreme_low_cov_val])
    }
}
