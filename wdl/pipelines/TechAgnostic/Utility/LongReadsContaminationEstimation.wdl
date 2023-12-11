version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/BAMutils.wdl" as BU
import "../../../tasks/QC/Contamination.wdl"

# this is a model that other sub-workflows can potentially follow,
# i.e. define a custom struct so that super workflows can use pre-defined JSON files
struct VBID2_config {
    File genotyping_sites

    Boolean is_hgdp_sites
    Boolean is_100k_sites

    Boolean disable_baq

    String? tech
}

workflow LongReadsContaminationEstimation {
    meta {
        desciption:
        "Estimate the cross-individual contamination level of a GRCh38 bam."
    }

    input {
        File bam
        File bai
        String tech
        File ref_map_file

        File gt_sites_bed
        Boolean is_hgdp_sites
        Boolean is_100k_sites

        Boolean disable_baq

        String disk_type
    }

    parameter_meta {
        # input:
        gt_sites_bed:     "Bed file holding the genotyping sites."
        is_hgdp_sites:    "Provided BED is HGDP genotyping sites."
        is_100k_sites:    "Provided BED is 100k genotyping sites, not 10k sites."
        disable_baq:      "If turned on, BAQ computation will be disabled (faster operation)."

        tech: "technology used for generating the data; accepted value: [ONT, Sequel, Revio]"
    }

    # if the coverage is too low, the tool errors out (and the data won't bring much value anyway)
    # here we guard against it by using bam file size, with a emperically-determined threshold
    Map[String, Int] bam_threshold_per_tech = {'ONT': 450, 'Revio': 150, 'Sequel': 250} # this value is technology dependent
    Int bam_file_threshold = bam_threshold_per_tech[tech]

    if (bam_file_threshold > ceil(size(bam, "MiB"))) {
        Float extreme_low_cov_val = 200  # completely arbitrary
    }

    if (bam_file_threshold <= ceil(size(bam, "MiB"))) {
        # quickly change to pileup
        Map[String, String] ref_map = read_map(ref_map_file)

        Int scaleup_factor = 20
        call BU.BamToRelevantPileup as Pileup {
            input:
                bam = bam,
                bai = bai,
                bed = gt_sites_bed,
                ref_fasta = ref_map['fasta'],
                disable_baq = disable_baq,
                disk_type = disk_type
        }

        call Contamination.VerifyBamID {
            input: pileup = Pileup.pileups, ref_fasta = ref_map['fasta'], is_hgdp_sites = is_hgdp_sites, is_100k_sites = is_100k_sites
        }
    }

    output {
        Float contamination_est = select_first([VerifyBamID.contamination_est, extreme_low_cov_val])
    }
}
