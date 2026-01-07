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

    Int max_retries

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
        Float bam_sz
        String tech
        File ref_map_file
        File standard_hg38_sq_header

        File gt_sites_bed
        Boolean is_hgdp_sites
        Boolean is_100k_sites

        Boolean disable_baq

        String disk_type

        Int max_retries
    }

    parameter_meta {
        # input:
        gt_sites_bed:     "Bed file holding the genotyping sites."
        is_hgdp_sites:    "Provided BED is HGDP genotyping sites."
        is_100k_sites:    "Provided BED is 100k genotyping sites, not 10k sites."
        disable_baq:      "If turned on, BAQ computation will be disabled (faster operation)."

        tech: "technology used for generating the data; accepted value: [ONT, Sequel, Revio]"

        max_retries: "Because of the strange samtools failures reading from NAS storage, we should make multiple attempts to get away from the trasient errors. If after the max retries, we still get those failures, this task will fail."
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

        call ReheaderBam { input:
            bam = bam,
            bai = bai,
            standard_hg38_sq_header = standard_hg38_sq_header,
            bam_sz = bam_sz,
            hack_readgroup_header = false
        }
        call BU.BamToRelevantPileup as Pileup {
            input:
                bam = ReheaderBam.reheadered_bam,
                bai = ReheaderBam.reheadered_bai,
                bed = gt_sites_bed,
                ref_fasta = ref_map['fasta'],
                disable_baq = disable_baq,
                disk_type = disk_type,
                max_retries = max_retries
        }

        call Contamination.VerifyBamID {
            input: pileup = Pileup.pileups, ref_fasta = ref_map['fasta'], is_hgdp_sites = is_hgdp_sites, is_100k_sites = is_100k_sites
        }
    }

    output {
        Float contamination_est = select_first([VerifyBamID.contamination_est, extreme_low_cov_val])
    }
}

# re-using a task to standardize BAM @SQ header
task ReheaderBam {
    meta {
        description: "Reheader a BAM file (the @SQ part) to facilitate downstream somalier fingerprinting."
    }
    parameter_meta {
        hack_readgroup_header:
        "If true, and there are no @RG lines in the input BAM header, a fake read group line will be added with sample name from sample_name_hack_hint."
    }
    input {
        File bam
        File bai
        Float bam_sz
        File standard_hg38_sq_header

        Boolean hack_readgroup_header = false
        String? sample_name_hack_hint
    }
    String prefix = basename(bam, ".bam")
    String out_prefix = "~{prefix}.sq_reheadered"

    output {
        File reheadered_bam = "~{out_prefix}.bam"
        File reheadered_bai = "~{out_prefix}.bam.bai"
    }

    Boolean fail = hack_readgroup_header && (!defined(sample_name_hack_hint))
    String use_this_sample_name = if (defined(sample_name_hack_hint)) then sample_name_hack_hint else "Null"
    command <<<
    set -eu
        if ~{fail}; then
            echo "You set hack_readgroup_header to true but did not provide sample_name_hack_hint" && exit 1
        fi

        # swap out the @SQ lines in the header
        samtools view -H ~{bam} \
        > old.header.txt

        start_line=$(grep -n "^@SQ" old.header.txt | head -1 | cut -d: -f1)
        end_line=$(grep -n "^@SQ" old.header.txt | tail -1 | cut -d: -f1)
        head -n $((start_line - 1)) old.header.txt > new.header.txt
        cat ~{standard_hg38_sq_header} >> new.header.txt
        tail -n +$((end_line + 1)) old.header.txt >> new.header.txt

        if ~{hack_readgroup_header}; then
            if grep -q "^@RG" old.header.txt; then
                echo "You asked for hack_readgroup_header but @RG lines exist in the header" && exit 1
            else
                echo ""
                printf "@RG\tID:FAKE_HACK\tSM:%s\tPL:LONGREAD" "~{use_this_sample_name}\n" \
                >> new.header.txt
                tail +1 new.header.txt
            fi
        fi

        diff new.header.txt old.header.txt || true

    set -x
        samtools reheader \
            new.header.txt \
            ~{bam} \
        > "~{out_prefix}.bam"

        samtools index -@1 "~{out_prefix}.bam"
    >>>
    #########################
    Int local_ssd_sz = if (bam_sz > 120) then 1500 else 375
    runtime {
        cpu:            2
        memory:         "8 GiB"
        disks:          "local-disk ~{local_ssd_sz} LOCAL"
        preemptible:    1
        docker:         "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}