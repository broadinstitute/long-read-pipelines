version 1.0

#######################################################
# This pipeline calls small variants using DeepVariant.
#######################################################

import "../../structs/Structs.wdl"
import "../Utility/VariantUtils.wdl"

workflow Run {
    meta {
        desciption:
        "Runs Clair3 on the input (sharded) BAM."
    }
    parameter_meta {
        how_to_shard_wg_for_calling: "An array of the BAM's shard; each element is assumed to be a tuple of (ID for the shard, (BAM of the shard, BAI of the shard))"
        prefix: "Prefix for output files"
    }

    input {
        Array[Pair[String, Pair[File, File]]] how_to_shard_wg_for_calling
        Boolean is_ont
        String prefix
        Map[String, String] ref_map

        # optimization
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    output {
        File clair_vcf = MergeAndSortClairVCFs.vcf
        File clair_tbi = MergeAndSortClairVCFs.tbi
        File clair_gvcf = MergeAndSortClair_gVCFs.vcf
        File clair_gtbi = MergeAndSortClair_gVCFs.tbi
    }

    #################################################################
    scatter (triplet in how_to_shard_wg_for_calling) {
        call Clair {
            input:
                bam = triplet.right.left,
                bai = triplet.right.right,

                ref_fasta     = ref_map['fasta'],
                ref_fasta_fai = ref_map['fai'],

                preset = if is_ont then "ONT" else "CCS",
                zones = zones
        }
    }
    call VariantUtils.MergeAndSortVCFs as MergeAndSortClairVCFs {
        input:
            vcfs = Clair.vcf,
            ref_fasta_fai = ref_map['fai'],
            prefix = prefix + ".clair"
    }
    call VariantUtils.MergeAndSortVCFs as MergeAndSortClair_gVCFs {
        input:
            vcfs = Clair.gvcf,
            ref_fasta_fai = ref_map['fai'],
            prefix = prefix + ".clair.g"
    }
}

task Clair {

    meta {
        description: "Call variants using Clair3."
    }

    parameter_meta {
        bam:             "input BAM from which to call variants"
        bai:             "index accompanying the BAM"

        ref_fasta:       "reference to which the BAM was aligned to"
        ref_fasta_fai:   "index accompanying the reference"

        sites_vcf:       "sites VCF"
        sites_vcf_tbi:   "sites VCF index"

        chr:             "chr on which to call variants"
        preset:          "calling preset (CCS, ONT)"

        zones:           "zones to run the task on"
        runtime_attr_override: "override the default runtime attributes"
    }

    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        File? sites_vcf
        File? sites_vcf_tbi

        String? chr
        String preset

        String zones

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(select_all([bam, bai, ref_fasta, ref_fasta_fai, sites_vcf]), "GB"))
    String platform = if preset == "CCS" then "hifi" else "ont"

    command <<<
        # avoid the infamous pipefail 141 https://stackoverflow.com/questions/19120263
        set -eux
        SM=$(samtools view -H ~{bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g')

        # example from https://github.com/HKU-BAL/Clair3#option-1--docker-pre-built-image
        set -euxo pipefail

        touch ~{bai}
        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        # --include_all_ctgs is turned on, as scatter-gather chops bam before Clair
        /opt/bin/run_clair3.sh ~{true='--vcf_fn=' false='' defined(sites_vcf)}~{select_first([sites_vcf, ""])} \
            --bam_fn=~{bam} \
            --ref_fn=~{ref_fasta} \
            --threads=${num_core} \
            --platform=~{platform} \
            --model_path="/opt/models/~{platform}" \
            --sample_name=$SM --gvcf ~{true='--ctg_name=' false='' defined(chr)}~{select_first([chr, "--include_all_ctgs"])} \
            --output="./"

        # for chrM, Clair3 creates a header only vcf, copy it to gVCF as-is
        if [[ ! -f merge_output.gvcf.gz ]]; then cp "merge_output.vcf.gz" "merge_output.gvcf.gz"; fi
    >>>

    output {
        File? pileup_vcf = "pileup.vcf.gz"
        File? pileup_vcf_tbi = "pileup.vcf.gz.tbi"
        File? full_alignment_vcf = "full_alignment.vcf.gz"
        File? full_alignment_tbi = "full_alignment.vcf.gz.tbi"

        # save both VCF and gVCF
        File vcf = "merge_output.vcf.gz"
        File? vcf_tbi = "merge_output.vcf.gz.tbi"
        File gvcf = "merge_output.gvcf.gz"
        File? gvcf_tbi = "merge_output.gvcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          36,
        mem_gb:             72,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/clair3:v1.0.4"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
