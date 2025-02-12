version 1.0

import "../Utility/Utils.wdl"
import "../Visualization/NanoPlot.wdl" as NP
import "../QC/AlignedMetrics.wdl" as AM


workflow SampleLevelAlignedMetrics {

    meta {
        description: "A utility (sub-)workflow to compute coverage on sample-leve BAM, and optionally over a provided BED file"
    }
    parameter_meta {
        aligned_bam: "Aligned BAM file"
        aligned_bai: "Index for the aligned BAM file"
        ref_fasta: "Reference FASTA file"
        bed_to_compute_coverage: "Optional BED file to compute coverage over"
    }

    input {
        File aligned_bam
        File aligned_bai

        File ref_fasta

        File? bed_to_compute_coverage
    }

    call Utils.ComputeGenomeLength { input: fasta = ref_fasta }
    call NP.NanoPlotFromBam { input: bam = aligned_bam, bai = aligned_bai }

    if (defined(bed_to_compute_coverage)) {
        call AM.MosDepthOverBed {
            input:
                bam = aligned_bam,
                bai = aligned_bai,
                bed = select_first([bed_to_compute_coverage])
        }

        call SummarizeDepthOverWholeBed as cov_over_region {
            input:
                mosdepth_output = MosDepthOverBed.regions
        }
    }

    output {

        File? bed_cov_summary = cov_over_region.cov_summary

        Float aligned_num_reads = NanoPlotFromBam.stats_map['number_of_reads']
        Float aligned_num_bases = NanoPlotFromBam.stats_map['number_of_bases_aligned']
        Float aligned_frac_bases = NanoPlotFromBam.stats_map['fraction_bases_aligned']
        Float aligned_est_fold_cov = NanoPlotFromBam.stats_map['number_of_bases_aligned']/ComputeGenomeLength.length

        Float aligned_read_length_mean = NanoPlotFromBam.stats_map['mean_read_length']
        Float aligned_read_length_median = NanoPlotFromBam.stats_map['median_read_length']
        Float aligned_read_length_stdev = NanoPlotFromBam.stats_map['read_length_stdev']
        Float aligned_read_length_N50 = NanoPlotFromBam.stats_map['n50']

        Float average_identity = NanoPlotFromBam.stats_map['average_identity']
        Float median_identity = NanoPlotFromBam.stats_map['median_identity']

        Map[String, Float] reads_stats = NanoPlotFromBam.stats_map
    }
}

task SummarizeDepthOverWholeBed {

    meta {
        description: "Summarize the output of MosDepth over a BED file"
    }
    parameter_meta {
        mosdepth_output: "Output of MosDepth over a BED file"
    }

    input {
        File mosdepth_output

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(mosdepth_output, "GB"))

    String prefix = sub(basename(mosdepth_output, ".regions.bed.gz"), "out.coverage.", "")

    command <<<
        set -euxo pipefail

        echo 'chr start stop gene cov_mean' | awk 'BEGIN {OFS="\t"} {print}' > ~{prefix}.summary.txt
        zcat ~{mosdepth_output} >> ~{prefix}.summary.txt
    >>>

    output {
        File cov_summary = "~{prefix}.summary.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
