version 1.0

import "../Utility/Utils.wdl"
import "../Visualization/NanoPlot.wdl" as NP


workflow SampleLevelAlignedMetrics {

    meta {
        description: "A utility (sub-)workflow to compute coverage on sample-leve BAM"
    }
    parameter_meta {
        aligned_bam: "Aligned BAM file"
        aligned_bai: "Index for the aligned BAM file"
    }

    input {
        File aligned_bam
        File aligned_bai
    }

    call NP.NanoPlotFromBam { input: bam = aligned_bam, bai = aligned_bai }

    output {
        Float aligned_num_reads = NanoPlotFromBam.stats_map['number_of_reads']
        Float aligned_num_bases = NanoPlotFromBam.stats_map['number_of_bases_aligned']
        Float aligned_frac_bases = NanoPlotFromBam.stats_map['fraction_bases_aligned']
        Float aligned_est_fold_cov = NanoPlotFromBam.stats_map['number_of_bases_aligned']/NanoPlotFromBam.stats_map['genome_length']

        Float aligned_read_length_mean = NanoPlotFromBam.stats_map['mean_read_length']
        Float aligned_read_length_median = NanoPlotFromBam.stats_map['median_read_length']
        Float aligned_read_length_stdev = NanoPlotFromBam.stats_map['read_length_stdev']
        Float aligned_read_length_N50 = NanoPlotFromBam.stats_map['n50']

        Float average_identity = NanoPlotFromBam.stats_map['average_identity']
        Float median_identity = NanoPlotFromBam.stats_map['median_identity']

        Map[String, Float] reads_stats = NanoPlotFromBam.stats_map
    }
}

task MosDepthOverBed {
    meta {
        description: "Compute coverage on sample-leve BAM over a provided BED file"
    }
    parameter_meta {
        bam: "Aligned BAM file"
        bai: "Index for the aligned BAM file"
        bed: "BED file to compute coverage over"
        outputBucket: "Cloud directory where output will be stored"
    }
    input {
        File bam
        File bai
        File bed
        String? outputBucket

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB"))
    String basename = basename(bam, ".bam")
    String bedname = basename(bed, ".bed")
    String prefix = "~{basename}.coverage_over_bed.~{bedname}"

    command <<<
        set -euxo pipefail

        mosdepth -t 4 -b ~{bed} -n -x -Q 1 ~{prefix} ~{bam}

        outFil="~{prefix}.regions.bed"
        echo 'chr	start	stop	gene	cov_mean' > "$outFil"
        zcat "~{prefix}.regions.bed.gz" >> "$outFil"

        if ~{defined(outputBucket)}; then
            outDir=$(echo "~{outputBucket}" | sed 's+/$++')
            gcloud storage cp "$outFil" "$outDir/$outFil"
            outFil="$outDir/$outFil"
        fi
    >>>

    output {
        File regions = "$outFil"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-mosdepth:0.3.1"
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
