version 1.0

import "Structs.wdl"

workflow LRMetrics {
    input {
        File unaligned_bam
        File aligned_bam
        File aligned_bai

        File ref_dict
        File ref_flat
    }

    # metrics on unaligned BAM
    call ReadMetrics as UnalignedReadMetrics {
        input:
            bam = unaligned_bam,
    }

    # metrics on aligned BAM
    call ReadMetrics as AlignedReadMetrics {
        input:
            bam = aligned_bam,
    }

    call MakeChrIntervalList {
        input:
            ref_dict = ref_dict
    }

    scatter (chr_info in MakeChrIntervalList.chrs) {
        call CoverageTrack {
            input:
                bam = aligned_bam,
                bai = aligned_bai,
                chr = chr_info[0],
                start = chr_info[1],
                end = chr_info[2]
        }
    }

    call FlagStats {
        input:
            bam = aligned_bam,
            bai = aligned_bai,
    }

    call RnaSeqMetrics {
        input:
            bam = aligned_bam,
            bai = aligned_bai,
            ref_flat = ref_flat
    }

    output {
        File flag_stats = FlagStats.flag_stats

        Array[File] coverage = CoverageTrack.coverage
        Array[File] coverage_tbi = CoverageTrack.coverage_tbi

        File aligned_np_hist = AlignedReadMetrics.np_hist
        File aligned_range_gap_hist = AlignedReadMetrics.range_gap_hist
        File aligned_zmw_hist = AlignedReadMetrics.zmw_hist
        File aligned_prl_counts = AlignedReadMetrics.prl_counts
        File aligned_prl_hist = AlignedReadMetrics.prl_hist
        File aligned_prl_nx = AlignedReadMetrics.prl_nx
        File aligned_prl_yield_hist = AlignedReadMetrics.prl_yield_hist
        File aligned_rl_counts = AlignedReadMetrics.rl_counts
        File aligned_rl_hist = AlignedReadMetrics.rl_hist
        File aligned_rl_nx = AlignedReadMetrics.rl_nx
        File aligned_rl_yield_hist = AlignedReadMetrics.rl_yield_hist

        File unaligned_np_hist = UnalignedReadMetrics.np_hist
        File unaligned_range_gap_hist = UnalignedReadMetrics.range_gap_hist
        File unaligned_zmw_hist = UnalignedReadMetrics.zmw_hist
        File unaligned_prl_counts = UnalignedReadMetrics.prl_counts
        File unaligned_prl_hist = UnalignedReadMetrics.prl_hist
        File unaligned_prl_nx = UnalignedReadMetrics.prl_nx
        File unaligned_prl_yield_hist = UnalignedReadMetrics.prl_yield_hist
        File unaligned_rl_counts = UnalignedReadMetrics.rl_counts
        File unaligned_rl_hist = UnalignedReadMetrics.rl_hist
        File unaligned_rl_nx = UnalignedReadMetrics.rl_nx
        File unaligned_rl_yield_hist = UnalignedReadMetrics.rl_yield_hist

        File rna_metrics = RnaSeqMetrics.rna_metrics
    }
}

task MakeChrIntervalList {
    input {
        File ref_dict

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    command <<<
        set -euxo pipefail

        grep '^@SQ' ~{ref_dict} | awk '{ print $2 "\t" 1 "\t" $3 }' | sed 's/[SL]N://g' | grep -v -e random -e chrUn -e decoy -e alt -e HLA -e EBV > chrs.txt
    >>>

    output {
        Array[Array[String]] chrs = read_tsv("chrs.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "kgarimella/lr-metrics:0.01.05"
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

task CoverageTrack {
    input {
        File bam
        File bai
        String chr
        String start
        String end

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB"))

    command <<<
        set -euxo pipefail

        samtools depth ~{bam} -r ~{chr}:~{start}-~{end} | bgzip > ~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz
        tabix -p bed ~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz
    >>>

    output {
        File coverage = "~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz"
        File coverage_tbi = "~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "kgarimella/lr-metrics:0.01.05"
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

task FlagStats {
    input {
        File bam
        File bai

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB"))

    command <<<
        set -euxo pipefail

        samtools flagstat ~{bam} > ~{basename}.flag_stats.txt
    >>>

    output {
        File flag_stats = "~{basename}.flag_stats.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "kgarimella/lr-metrics:0.01.05"
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


task RnaSeqMetrics {
    input {
        File bam
        File bai
        File ref_flat

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB") + size(ref_flat, "GB"))

    command <<<
        set -euxo pipefail

        java -jar /picard.jar CollectRnaSeqMetrics \
            I=~{bam} \
            REF_FLAT=~{ref_flat} \
            STRAND=NONE \
            O=~{basename}.rna_metrics.txt
    >>>

    output {
        File rna_metrics = "~{basename}.rna_metrics.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "kgarimella/lr-metrics:0.01.05"
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

task ReadMetrics {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        java -jar /gatk.jar ComputeLongReadMetrics -I ~{bam} -O ~{basename}.read_metrics
    >>>

    output {
        File np_hist = "~{basename}.read_metrics.np_hist.txt"
        File range_gap_hist = "~{basename}.read_metrics.range_gap_hist.txt"
        File zmw_hist = "~{basename}.read_metrics.zmw_hist.txt"
        File prl_counts = "~{basename}.read_metrics.prl_counts.txt"
        File prl_hist = "~{basename}.read_metrics.prl_hist.txt"
        File prl_nx = "~{basename}.read_metrics.prl_nx.txt"
        File prl_yield_hist = "~{basename}.read_metrics.prl_yield_hist.txt"
        File rl_counts = "~{basename}.read_metrics.rl_counts.txt"
        File rl_hist = "~{basename}.read_metrics.rl_hist.txt"
        File rl_nx = "~{basename}.read_metrics.rl_nx.txt"
        File rl_yield_hist = "~{basename}.read_metrics.rl_yield_hist.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             40,
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "kgarimella/lr-metrics:0.01.05"
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

# yield
# read length histogram (https://www.pacb.com/wp-content/uploads/HiFi-Base-Yield-Dens-plot_Sequel-II-System.png)
# read length N50
# number of zmws
# 3000x3000 zmw yield grid

# depth of coverage (autosome)
# depth of coverage (X)
# depth of coverage (Y)
# depth of coverage (MT)
# read lengths
# alignment stats
# yield
# number of passes
# error rate
# error rate / number of passes
# indel error lengths
# contamination
