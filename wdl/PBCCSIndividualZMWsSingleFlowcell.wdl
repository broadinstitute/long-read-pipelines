version 1.0

##########################################################################################
## A workflow that performs per-read CCS correction on PacBio HiFi reads from a single
## flow cell. The workflow first shards the subreads into clusters, then into files
## containing subreads from a single ZMW.  The ZMW bam is then corrected via CCS.  The
## per-read correction status, read length stats, software run time, and memory use are
## individually generated and later combined into a single report.  The corrected data
## itself is not stored.
##########################################################################################

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF

workflow PBCCSIndividualZMWsSingleFlowcell {
    input {
        String raw_reads_gcs_bucket
        Int num_shards = 5000

        String gcs_out_root_dir
    }

    parameter_meta {
        raw_reads_gcs_bucket: "GCS bucket holding subreads BAMs (and other related files) holding the sequences to be CCS-ed"
        num_shards:           "[default-valued] number of sharded BAMs to create (tune for performance)"
        gcs_out_root_dir:     "GCS bucket to store the metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = raw_reads_gcs_bucket }

    # double scatter: one FC may generate multiple raw BAMs, we perform another layer scatter on each of these BAMs
    scatter (subread_bam in FindBams.subread_bams) {
        # break one raw BAM into fixed number of shards
        File subread_pbi = sub(subread_bam, ".bam$", ".bam.pbi")
        call Utils.ShardLongReads { input: unaligned_bam = subread_bam, unaligned_pbi = subread_pbi, num_shards = num_shards }

        # then perform correction on each of the shard
        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PerReadCCS { input: subreads = subreads }
        }

        call MergeReports {
            input:
                files = PerReadCCS.per_read_report,
                out = basename(subread_bam, ".bam") + ".per_read_report.txt"
        }

        ##########
        # store the results into designated bucket
        ##########

        call FF.FinalizeToDir as FinalizeReport {
            input:
                files = [ MergeReports.read_report ],
                outdir = outdir
        }
    }
}

task PerReadCCS {
    input {
        File subreads

        Int min_passes = 3
        Float min_snr = 2.5
        Int min_length = 10
        Int max_length = 50000
        Float min_rq = 0.99
        Boolean by_strand = false

        Int cpus = 4

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(subreads, "GB"))

    command <<<
        set -x

        NUM_ZMWS=$(samtools view ~{subreads} | awk -F"/" '{ print $2 }' | uniq | wc -l)

        pbindex ~{subreads}

        mkdir shards reports reads times
        python /usr/local/bin/shard_bam.py -p shards/shard -n $NUM_ZMWS ~{subreads}

        for ((i = 0 ; i < $NUM_ZMWS; i++)); do
            bam=shards/shard$i.bam
            ZMW=$(samtools view $bam | head -1 | awk -F"/" '{ print $1 "/" $2 }' | sed 's/\//./g')

            # Run CCS:
            /usr/bin/time --verbose \
                ccs --min-passes ~{min_passes} \
                    --min-snr ~{min_snr} \
                    --min-length ~{min_length} \
                    --max-length ~{max_length} \
                    --min-rq ~{min_rq} \
                    --num-threads ~{cpus} \
                    --report-file reports/ccs_report.$ZMW.txt \
                    ~{if by_strand then "--by-strand" else ""} $bam ccs_unmapped.bam \
                > times/timing.$ZMW.txt 2>&1

            samtools view $bam | awk '{ print $1, length($10) }' > reads/read_info.$ZMW.txt

            rm ccs_unmapped.bam

            CCS_STATUS=$(cat reports/ccs_report.$ZMW.txt | sed 's/ (.*//' | sed 's/ : / /' | sed 's/ *1/\t1/' | sed 's/ *0/\t0/' | sed 's/ /_/g' | grep -v -e Exclusive_counts_for_ZMWs_failing_filters -e ZMWs_fail_filters -e ZMWs_shortcut_filters -e ZMWs_input | sort -n -k2 -r | head -1 | awk '{ print $1 }')
            NUM_SUBREADS=$(cat reads/read_info.$ZMW.txt | wc -l)
            USER_TIME=$(cat times/timing.$ZMW.txt | grep 'User time' | awk -F": " '{ print $2 }')
            MAX_RSS=$(cat times/timing.$ZMW.txt | grep 'Maximum resident set size' | awk -F": " '{ print $2 }')
            RES=$(awk '{ print $2 }' reads/read_info.$ZMW.txt | datamash q1 1 median 1 q3 1 iqr 1 mean 1 sstdev 1 min 1 max 1)
            RLS=$(awk '{ print $2 }' reads/read_info.$ZMW.txt | tr '\n' ',' | sed 's/,$//')

            echo $ZMW $CCS_STATUS $NUM_SUBREADS $USER_TIME $MAX_RSS $RES $RLS >> per_read_report.txt
        done
    >>>

    output {
        File per_read_report = "per_read_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.20"
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

task MergeReports {
    input {
        Array[File] files
        String out

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(files, "GB"))

    command <<<
        set -x

        echo 'zmw ccs_status num_subreads user_time max_rss length_q1 length_median length_q3 length_iqr length_mean length_sd length_min length_max read_lengths' > ~{out}
        cat ~{sep=' ' files} >> ~{out}
    >>>

    output {
        File read_report = out
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.20"
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
