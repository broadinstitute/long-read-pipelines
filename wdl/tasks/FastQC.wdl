version 1.0

import "Structs.wdl"

task FastQC {
    input {
        File bam
        File bai

        Int num_cpus = 4

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        export MONITOR_MOUNT_POINT="/cromwell_root"
        wget https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/jts_kvg_sp_malaria/scripts/monitor/legacy/vm_local_monitoring_script.sh -O monitoring_script.sh
        chmod +x monitoring_script.sh
        ./monitoring_script.sh &> resources.log &
        monitoring_pid=$!

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        fastqc -t $num_core --extract ~{bam}

        find . -name 'fastqc_data.txt' -exec mv {} fastqc_data.txt \;
        find . -name 'fastqc_report.html' -exec mv {} fastqc_report.html \;

        number_of_reads=$(grep 'Total Sequences' fastqc_data.txt | awk '{ print $3 }')
        read_length=$(grep 'Sequence length' fastqc_data.txt | awk '{ print $3 }' | cut -f2 -d'-')

        echo $number_of_reads | awk '{ print "number_of_reads\t" $1 }' >> map.txt
        echo $read_length | awk '{ print "read_length\t" $1 }' >> map.txt
        echo $number_of_reads $read_length | awk '{ print "number_of_bases\t" $1*$2 }' >> map.txt

        mean_qual=$(sed -n '/Per base sequence quality/,/END_MODULE/p' fastqc_data.txt | \
            grep -v '^#' | \
            grep -v '>>' | \
            awk '{ print $2 }' | \
            awk '{x+=$1; next} END{print x/NR}')

        echo $mean_qual | awk '{ print "mean_qual\t" $1 }' >> map.txt

        median_qual=$(sed -n '/Per base sequence quality/,/END_MODULE/p' fastqc_data.txt | \
            grep -v '^#' | \
            grep -v '>>' | \
            awk '{ print $2 }' | \
            awk '{x+=$1; next} END{print x/NR}' | \
            sort -n | \
            awk '{ a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')

        echo $median_qual | awk '{ print "median_qual\t" $1 }' >> map.txt

        kill $monitoring_pid
    >>>

    output {
        Map[String, Float] stats_map = read_map("map.txt")

        File stats = "fastqc_data.txt"
        File report = "fastqc_report.html"
        File monitoring_log = "resources.log"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "staphb/fastqc:latest"
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