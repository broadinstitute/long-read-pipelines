version 1.0

import "../../structs/Structs.wdl"

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

        num_core=$(awk '/^processor/{print $3}' /proc/cpuinfo | wc -l)

        fastqc -t $num_core --extract ~{bam}

        find . -name 'fastqc_data.txt' -exec mv {} fastqc_data.txt \;
        find . -name 'fastqc_report.html' -exec mv {} fastqc_report.html \;

        number_of_reads=$(grep 'Total Sequences' fastqc_data.txt | awk '{ print $3 }')
        read_length=$(grep 'Sequence length' fastqc_data.txt | awk '{ print $3 }' | cut -f2 -d'-')

        echo $number_of_reads | awk '{ print "number_of_reads\t" $1 }' >> map.txt
        echo $read_length | awk '{ print "read_length\t" $1 }' >> map.txt
        echo $number_of_reads $read_length | awk '{ print "number_of_bases\t" $1*$2 }' >> map.txt

        echo "Calculating number of base qualities:"
        echo -n "Should be: "
        sed -n '/Per base sequence quality/,/END_MODULE/p' fastqc_data.txt | grep -v '^#' | grep -v '>>' | wc -l

        # Need to unset `pipefail` for these commands:
        set -euxo

        # Get the mean and median base quality scores:
        num_base_qualities=$(sed -n '/Per base sequence quality/,/END_MODULE/p' fastqc_data.txt | \
            grep -v '^#' | \
            grep -v '>>' | \
            wc -l)

        echo "Num Base Qualities: ${num_base_qualities}"

        if [[ ${num_base_qualities} -eq 0 ]]; then
            echo "WARNING: No base quality information found in FastQC report.  Setting mean_qual and median_qual to 0."
            mean_qual=0
            median_qual=0
        else
            mean_qual=$(sed -n '/Per base sequence quality/,/END_MODULE/p' fastqc_data.txt | \
                grep -v '^#' | \
                grep -v '>>' | \
                awk '{ print $2 }' | \
                awk '{x+=$1; next} END{print x/NR}')

            median_qual=$(sed -n '/Per base sequence quality/,/END_MODULE/p' fastqc_data.txt | \
                grep -v '^#' | \
                grep -v '>>' | \
                awk '{ print $2 }' | \
                awk '{x+=$1; next} END{print x/NR}' | \
                sort -n | \
                awk '{ a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')
        fi

        echo $mean_qual | awk '{ print "mean_qual\t" $1 }' >> map.txt
        echo $median_qual | awk '{ print "median_qual\t" $1 }' >> map.txt
    >>>

    output {
        Map[String, Float] stats_map = read_map("map.txt")

        File stats = "fastqc_data.txt"
        File report = "fastqc_report.html"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "staphb/fastqc:0.12.1"
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
