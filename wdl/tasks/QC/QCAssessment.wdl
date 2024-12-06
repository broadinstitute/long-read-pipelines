version 1.0

import "../../structs/Structs.wdl"

task AssessQualityMetrics {
    meta {
        description: "Assess quality metrics from mosdepth coverage and callable loci data to determine pass/fail status.  To pass, the fraction of callable bases must be greater than min_callable_fraction and at least min_coverage_threshold_regions fraction of the genome must have a coverage depth greater than min_coverage."
    }

    parameter_meta {
        callable_loci_summary: "Summary file from GATK CallableLoci"
        mosdepth_region_bed: "Region bed file from mosdepth"
        min_coverage: "Minimum required mean coverage depth"
        min_coverage_threshold_regions: "Minimum required fraction of genome that must have a coverage depth greater than min_coverage"
        min_callable_fraction: "Minimum required fraction of genome that is callable"
        prefix: "Prefix for output files"
    }

    input {
        File callable_loci_summary
        File mosdepth_region_bed

        Float min_coverage = 5
        Float min_coverage_threshold_regions = 0.2
        Float min_callable_fraction = 0.50
        
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    command <<<
        set -euxo pipefail

        # Extract the coverage information from the mosdepth region bed file
        echo "~{mosdepth_region_bed}" | grep -q '\.bed.gz$'
        r=$?
        if [ $r -eq 0 ]; then
            cat ~{mosdepth_region_bed} | gunzip > tmp.bed
        else
            cp ~{mosdepth_region_bed} tmp.bed
        fi

        awk '
        BEGIN { count = 0; total = 0 }
        {
            total++
            if ($NF > ~{min_coverage}) {
                count++
            }
        }
        END {
            is_good="false"
            if (count / total > ~{min_coverage_threshold_regions}) {
                is_good="true"
            }
            printf("%d\t%d\t%f\t%s", count, total, count / total, is_good)
        }
        ' tmp.bed > mosdepth_qc_status.txt

        # Store the mosdepth qc status in a variable so we can use it in an if statement below:
        mosdepth_qc_status=$(cat mosdepth_qc_status.txt | awk '{print $NF}')

        # Extract callable fraction from CallableLoci summary
        callable_bases=$(grep CALLABLE ~{callable_loci_summary} | awk '{print $2}')
        total_bases=$(awk '{sum+=$2} END {print sum}' ~{callable_loci_summary})
        callable_frac=$(awk "BEGIN {printf(\"%.8f\", ${callable_bases}/${total_bases})}")
        callable_status=$(awk "BEGIN { if (${callable_frac} > ~{min_callable_fraction}) { print \"true\" } else { print \"false\" }}")

        # Set blank default message:
        message=""

        # Determine pass/fail status
        if $mosdepth_qc_status && $callable_status; then
            qc_status="Pass"
        else
            qc_status="Fail"
            # Now make a message so the user knows why QC failed:
            if ! $mosdepth_qc_status; then
                message="Low region coverage ($(awk '{print $3}' mosdepth_qc_status.txt) <= ~{min_coverage_threshold_regions}), "
            fi
            if ! $callable_status; then
                message="${message}Low callable fraction (${callable_frac}), "
            fi
            message=$(echo $message | sed 's/,$//')
        fi

        echo "${qc_status}" > qc_status.txt
        echo "${message}" > qc_message.txt
    >>>

    output {
        String qc_status = read_string("qc_status.txt")
        String qc_message = read_string("qc_message.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "ubuntu:22.04"
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
