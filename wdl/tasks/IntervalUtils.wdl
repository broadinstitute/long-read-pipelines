version 1.0

import "Structs.wdl"

task GetContigNames {
    input {
        File ref_dict

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(ref_dict, "GB"))

    command <<<
        set -euxo pipefail

        grep '^@SQ' ~{ref_dict} | awk '{ print $2 }' | sed 's/SN://' > contig_names.txt
    >>>

    output {
        Array[String] contig_names = read_lines("contig_names.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
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

task GenerateIntervals {
    input {
        File ref_dict

        String selected_contig

        Int chunk_bp  = 6000000
        Int stride_bp = 2000000
        Int buffer_bp = 0
        Boolean buffer_start = false
        Boolean buffer_end = false

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(ref_dict, "GB"))

    command <<<
        set -euxo pipefail

        python3 <<EOF

        import re

        chunk_bp     = ~{chunk_bp}
        stride_bp    = ~{stride_bp}
        buffer_bp    = ~{buffer_bp}
        buffer_start = ~{if buffer_start then "True" else "False"}
        buffer_end   = ~{if buffer_end then "True" else "False"}

        with open("~{ref_dict}", "r") as rd:
            for line in rd:
                if line.startswith("@SQ"):
                    pieces = line.split("\t")
                    contig = re.sub("SN:", "", pieces[1])
                    length = int(re.sub("LN:", "", pieces[2]))

                    if contig == "~{selected_contig}":
                        with open("intervals.txt", "w") as rw:
                            if length < chunk_bp:
                                rw.write(f'{contig}:1-{length}\n')
                            else:
                                start = 1
                                end = 1
                                while end <= length - 1:
                                    end = min(start + chunk_bp - 1, length)

                                    padded_start = start + buffer_bp
                                    if start == 1 and not buffer_start:
                                        padded_start = start

                                    padded_end = end - buffer_bp
                                    if end == length and not buffer_end:
                                        padded_end = end

                                    rw.write(f'{contig}:{padded_start}-{padded_end}\n')

                                    start += stride_bp

                        break

        EOF
    >>>

    output {
        Array[String] intervals = read_lines("intervals.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
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