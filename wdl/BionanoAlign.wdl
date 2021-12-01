version 1.0

import "tasks/Structs.wdl"

workflow CallAlign {
    input {
        File ref_cmap
        File query_cmap
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    call Align {
        input:
            ref_cmap = ref_cmap,
            query_cmap = query_cmap,
            output_prefix = output_prefix,
    }

    output {
        File align_output = Align.mapfiles
    }
}

task Align {
    input {
        File ref_cmap
        File query_cmap
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3 * (ceil(size([ref_cmap], "GB")) + ceil(size([query_cmap], "GB")))

    command <<<
        source activate bionano_minimal
        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        /home/bionano2/tools/pipeline/1.0/RefAligner/1.0/RefAligner -ref ~{ref_cmap} -i ~{query_cmap} \
        -maxthreads $num_core -o ~{output_prefix} -f -stdout -stderr  -endoutlier 1e-2 -outlier 1e-4 \
        -A 5 -M 1 -Mfast 1 -biaswt 0 -T 1e-12 -nosplit 2 -MultiMatches 5

        tar czf output.tar.gz ~{output_prefix}*[cx]map
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/ymostovoy/bionano-pipeline:latest"
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

    output {
        File stdout = output_prefix+".stdout"
        File mapfiles = "output.tar.gz"
    }
}
