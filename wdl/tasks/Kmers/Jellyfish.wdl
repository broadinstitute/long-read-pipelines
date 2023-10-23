version 1.0

import "../../structs/Structs.wdl"

task KmerCounts {

    meta {
        description: "Generates kmer counts table from single sample."
    }


    input {
        File fasta
        Int kmer_size
        String prefix
        Int hash_size = "100M"
        Int threads = 16
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        RuntimeAttr? runtime_attr_override
    }

    Int file_size = ceil(size(fasta, "GB"))
	Int disk_size = if file_size > 200 then 2*file_size else file_size + 200

    command <<<
        set -euxo pipefail
        jellyfish count -m ~{kmer_size} -s ~{hash_size} -t ~{threads} -C ~{fasta} -o ~{prefix}.jf
        
    >>>

    output {
        File kmercount = "~{prefix}.jf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "hangsuunc/jellyfish:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MergeCounts {
    meta {
        description: "merge kmer counts table from multiple sample."
    }  

    input {
        Array[File] jf
        String outputprefix
        Int threads = 16
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
        RuntimeAttr? runtime_attr_override
    }
    Int file_size = ceil(size(jf, "GB"))
	Int disk_size = if file_size > 200 then 2*file_size else file_size + 200

    command <<<
        set -euxo pipefail
        jellyfish merge -o ~{outputprefix}.merged.jf ~{sep=" " jf}
        
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "hangsuunc/jellyfish:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

}

task findCommonKmers {
    meta {
        description: "find common kmers"
    }  

    input {
        File jf
        String outputprefix
        Int threads = 16
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
        RuntimeAttr? runtime_attr_override
    }
    Int file_size = ceil(size(jf, "GB"))
	Int disk_size = if file_size > 200 then 2*file_size else file_size + 200

    command <<<
        set -euxo pipefail
        jellyfish dump -c ~{jf} 
    >>>
    
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "hangsuunc/jellyfish:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

}