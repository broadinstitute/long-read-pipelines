version 1.0

import "Structs.wdl"

task RepeatMasker {
    input {
        File fasta
        File famdb

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(fasta, "GB"))+90 # unzipped famdb is 67gb in current release (3.8)
    String famdb_name = basename(famdb)

    command <<<
        set -euxo pipefail

        # unzip famdb if needed, and put in the right folder for RM to find
        mv ~{famdb} RepeatMasker/Libraries/famdb/
        cd RepeatMasker/Libraries/famdb/

        if [[ ~{famdb_name} == *.gz ]]
        then
            gunzip ~{famdb_name}
        fi

        cd ../../../

        ./RepeatMasker/RepeatMasker -e rmblast -pa 4 -s -species human ~{fasta} -rmblast_dir / -trf_prgm /trf
     >>>

     output {
         File RMout = '~{fasta}.out'
     }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "quay.io/ymostovoy/lr-repeatmasker:latest"
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
