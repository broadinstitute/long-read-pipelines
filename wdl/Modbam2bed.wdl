version 1.0

import "tasks/Structs.wdl"

workflow CallModbam2bed {
    input {
        File ref_map_file
        File aligned_mod_bam
        File aligned_mod_bai
        String prefix
        String mod="5mC"
        Int min_reads=6
        File? region_bed
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    call Modbam2bed {
        input:
            ref_fasta = ref_map['fasta'],
            aligned_mod_bam = aligned_mod_bam,
            aligned_mod_bai = aligned_mod_bai,
            prefix=prefix,
            mod = mod,
            min_reads = min_reads,
            region_bed = region_bed
    }

    output {
        File mod_bed = Modbam2bed.mod_bed
        File mod_bedgraph = Modbam2bed.mod_bedgraph
    }
}

task Modbam2bed {
    input {
        File ref_fasta
        File aligned_mod_bam
        File aligned_mod_bai
        String prefix
        String mod
        Int min_reads
        File? region_bed

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(1.5*(size(ref_fasta, "GB") + size(aligned_mod_bam, "GB")))

    command <<<
        set -euxo pipefail

        num_cores=$(grep -c '^processor' /proc/cpuinfo | awk '{ print $1 - 1 }')

        modbam2bed -m ~{mod} -e ~{ref_fasta} ~{aligned_mod_bam} > ~{prefix}.mod_mapped.bed

        if test -f ~{region_bed}; then bedtools intersect -wa -a ~{prefix}.mod_mapped.bed -b ~{region_bed} > tmp; mv tmp ~{prefix}.mod_mapped.bed; fi

        awk '$5>=MIN {OFS="\t"; print $1,$2,$3,$11}' MIN="~{min_reads}" ~{prefix}.mod_mapped.bed > ~{prefix}.mod_mapped.min~{min_reads}.bedgraph
    >>>

    output {
        File mod_bed = "~{prefix}.mod_mapped.bed"
        File mod_bedgraph = "~{prefix}.mod_mapped.min~{min_reads}.bedgraph"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/ymostovoy/modbam2bed:latest"
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