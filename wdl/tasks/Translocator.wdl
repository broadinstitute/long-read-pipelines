version 1.0

##########################################################################################
# A workflow that runs Translocator
##########################################################################################

import "Structs.wdl"

task Translocator {
    input {
        File aligned_bam
        File ref_fasta

        String prefix
        Float? min_het_af = 0.2

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        aligned_bam: "sorted input BAM, preferably aligned with NGMLR or minimap2"
        ref_fasta: "reference to which the BAM was aligned"

        prefix: "prefix for output"
        min_het_af: "Minimum het allele frequency (default 0.2)"
    }

    Int disk_size = 2 * ceil(size([aligned_bam, ref_fasta], "GB"))

    command <<<
<<<<<<< HEAD
         translocator -m ~{aligned_bam} -a ~{ref_fasta} -v ~{prefix}.translocator.vcf --min_het_af ~{min_het_af} --global_remap
=======
        translocator -m ~{aligned_bam} -a ~{ref_fasta} -v ~{prefix}.translocator.vcf --min_het_af ~{min_het_af} --global_remap
>>>>>>> 1759fc26314856e93b45933cc6b8d62c03e3d8d2
    >>>

    output {
        File vcf="~{prefix}.translocator.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          3,
        mem_gb:             60,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:   1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-translocator:0.1"
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
