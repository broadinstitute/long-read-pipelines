version 1.0

##########################################################################################
# This pipeline calls SVs on an input LR BAM using various known SV algorithms
# that are specifically designed to work with long read data.
# Each individual task/algo. is directly callable, if so desired.
##########################################################################################

import "../../structs/Structs.wdl"

task CuteSV {
    input {
        File bam
        File bai
        File ref_fasta

        String preset

        Boolean report_readid = true
        Int min_support	 = 1
        Int min_size = 30
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:           "input BAM from which to call SVs"
        bai:           "index accompanying the BAM"
        ref_fasta:     "reference to which the BAM was aligned to"

        preset:        "calling preset (CLR, HIFI, or ONT)"

        report_readid: "report supporting read ids for each SV"
        min_support:   "minimum number of reads that support a SV to be reported"
        min_size:      "minimum length of SV to be reported"
        prefix:        "prefix for output"
    }

    Int disk_size = 2 * ceil(size([bam, bai, ref_fasta], "GB"))

    Map[String, String] preset_values = {
        "CLR": "--max_cluster_bias_INS 100  --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200  --diff_ratio_merging_DEL 0.5",
        "CCS": "--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5",
        "ONT": "--max_cluster_bias_INS 100  --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100  --diff_ratio_merging_DEL 0.3",
    }

    command <<<
        set -euxo pipefail

        mkdir work

        SM=$(samtools view -H ~{bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g')

        cuteSV --sample $SM \
               --genotype \
               ~{if report_readid then "--report_readid" else ""} \
               --min_support ~{min_support} \
               --min_size ~{min_size} \
               ~{preset_values[preset]} \
               ~{bam} \
               ~{ref_fasta} \
               out.vcf \
               work

        grep -v -e '##fileDate' -e '##CommandLine' out.vcf > ~{prefix}.cutesv.vcf
    >>>

    output {
        File vcf = "~{prefix}.cutesv.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sv:0.1.8"
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
