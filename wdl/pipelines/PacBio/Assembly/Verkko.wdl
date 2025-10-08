version 1.0

import "../../../structs/Structs.wdl"

workflow Verkko {

    meta {
        description: "Perform a genome assembly using Verkko2."
    }
    parameter_meta {
        reads:    "reads (in fasta or fastq format, compressed or uncompressed)"
        prefix:   "prefix to apply to assembly output filenames"
    }
}

task VerkoAssemble {
    input {
        File pacbio_hifi_reads
        File? nanopore_scaffolding_reads
        String prefix

        Boolean is_haploid = true

        File? hap_kmers_file_1
        File? hap_kmers_file_2
        String? hap_kmers_type

        String extra_args = ""

        Int default_preemptible_tries = 0
        Int default_max_retries = 0
        RuntimeAttr? runtime_attr_override
    }

    String out_folder_name="~{prefix}_assembly_verkko"

    String nanopore_scaffolding_arg = if defined(nanopore_scaffolding_reads) then "--nano ~{nanopore_scaffolding_reads}" else ""
    String hap_kmers_arg = if (defined(hap_kmers_file_1) && defined(hap_kmers_file_2) && defined(hap_kmers_type)) then "--hap-kmers ~{hap_kmers_file_1} ~{hap_kmers_file_2} ~{hap_kmers_type}" else ""

    Int disk_size = 10 + 2*(2 * ceil(size(pacbio_hifi_reads, "GB")) + 2 * ceil(size(nanopore_scaffolding_reads, "GB")) + 2 * ceil(size(hap_kmers_file_1, "GB")) + 2 * ceil(size(hap_kmers_file_2, "GB")))

    command <<<
        set -euxo pipefail

        time verkko \
        -d ~{out_folder_name} \
        --hifi ~{pacbio_hifi_reads} \
        ~{nanopore_scaffolding_arg} \
        ~{true="--haploid" false="" is_haploid} \
        ~{hap_kmers_arg}

        tar -czf ~{out_folder_name}.tar.gz ~{out_folder_name}
    >>>

    output {
        File assembly_fasta = "~{out_folder_name}/~{prefix}.fasta"
        File output_tar = "~{prefix}_assembly_verkko.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          12,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       15,
        preemptible_tries:  default_preemptible_tries,
        max_retries:        default_max_retries,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-verkko:2.2.1"
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