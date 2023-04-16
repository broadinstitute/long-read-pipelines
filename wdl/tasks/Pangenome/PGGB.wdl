version 1.0

import "../../structs/Structs.wdl"

task PrepareForPGGB {
    input {
        String output_fname_prefix = "pangenome"
        Array[File] input_genomes

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us-central1-docker.pkg.dev/broad-dsp-lrma/general/pggb:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euxo pipefail

        prepare_pggb.py ~{sep=' ' input_genomes} > "~{output_fname_prefix}.fasta"
        bgzip "~{output_fname_prefix}.fasta"
        samtools faidx "~{output_fname_prefix}.fasta.gz"
    >>>

    output {
        File pangenome_fasta_gz = "~{output_fname_prefix}.fasta.gz"
        File pangenome_fai = "~{output_fname_prefix}.fasta.gz.fai"
    }

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task PGGB {
    input {
        String pangenome_name
        File pangenome_fasta_gz
        File pangenome_fai
        String extra_pggb_args = ""

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        2,
        docker:             "us-central1-docker.pkg.dev/broad-dsp-lrma/general/pggb:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int num_cpu = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    command <<<
        set -euxo pipefail

        # Count number of haplotypes in pangenome
        hap_count=$(zcat ~{pangenome_fasta_gz} | grep -c "^>")

        pggb -i ~{pangenome_fasta_gz} -n "${hap_count}" -o output -t ~{num_cpu} ~{extra_pggb_args}

        mv output/*.final.gfa "~{pangenome_name}.pggb.gfa"
        mv output/*.wfmash.paf "~{pangenome_name}.wfmash.paf"

        mv output/*.final.og "~{pangenome_name}.og"
        mv output/*.final.og.viz_multiqc.png "~{pangenome_name}.odgi_layout1d.png"

        mv output/*.final.og.lay "~{pangenome_name}.og.lay"
        mv output/*.final.og.lay.tsv "~{pangenome_name}.og.lay.tsv"
        mv output/*.final.og.lay.draw.png "~{pangenome_name}.odgi_layout2d.png"
    >>>

    output {
        File pangenome_gfa = "~{pangenome_name}.gfa"
        File pangenome_wfmash_paf = "~{pangenome_name}.wfmash.paf"

        File pangenome_odgi = "~{pangenome_name}.og"
        File pangenome_viz_1d = "~{pangenome_name}.odgi_layout1d.png"

        File pangenome_odgi_layout = "~{pangenome_name}.og.lay"
        File pangenome_odgi_layout_tsv = "~{pangenome_name}.og.lay.tsv"
        File pangenome_viz_2d = "~{pangenome_name}.odgi_layout2d.png"
    }

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}