version 1.0

import "../../structs/Structs.wdl"
task DrugResIGV {

    meta {
        description: "DrugResIGV generates IGV snapshots for the different drug resistance loci for malaria."
    }

    parameter_meta {
        runtime_attr_override: "Override the default runtime attributes"

        regions_bed:        "GCS path to bed file containing drug resistance loci"
        fasta_path:         "GCS path to fasta file"
        fasta_index_path:   "GCS path to fasta.fai file"
        aligned_bam:        "GCS path to aligned bam"
        aligned_bai:        "GCS path to aligned bai"
        sample_name:        "name of the sample sequenced"
    }

    input {
        RuntimeAttr? runtime_attr_override

        File aligned_bam
        File aligned_bai

        File fasta_path
        File fasta_index_path
        File regions_bed

        String sample_name
    }

    Int disk_size_gb = 20 + ceil(size(aligned_bam, "GB"))

    
    command <<<
        set -euxo pipefail
        pwd
        ls

        mkdir -p out

        echo "CREATING SNAPSHOTS..."        
        make_IGV_snapshots.py \
            -g ~{fasta_path} \
            ~{aligned_bam} \
            -r ~{regions_bed} \
            -nf4 -o out
    >>>

    output {
        Array[File] snapshots = glob("out/*.png")
        File batch_script = "out/IGV_snapshots.bat"
    }
    

    # ----------------------------------------------------------- #
    # Runtime Config

    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             16,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "stevekm/igv-snapshot-automator"
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
