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
        gcs_out_root_dir:   "GCS bucket to store the results"
        aligned_bam:        "GCS path to aligned bam"
        aligned_bai:        "GCS path to aligned bai"
        sample_name:        "name of the sample sequenced"
    }

    input {
        RuntimeAttr? runtime_attr_override

        String aligned_bam
        String aligned_bai

        String fasta_path
        String fasta_index_path
        String regions_bed

        String sample_name
        String gcs_out_root_dir
    }

    Int disk_size_gb = 20 + ceil(size(regions_bed, "GB"))

    String gcs_out_dir = sub(gcs_out_root_dir, "/$", "") + "/MalariaReports/snapshots"
    
    command <<<
        set -euxo
        pwd
        ls

        mkdir -p out
        echo "RETRIEVING BAM AND REF FILES..."
        gsutil cp ~{fasta_path}
        gsutil cp ~{fasta_index_path}
        gsutil cp ~{aligned_bam}

        echo "CREATING SNAPSHOTS..."        
        python3 make_IGV_snapshots.py \
            -g ~{fasta_path} \
            ~{aligned_bam} \
            -r ~{regions_bed} \
            -onlysnap -nf4 -o out
        
        ls out > snapshots.txt

        cat snapshots.txt | gsutil -m cp -I "~{gcs_out_dir}"
    >>>

    Array[File] snapshots = glob("/out/*.png")

    output {
        String gcs_dir = gcs_out_dir
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
