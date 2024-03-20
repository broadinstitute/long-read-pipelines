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
        
        fastqc_path:        "directory of fastqc_report used for finding BAM files"
        sample_name:        "name of the sample sequenced"
    }

    input {
        RuntimeAttr? runtime_attr_override

        String fastqc_path

        String fasta_path
        String fasta_index_path
        String regions_bed

        String sample_name
        String gcs_out_root_dir
    }

    Int disk_size_gb = 20 + ceil(size(regions_bed, "GB"))

    String results_dir = if (defined(fastqc_path)) then sub(select_first([fastqc_path]), "results\/.*", "") else ""
    String coverage_dir = "~{results_dir}results/SRWholeGenome/~{sample_name}/alignments/"

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/MalariaReports/~{sample_name}/snapshots"
    
    command <<<
        set -euxo
        pwd
        ls

        mkdir -p data
        mkdir -p out
        echo "RETRIEVING REFERENCE FILES..."
        gsutil cp ~{fasta_path} /data
        gsutil cp ~{fasta_index_path} /data

        echo "RETRIEVING BAM AND BAI FILES..."
        echo ~{coverage_dir}
        gsutil ls ~{coverage_dir}  > filelist.txt
        cat filelist.txt | gsutil -m cp -I /data

        echo "CREATING SNAPSHOTS..."        
        bam_file=$(gsutil -m cp -I /data/*.bam -)
        
        python3 make_IGV_snapshots.py \
            -g ~{fasta_path} \
            /data/FP0008-C.bam \
            -r /data/FP0008-C_trim.bed \
            -onlysnap -nf4 -o ~{outdir}
        
    >>>

    output {
        
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
