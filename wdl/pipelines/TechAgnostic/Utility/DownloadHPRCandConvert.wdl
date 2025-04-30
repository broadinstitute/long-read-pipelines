version 1.0


# This WDL workflow downloads BAM files from HPRC and converts them to fastq format.
workflow HPRCDownloadandCovert {
    input {
        String sample_id
        Array[String] addresses
    }
    scatter (i in range(length(addresses))) {
        call HPRCDownloadandConvert {
            input:
                sample_id = i,
                address = addresses[i],
        }
    }

    call MergeFastqFile {
        input:
            fastq_files = HPRCDownloadandConvert.fastq_file
    }

    
    output {
        File merged_fastq = MergeFastqFile.fastq_file

    }
}


task HPRCDownloadandConvert {
    input {
        String sample_id
        String address

    }
    parameter_meta {
    }
    
    Int disk_size_gb = 2048
    Int mem_gb = 128
    
    command <<<
        set -euxo pipefail
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        
        BILLING_PROJECT="broad-firecloud-dsde-methods"
        df -h

        # 1. Downloading the files
        if [[ ~{address} == gs://* ]]; then
            SUCCESS=$(gsutil -u ${BILLING_PROJECT} cp ~{address} . && echo 1 || echo 0)
        else
            SUCCESS=$(wget ~{address} && echo 1 || echo 0)
        fi

        if [[ ${SUCCESS} -ne 1 ]]; then
            echo "Download failed. Exiting."
            touch $(basename ~{address})
        fi

        # 2. Converting the files to fastq format
        FILE_NAME=$(basename ~{address})
        if [[ ${FILE_NAME} == *.bam ]]; then
            samtools fastq -@ 8 -n ${FILE_NAME} >> ~{sample_id}.fastq 
        elif [[ ${FILE_NAME} == *.fastq.gz ]]; then
            gunzip -c ${FILE_NAME} >> ~{sample_id}.fastq 
        elif [[ ${FILE_NAME} == *.fastq ]]; then
            cat ${FILE_NAME} >> ~{sample_id}.fastq 
        else
            touch ~{sample_id}.fastq
            echo "File format not recognized. Exiting."
        fi
        
    >>>
    
    output {
        File fastq_file = "~{sample_id}.fastq"
    }

    runtime {
        docker: "fcunial/callset_integration"
        cpu: 32
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}

task MergeFastqFile{
    input {
        Array[File] fastq_files
    }
    parameter_meta {
    }
    Int disk_size_gb = 2048
    Int mem_gb = 128
    command <<<
        set -euxo pipefail
        cat ~{sep=" " fastq_files} > merged.fastq
        gzip merged.fastq
    >>>
    
    output {
        File fastq_file = "merged.fastq.gz"
    }

    runtime {
        docker: "fcunial/callset_integration"
        cpu: 32
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}