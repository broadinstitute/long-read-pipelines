# download hd5 files and convert to BAM

workflow DownloadPacBioReadsWorkflow {
    Array[String] hdf5_file_urls
    File ref_fasta
    File ref_fasta_fai
    String output_prefix
    String docker_image="weisburd/download-pacbio-study@sha256:96e928a08fef30cf76ccdd1306b6622e9b96dc029b952859ebffd640ecb84308"

    scatter (hdf5_file_url in hdf5_file_urls) {
        call DownloadHdf5File {
            input:
                hdf5_file_url=hdf5_file_url,
                docker_image=docker_image
        }

        call AlignHdf5File {
            input:
                reads_fastq_gz=DownloadHdf5File.reads_fastq_gz,
                ref_fasta=ref_fasta,
                ref_fasta_fai=ref_fasta_fai,
                docker_image=docker_image
        }
    }

    call ComputeTotalBamSize {
        input:
            bam_file_sizes=AlignHdf5File.output_bam_size,
            docker_image=docker_image
    }

    call GatherBams {
        input:
            bam_files=AlignHdf5File.output_bam,
            total_bam_size=ComputeTotalBamSize.total_bam_size,
            ref_fasta=ref_fasta,
            ref_fasta_fai=ref_fasta_fai,
            docker_image=docker_image,
            output_prefix=output_prefix
    }
}

task DownloadHdf5File {
    String hdf5_file_url
    String docker_image
    Int disk_size = 50

    command {
        set -xe

        echo ${hdf5_file_url}

        wget ${hdf5_file_url} -O hdf5.tgz

        tar xzf *.tgz
        rm *.tgz

        python /bash5tools.py --readType subreads --outType fastq *.bas.h5

        df -kh
        ls -lh1
        rm -f *.h5
        mv *.fastq reads.fastq
        gzip reads.fastq
    }

    output {
        File reads_fastq_gz="reads.fastq.gz"
    }

    runtime {
        docker: "${docker_image}"
        disks: "local-disk ${disk_size} HDD"
        preemptible: 3
        maxRetries: 3
    }
}

task AlignHdf5File {
    File reads_fastq_gz
    File ref_fasta
    File ref_fasta_fai
    String docker_image

    Int disk_size = ceil(size(ref_fasta, "GB") + 2 * size(reads_fastq_gz, "GB") + 10)


    command {
        set -xe

        /minimap2 -ax map-pb ${ref_fasta} ${reads_fastq_gz} > aln.sam

        samtools view  -T ${ref_fasta} -b aln.sam > aln.bam
        rm aln.sam
        samtools view -c aln.bam
        samtools sort --reference=${ref_fasta} --output-fmt=BAM -o aln.sorted.bam aln.bam
        samtools view -c aln.sorted.bam
    }

    output {
        File output_bam="aln.sorted.bam"
        Float output_bam_size=size(output_bam, "GB")
    }

    runtime {
        docker: "${docker_image}"
        memory: "16 GB"
        disks: "local-disk ${disk_size} SSD"
        preemptible: 1
    }
}

task ComputeTotalBamSize {
    Array[Float] bam_file_sizes
    String docker_image

    command {
        set -xe

        python -c 'print(int(sum([ ${sep=', ' bam_file_sizes} ])))'
    }

    output {
        Int total_bam_size = read_int(stdout()) + 1
    }

    runtime {
        docker: "${docker_image}"
    }
}

task GatherBams {
    Array[File] bam_files
    Float total_bam_size
    File ref_fasta
    File ref_fasta_fai
    String docker_image
    String output_prefix

    Int disk_size = ceil(3 * total_bam_size + size(ref_fasta, "GB") + 20)

    command {
        set -xe

        java -Xmx15G -jar /picard.jar MergeSamFiles \
            I=${sep=' I=' bam_files} \
            O=${output_prefix}.bam \
            SORT_ORDER=coordinate \
            ASSUME_SORTED=true \
            USE_THREADING=true \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=LENIENT

        mv ${output_prefix}.bai ${output_prefix}.bam.bai
    }

    output {
        File output_bam = "${output_prefix}.bam"
        File output_bam_bai = "${output_prefix}.bam.bai"
    }

    runtime {
        docker: "${docker_image}"
        memory: "16 GB"
        disks: "local-disk ${disk_size} SSD"
    }
}
