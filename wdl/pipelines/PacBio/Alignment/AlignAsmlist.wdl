version 1.0

workflow AlignAsmtoGenome{
    meta{
        description: "a workflow that extract multiple sample in a genomic interval and construct a graph"
    }
    input{
        File reference
        String prefix
        String sampleFolder
        String output_directory
    }

    call FindFiles{input: sampleFolder = sampleFolder}

    scatter (fasta in FindFiles.fastaFiles)  {
        call minimap_align{input: query_fasta=fasta, reference_fasta= reference, pref= basename(fasta)}

    }

    call FinalizeToDir as F1 { input: files = minimap_align.bam, outdir = output_directory}
    call FinalizeToDir as F2 { input: files =  minimap_align.bai,outdir = output_directory }


    
    output{
        # File merged_fa = catfasta.fasta
        Array[File] bam_file = minimap_align.bam
        Array[File] bai_file = minimap_align.bai
    }
}

task FindFiles {
  input {
    String sampleFolder
  }

command {
    # Use GSUtil to list all files in the given directory
    gsutil ls "${sampleFolder}" > files.txt

    # Filter the lines with ".bam" extension and store the result in "bam_files.txt"
    grep -E "\.fasta$" files.txt > fasta_files.txt

    cat fasta_files.txt
  }

  output {
    # Output the list of .bam files
    Array[String] fastaFiles = read_lines("fasta_files.txt")
  }

    runtime {
        docker: "broadinstitute/gatk:4.4.0.0"
        disks: "local-disk 100 HDD"
    }
}



task minimap_align{
    input{
        File query_fasta
        File reference_fasta
        String pref
    }
    command <<<
    minimap2 -ax asm20 ~{reference_fasta} ~{query_fasta} > ~{pref}.aln.sam

    samtools view -H ~{pref}.aln.sam

    samtools sort ~{pref}.aln.sam -o ~{pref}.aln.sorted.bam
    samtools index ~{pref}.aln.sorted.bam ~{pref}.aln.sorted.bam.bai

    >>>

    output{
        File sam = "~{pref}.aln.sam"
        File bam="~{pref}.aln.sorted.bam"
        File bai= "~{pref}.aln.sorted.bam.bai"
    }

    Int disk_size = 100

    runtime {
        cpu: 1
        memory: "128 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
    }
}

task FinalizeToDir {

    meta {
        description: "Copies the given file to the specified bucket."
    }

    parameter_meta {
        files: {
            description: "files to finalize",
            localization_optional: true
        }
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
        outdir: "directory to which files should be uploaded"
    }

    input {
        Array[File] files
        String outdir

        File? keyfile
    }

    String gcs_output_dir = sub(outdir, "/+$", "")

    command <<<
        set -euxo pipefail

        cat ~{write_lines(files)} | gsutil -m cp -I "~{gcs_output_dir}"
    >>>

    output {
        String gcs_dir = gcs_output_dir
    }

    #########################
    runtime {
        cpu: 1
        memory: "1 GiB"
        disks: "local-disk " + 10 + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
}