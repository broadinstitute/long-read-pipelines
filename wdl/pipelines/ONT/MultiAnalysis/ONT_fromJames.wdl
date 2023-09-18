# a task that takes the ONT summary folder and outputs a string array of all of the bams in the basecalling directory and a file handle for the summary file
version 1.0

workflow ONT_samples{
    meta{
        description: "a workflow that get a bam file from a list of base caller"
    }
    input{
        String SampleBasecallingFolder
        String Samplename
    }
    call FindONTResultsFiles{input: sampleBasecallingFolder = SampleBasecallingFolder}
    call MergeSamFiles{input: inputBams = FindONTResultsFiles.bamFiles}
    call AddOrReplaceReadGroups{input: inputBam=MergeSamFiles.outputBam, samplename=Samplename}

    output{
        File final_bam = AddOrReplaceReadGroups.outputBam
        File final_bai = AddOrReplaceReadGroups.outputBamIndex
    }
}

task FindONTResultsFiles {
  input {
    String sampleBasecallingFolder
  }

command {
    # Use GSUtil to list all files in the given directory
    gsutil ls "${sampleBasecallingFolder}/pass" > basecalled_files.txt
    #gsutil ls "${sampleBasecallingFolder}" > basecallroot_files.txt

    # Filter the lines with ".bam" extension and store the result in "bam_files.txt"
    grep -E "\.bam$" basecalled_files.txt > bam_files.txt

    # Filter the lines with "sequencing_summary_......cnv" extension and store the result in "sequencing_summary_files.txt"
    #grep -E "sequencing_summary.*\.txt" basecallroot_files.txt | head -n 1  > sequencing_summary_files.txt

    cat bam_files.txt
    #cat sequencing_summary_files.txt
  }

  output {
    # Output the list of .bam files
    Array[String] bamFiles = read_lines("bam_files.txt")

    # Output the list of sequencing summary files
    #String ONTSummaryFile = "gs://lrma-targeted-grid-upload/experimentdata/hg002_svs_list1_target_control/HMW_HG002_FGHI_Preshear_Gtube/20230823_1508_1H_PAS12472_db5685df/basecalling/sequencing_summary.txt"
    #String ONTSummaryFile = read_string("sequencing_summary_files.txt")
  }

    runtime {
        docker: "broadinstitute/gatk:4.4.0.0"
        disks: "local-disk 100 HDD"
    }
}


task AddOrReplaceReadGroups {
  input {
    File inputBam
    String samplename
  }

  String name=basename(inputBam, ".bam")
  command {
    gatk AddOrReplaceReadGroups \
      INPUT=${inputBam} \
      OUTPUT=${name}.reheader.bam \
      VALIDATION_STRINGENCY=LENIENT \
      RGID=4 \
      RGLB=lib1 \
      RGPL=ONT \
      RGPU=unit1 \
      RGSM=${samplename}

    gatk BuildBamIndex \
      INPUT=${name}.reheader.bam
  }

   output {
        File outputBam = "~{name}.reheader.bam"
        File outputBamIndex = "~{name}.reheader.bai"
    }

    runtime {
        docker: "broadinstitute/gatk:4.4.0.0"
        disks: "local-disk 100 HDD"
    }
}

task MergeSamFiles {
  input {
    Array[File] inputBams
  }
  String name=basename(inputBams[0], ".bam")
  command {
     echo "~{sep='" >> inputFiles.list\necho "' inputBams}" >> inputFiles.list

     wc -l inputFiles.list

    gatk MergeSamFiles \
      -I inputFiles.list \
      -O ${name}.merged.bam \
      -SO coordinate

    gatk BuildBamIndex \
       INPUT=${name}.merged.bam
  }

  output {
    File outputBam = "~{name}.merged.bam"
    File outputBamIndex = "~{name}.merged.bai"
  }
  runtime {
      docker: "broadinstitute/gatk:4.4.0.0"
      disks: "local-disk 300 HDD"
      memory: "40 GB"
  }
}