version 1.0



workflow MergePhasedVCF{
    meta{
        description: "a workflow that get a bam file from a list of base caller"
    }
    input{
        String SampleFolder
        String Samplename
    }
    call FindFiles{input: sampleFolder = SampleFolder}
    call MergeVcf{input: vcf_files = FindFiles.vcfFiles}
    
    output{
        File wholegenomevcf = MergeVcf.whole_genome_vcf
    }
}


task FindFiles {
  input {
    String sampleFolder
  }

command {
    # Use GSUtil to list all files in the given directory
    gsutil ls "${sampleFolder}" > vcf_files.txt
    # Filter the lines with ".bam" extension and store the result in "bam_files.txt"   
    grep -E "\.bcf$" vcf_files.txt > output.txt

    cat output.txt
  }

  output {
    # Output the list of .bam files
    Array[String] vcfFiles = read_lines("output.txt")

  }

    runtime {
        docker: "broadinstitute/gatk:4.4.0.0"
        disks: "local-disk 100 HDD"
    }
}

task MergeVcf {
  input {
    Array[File] vcf_files
  }

command {
    bcftools concat -o whole_genome.vcf.gz ~{sep=" " vcf_files}
    # bcftools sort -o whole_genome_sorted.vcf.gz whole_genome.vcf.gz
    
  }

  output {
    # Output the list of .bam files
    File whole_genome_vcf = "whole_genome.vcf.gz"

  }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
        disks: "local-disk 100 HDD"
    }
}
