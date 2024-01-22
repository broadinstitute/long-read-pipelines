version 1.0

import "../../../tasks/Utility/Hail.wdl" as Hail
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ConvertToHailMTT2T {

    meta {
        description: "Convert a VCF to a Hail MatrixTable"
    }
    parameter_meta {
    }

    input {
        File bcf
        File referencefa
        File referencefai
        String prefix

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/Hail/~{prefix}"

    call FindFiles{input: sampleFolder = SampleFolder}
    scatter (bcf_file in FindFiles.vcfFiles){
        call preprocess{input: bcf = bcf_file, prefix = basename(bcf_file, ".bcf")}
    }
    
    call Hail.ConvertToHailMT as RunConvertToHailMT {
        input:
            gvcf = preprocess.whole_genome_vcf,
            tbi = preprocess.whole_genome_vcf_tbi
            prefix = prefix,
            outdir = outdir,
            reference = "chm13v2.0",
            ref_fasta = referencefa,
            ref_fai = referencefai,

    }

    ##########
    # store the results into designated bucket
    ##########

    output {
        String joint_mt = RunConvertToHailMT.gcs_path
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

task preprocess {
  input {
    File bcf
    String prefix
  }

command {
    bcftools view -Oz -o ~{prefix}.tmp.vcf.gz ~{bcf}
    # bcftools sort -o whole_genome_sorted.vcf.gz tmp.vcf.gz
    bcftools index --tbi --force ~{prefix}.tmp.vcf.gz


    
  }

  output {
    # Output the list of .bam files
    File vcf = "~{prefix}.tmp.vcf.gz"
    File tbi = "~{prefix}.tmp.vcf.gz.tbi"

  }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
        disks: "local-disk 100 HDD"
    }
}
