version 1.0

import "../../../tasks/Utility/Hail.wdl" as Hail
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ConvertToHailMTT2T {

    meta {
        description: "Convert a VCF to a Hail MatrixTable"
    }
    parameter_meta {
        whole_genome_vcf:       "VCF file"
        prefix:           "prefix for output Hail MatrixTable"
        gcs_out_root_dir: "GCS bucket in which to store the Hail MatrixTable"
    }

    input {
        File bcf
        File referencefa
        File referencefai
        String prefix

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/Hail/~{prefix}"
    call preprocess{
        input:
            bcf = bcf
    }
    # Gather across multiple input gVCFs
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

task preprocess {
  input {
    File bcf
    String prefix
  }

command {
    bcftools view -Oz -o ~{prefix}.vcf.gz ~{bcf}
    # bcftools sort -o whole_genome_sorted.vcf.gz tmp.vcf.gz
    bcftools index --tbi --force ~{prefix}.vcf.gz


    
  }

  output {
    # Output the list of .bam files
    File whole_genome_vcf = "~{prefix}.vcf.gz"
    File whole_genome_vcf_tbi = "~{prefix}.vcf.gz.tbi"

  }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
        disks: "local-disk 100 HDD"
    }
}