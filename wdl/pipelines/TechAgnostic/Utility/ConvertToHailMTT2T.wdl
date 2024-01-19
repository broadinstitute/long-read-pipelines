version 1.0

import "../../../tasks/Utility/Hail.wdl" as Hail
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ConvertToHailMTT2T {

    meta {
        description: "Convert a gVCF to a Hail MatrixTable"
    }
    parameter_meta {
        phased_vcf:       "joint-called gVCF file"
        phased_gvcf_tbi:   ".tbi index for joint-called gVCF file"
        prefix:           "prefix for output Hail MatrixTable"
        gcs_out_root_dir: "GCS bucket in which to store the Hail MatrixTable"
    }

    input {
        File bcf_file
        File referencefa
        File referencefai
        String prefix

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/Hail/~{prefix}"

    # Gather across multiple input gVCFs
    call preprocess {
        input:
            bcf = bcf_file
    }
    call Hail.ConvertToHailMT as RunConvertToHailMT {
        input:
            gvcf = preprocess.whole_genome_vcf,
            tbi = preprocess.whole_genome_vcf_tbi,
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
  }

command {
    bcftools view -Oz -o tmp.vcf.gz ~{bcf}
    bcftools sort -o whole_genome_sorted.vcf.gz tmp.vcf.gz
    tabix -p vcf whole_genome_sorted.vcf.gz


    
  }

  output {
    # Output the list of .bam files
    File whole_genome_vcf = "whole_genome_sorted.vcf.gz"
    File whole_genome_vcf_tbi = "whole_genome_sorted.vcf.gz.tbi"

  }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
        disks: "local-disk 100 HDD"
    }
}