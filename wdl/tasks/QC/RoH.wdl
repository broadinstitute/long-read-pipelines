version 1.0

workflow RoH {
    input {
        String sample_id
        File eval_vcf
        File eval_vcf_index
    }
    # call SplitVCFbySample { input:
    #     joint_vcf = eval_vcf,
    #     reference_index = reference_fai,
    #     region = region,
    #     samplename = sample_id
    # }
    call RoH { input:
        vcf = eval_vcf,
        vcf_index = eval_vcf_index,
        samplename = sample_id
    }


     output {
        File result = RoH.roh_result
        
    }
}

task RoH {
    input{       
        File vcf
        File vcf_index
        String samplename
    }
    
    command <<<
        set -x pipefail

        bcftools roh -G30 --AF-dflt 0.4 -s ~{samplename} ~{vcf} | bgzip -c > roh.txt.gz

    >>>
    
    output {
		File roh_result = "roh.txt.gz"

    }


    Int disk_size = 1 + ceil(2 * (size(vcf, "GiB")))

    runtime {
        cpu: 4
        memory: "24 GiB"
        disks: "local-disk " + disk_size + " SSD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 1
        maxRetries: 0
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
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

task SplitVCFbySample {
    input{       
        File joint_vcf
        File reference_index
        String region
        String samplename
    }
    
    command <<<
        set -x pipefail

        time bcftools index ~{joint_vcf}

        time bcftools view --threads 4 -r ~{region} ~{joint_vcf} -o ~{samplename}.subset.g.vcf.gz

        time bcftools reheader --fai ~{reference_index} ~{samplename}.subset.g.vcf.gz > ~{samplename}.subset.reheadered.g.vcf.gz

        tabix -p vcf ~{samplename}.subset.reheadered.g.vcf.gz

    >>>
    
    output {
		File single_sample_vcf = "~{samplename}.subset.reheadered.g.vcf.gz"
        File single_sample_vcf_tbi = "~{samplename}.subset.reheadered.g.vcf.gz.tbi"
    }


    Int disk_size = 1 + ceil(2 * (size(joint_vcf, "GiB")))

    runtime {
        cpu: 4
        memory: "24 GiB"
        disks: "local-disk " + disk_size + " SSD" #"local-disk 100 HDD"
        preemptible: 1
        maxRetries: 0
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}
