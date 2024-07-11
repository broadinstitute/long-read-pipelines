version 1.0

workflow VcfdistEval {
    input {
        Array[String] samples
        File truth_vcf
        File eval_vcf
        File bed_file
        File fasta_file
        File fai
        String docker = "timd1/vcfdist:latest"
        Int verbosity = 1
        String region = "chr1:100000000-110000000" 
        String output_directory
    }

    scatter (sample_id in samples)  {

        call SplitVCFbySample as SP_eval { input:
            joint_vcf = eval_vcf,
            reference_index = fai,
            region = region,
            samplename = sample_id
        }
        call SplitVCFbySample as SP_truth { input:
            joint_vcf = truth_vcf,
            reference_index = fai,
            region = region,
            samplename = sample_id,
            
        }

        call RunVcfdistTask {
            input:
                truth_vcf = SP_truth.single_sample_vcf,
                eval_vcf = SP_eval.single_sample_vcf,
                bed_file = bed_file,
                fasta_file = fasta_file,
                docker = docker,
                verbosity = verbosity
        }
    
    }

    call FinalizeToDir as F1 { input: files = RunVcfdistTask.prs_tsv, outdir = output_directory}

    output {
        Array[File] prs_tsv = RunVcfdistTask.prs_tsv
        Array[File] vcfdistsummary = RunVcfdistTask.summary
        Array[File] precrec_tsv = RunVcfdistTask.precrec_tsv
        Array[File] query_tsv = RunVcfdistTask.query_tsv
        Array[File] truth_tsv = RunVcfdistTask.truth_tsv
        
    }
}

task RunVcfdistTask {
    input {
        File truth_vcf
        File eval_vcf
        File bed_file
        File fasta_file
        Int verbosity
        String sampleid
        
        String docker
        Int disk_size_gb = ceil(size(truth_vcf, "GiB") + 10)
        Int mem_gb = 16
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        vcfdist \
            ~{eval_vcf} \
            ~{truth_vcf} \
            ~{fasta_file} \
            -b ~{bed_file} \
            -v ~{verbosity}

        mv "precision-recall-summary.tsv" "~{sampleid}-precision-recall-summary.tsv"
        mv "phasing-summary.tsv" "~{sampleid}-phasing-summary.tsv"
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File prs_tsv = "~{sampleid}-precision-recall-summary.tsv"
        File summary = "summary.vcf"
        File precrec_tsv = "precision-recall.tsv"
        File query_tsv = "query.tsv"
        File truth_tsv = "truth.tsv"
        File phasing_tsv = "~{sampleid}-phasing-summary.tsv"
        File switch_flip_tsv = "switchflips.tsv"
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

        bcftools index ~{joint_vcf}

        bcftools view -s ~{samplename} ~{joint_vcf} -r ~{region} -o ~{samplename}.subset.g.vcf.gz

        bcftools reheader --fai ~{reference_index} ~{samplename}.subset.g.vcf.gz > ~{samplename}.subset.reheadered.g.vcf.gz

        tabix -p vcf ~{samplename}.subset.reheadered.g.vcf.gz

    >>>
    
    output {
		File single_sample_vcf = "~{samplename}.subset.reheadered.g.vcf.gz"
        File single_sample_vcf_tbi = "~{samplename}.subset.reheadered.g.vcf.gz.tbi"
    }


    Int disk_size = 1 + ceil(2 * (size(joint_vcf, "GiB")))

    runtime {
        cpu: 1
        memory: "64 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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