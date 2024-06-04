version 1.0



workflow ConsensusAsm{
    meta{
        description: "a workflow that construct consensus assemblies from shapeit4 scaffold"
    }
    input{
        Array[String] SampleIDs
        File final_sv_scaffold
        File final_snp_scaffold
        File reference_fasta
        String output_directory
    }

    call merge_scaffold{ input:
        input_sv_scaffold = final_sv_scaffold,
        input_snp_scaffold = final_snp_scaffold,
        prefix = "chr6"
    }

    scatter (sampleid in SampleIDs){
        call SplitVCFbySample as SP { input:
            joint_vcf = merge_scaffold.combined_scaffold,
            region = "chr6:28381449-33301940", # MHC region in T2T
            samplename = sampleid
        }
        call ConstructConsensus { input:
            single_sample_vcf = SP.single_sample_vcf,
            reference = reference_fasta,
            samplename = sampleid
        }
    }
    call FinalizeToDir as F1 { input:
        files = ConstructConsensus.consensus_hap1,
        outdir = output_directory
    }
    call FinalizeToDir as F2 { input:
        files = ConstructConsensus.consensus_hap2,
        outdir = output_directory
    }
    output{
        Array[File] Haplotype1 = ConstructConsensus.consensus_hap1
        Array[File] Haplotype2 = ConstructConsensus.consensus_hap2

    }
}

task merge_scaffold{
    input{
        File input_sv_scaffold
        File input_snp_scaffold
        String prefix
    }
    command <<<
        bcftools index ~{input_sv_scaffold}
        bcftools index ~{input_snp_scaffold}
        bcftools concat -a -D ~{input_sv_scaffold} ~{input_snp_scaffold} -o ~{prefix}.combined_scaffold.bcf
        bcftools index ~{prefix}.combined_scaffold.bcf
    >>>
   

    output{
        File combined_scaffold = "~{prefix}.combined_scaffold.bcf"
        File combined_scaffold_index = "~{prefix}.combined_scaffold.bcf.csi"

    }

    Int disk_size = 100 

    runtime {
        cpu: 8
        memory: "32 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}


task SplitVCFbySample {
    input{       
        File joint_vcf
        String region
        String samplename
    }
    
    command <<<
        set -x pipefail

        bcftools index ~{joint_vcf}

        bcftools view -s ~{samplename} ~{joint_vcf} -r ~{region} -o ~{samplename}.subset.g.vcf.bgz

        tabix -p vcf ~{samplename}.subset.g.vcf.bgz

    >>>
    
    output {
		File single_sample_vcf = "~{samplename}.subset.g.vcf.bgz"
        File single_sample_vcf_tbi = "~{samplename}.subset.g.vcf.bgz.tbi"
    }


    Int disk_size = 1 + ceil(2 * (size(joint_vcf, "GiB")))

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task ConstructConsensus {
    input{       
        File single_sample_vcf
        File reference
        String samplename
    }
    
    command <<<
        set -x pipefail
        bcftools index ~{single_sample_vcf}
        bcftools consensus -e 'ALT~"<.*>"' -H 1 -f ~{reference} ~{single_sample_vcf} > ~{samplename}_MHC_hap1.fasta
        bcftools consensus -e 'ALT~"<.*>"' -H 2 -f ~{reference} ~{single_sample_vcf} > ~{samplename}_MHC_hap2.fasta

    >>>
    
    output {
		File consensus_hap1 = "~{samplename}_MHC_hap1.fasta"
        File consensus_hap2 = "~{samplename}_MHC_hap2.fasta"
    }


    Int disk_size = 1 + ceil(2 * (size(single_sample_vcf, "GiB")))

    runtime {
        cpu: 1
        memory: "4 GiB"
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