version 1.0

workflow pbCpG{
    meta{
        description: "a workflow that using pbCpG to extract methylated signal from MM/ML tagged bams"
    }

    input{
        File bam
        File bai
        String output_prefix
        Int nthreads
    }
    call AlignedBamToCpGScores{input: bam=bam, bai=bai, outputprefix = output_prefix, num_threads = nthreads}
    
    output{
        Array[File] output_file = AlignedBamToCpGScores.output_files
    }
}

task AlignedBamToCpGScores{
    input{
        File bam
        File bai
        String outputprefix
        Int num_threads
    }

    String model_dir = "/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite"

    command <<<
        aligned_bam_to_cpg_scores \
            --bam ~{bam} \
            --output-prefix ~{outputprefix} \
            --model ~{model_dir} \
            --threads ~{num_threads}
    >>>

    output{
        Array[File] output_files = glob("*")
    }

    Int disk_size = 100 + ceil(2 * (size(bam, "GiB")))

    runtime {
        cpu: 4
        memory: "16 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "hangsuunc/pb-cpg:v1"
    }
}