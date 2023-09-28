version 1.0

workflow CalculateCoverage{
    meta{
        description: "a workflow to calculate coverage in a genomic interval"
    }
    input{
        File wholegenomebam
        File wholegenomebai
        String genomeregion
        String prefix
        
    }
    call calculatecoverage{input: inputbam=wholegenomebam, inputindex=wholegenomebai, region=genomeregion, prefix=prefix}
    output{
        File coveragefile = calculatecoverage.depth
    }
}

task calculatecoverage{
    input{
        File inputbam
        File inputindex
        String region
        String prefix
    }
    command <<<
        # samtools coverage -r ~{region} ~{inputbam} -o ~{prefix}.out
        samtools depth -r ~{region} ~{inputbam} -o ~{prefix}.depth
    >>>

    output{
        #File coverage="~{prefix}.out"
        File depth = "~{prefix}.depth"
    }

    Int disk_size = 1 + ceil(2 * size(inputbam, "GiB"))

    runtime {
        cpu: 2
        memory: "8 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        # bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}