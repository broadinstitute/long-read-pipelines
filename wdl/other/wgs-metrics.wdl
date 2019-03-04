
workflow WgsMetricsWorkflow {
    File ref_fasta
    String input_bam
    String input_bam_bai

    call WgsMetrics {
        input:
            ref_fasta=ref_fasta,
            input_bam=input_bam,
            input_bam_bai=input_bam_bai
    }
}

task WgsMetrics {
    File ref_fasta
    String input_bam
    String input_bam_bai
    Int num_reads = 100000000
    Int disk_size = ceil(0.02*size(input_bam, "GB") + 10)

    String interval_chrom="1"
    String interval_start="25000000"
    String interval_end  ="120000000"

    #samtools view -h ${input_bam} | head -n ${num_reads} | samtools view -bS - > small.bam \
            #&& samtools index small.bam

    command {
        set -xe

        if echo "${ref_fasta}" | grep -q "hg38"; then
            CHROM="chr${interval_chrom}"
        else
            CHROM="${interval_chrom}"
        fi
        java -Xms2g -jar /gatk.jar PrintReads -R ${ref_fasta} -I ${input_bam} -L "$CHROM":${interval_start}-${interval_end} -O small.bam
        samtools view -H small.bam  > subset.interval_list
        echo "$CHROM	${interval_start}	${interval_end}	+	interval1" >> subset.interval_list
        java -Xmx2G -jar /picard.jar CollectWgsMetrics \
            USE_FAST_ALGORITHM=true \
            INTERVALS=subset.interval_list \
            REFERENCE_SEQUENCE=${ref_fasta} \
            INPUT=small.bam \
            OUTPUT=metrics.output \
            STOP_AFTER=${num_reads} \
            SAMPLE_SIZE=1
        cat metrics.output | grep -v ^# | grep -v ^$ | head -n 2 | verticalize | tee metrics.txt
    }

    output {
        File output_metrics="metrics.txt"
    }

    runtime {
        docker: "weisburd/base-image-for-str-tools@sha256:bc02c67c69bbd13165ef2cd83f7b2ec87814d99fc007070cc7cd8865b29998a9"
        disks: "local-disk ${disk_size} SSD"
        memory: "3 GB"
    }
}
