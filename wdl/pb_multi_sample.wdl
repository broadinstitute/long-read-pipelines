workflow PBMultiSampleWorkflow {
    Array[String] input_bams
    Array[String] input_remaining_bams
    String name
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    String bam_out
    String bam_remaining_out
    String rename

    String output_prefix="."
    String base_image="kgarimella/pbtools@sha256:304d90d8130fc8a8704eef566a7de7b614e386e1c02964ec7821f9bd0c25a046"

    call MergeBams as MergeCorrected {
        input:
            bam_outs=input_bams,
            merged_bam=name + ".merged.bam",
            docker_image=base_image
    }

    call Reheader as ReheaderCorrected {
        input:
            bam_in=MergeCorrected.merged,
            bam_out=bam_out,
            old_sm=name,
            new_sm=rename,
            docker_image=base_image
    }

    call MergeBams as MergeRemaining {
        input:
            bam_outs=input_remaining_bams,
            merged_bam=name + ".merged.bam",
            docker_image=base_image
    }

    call Reheader as ReheaderRemaining {
        input:
            bam_in=MergeRemaining.merged,
            bam_out=bam_remaining_out,
            old_sm=name,
            new_sm=rename,
            docker_image=base_image
    }

#    call ReadLengths as ReadLengthsCorrected {
#        input:
#            input_bam=MergeCorrected.merged,
#            output_lengths=basename(MergeCorrected.merged, ".bam") + ".readlengths.txt",
#            docker_image=base_image
#    }
#
#    call AlignmentStats as AlignmentStatsUncorrected {
#        input:
#            bam_file=MergeCorrected.merged,
#            ref_fasta=ref_fasta,
#            ref_fasta_fai=ref_fasta_fai,
#            ref_dict=ref_dict,
#            bam_report=basename(MergeCorrected.merged, ".bam") + ".alignment.report.txt",
#            docker_image=base_image
#    }
#
#    call Depth as DepthCorrected {
#        input:
#            input_bam=MergeCorrected.merged,
#            input_bai=MergeCorrected.merged_bai,
#            ref_fasta=ref_fasta,
#            ref_fasta_fai=ref_fasta_fai,
#            docker_image=base_image
#    }
#
#    call HaplotypeCaller {
#        input:
#            input_bam=MergeCorrected.merged,
#            input_bai=MergeCorrected.merged_bai,
#            ref_fasta=ref_fasta,
#            ref_fasta_fai=ref_fasta_fai,
#            ref_dict=ref_dict,
#            calls=basename(MergeCorrected.merged, ".bam") + ".hc.vcf",
#            docker_image=base_image
#    }
#
#    call PBSV {
#        input:
#            input_bam=MergeCorrected.merged,
#            input_bai=MergeCorrected.merged_bai,
#            ref_fasta=ref_fasta,
#            ref_fasta_fai=ref_fasta_fai,
#            calls=basename(MergeCorrected.merged, ".bam") + ".pbsv.vcf",
#            docker_image=base_image
#    }
#
#    call Sniffles {
#        input:
#            input_bam=MergeCorrected.merged,
#            input_bai=MergeCorrected.merged_bai,
#            calls=basename(MergeCorrected.merged, ".bam") + ".sniffles.vcf",
#            docker_image=base_image
#    }
}

task Reheader {
    File bam_in
    String bam_out
    String old_sm
    String new_sm
    String docker_image

    Int cpus = 2
    Int disk_size = ceil(2.1*(size(bam_in, "GB")))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        samtools view -H ${bam_in} | sed "s/Genetech_Pilot_CCS/${new_sm}/" | samtools reheader - ${bam_in} > ${bam_out}
        samtools index ${bam_out}

        df -h .
        tree -h
    >>>

    output {
        File b_out = "${bam_out}"
        File b_idx = "${bam_out}.bai"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "4G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        preemptible: 3
        maxRetries: 3
    }
}

task ReadLengths {
    String input_bam
    String output_lengths
    String docker_image

    Int cpus = 1
    Int disk_size = 2

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        samtools view ${input_bam} | awk -F"[\t/]" '{ print $2, $3, $7, length($12) }' > ${output_lengths}

        df -h .
        tree -h
    >>>

    output {
        File listing = "${output_lengths}"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "1G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        preemptible: 3
        maxRetries: 3
    }
}

task SplitIntervalsByChr {
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    String docker_image

    Int cpus = 1
    Int disk_size = ceil((size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(ref_dict, "GB")) * 1.1)

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        java -jar /gatk.jar SplitIntervals -R ${ref_fasta} -mode INTERVAL_COUNT -O ./ -scatter 100000000

        df -h .
        tree -h
    >>>

    output {
        Array[File] interval_files = glob("*.interval_list")
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "2G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        preemptible: 3
        maxRetries: 3
    }
}

task MergeBams {
    Array[File] bam_outs
    String merged_bam
    String docker_image
    Int disk_inflate_factor = 3

    Int cpus = 2
    Int disk_size = disk_inflate_factor*ceil(length(bam_outs)*size(bam_outs[0], "GB"))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        java -Xmx4g -jar /gatk.jar MergeSamFiles -I ${sep=" -I " bam_outs} -O ${merged_bam} -AS --CREATE_INDEX --USE_JDK_DEFLATER --USE_JDK_INFLATER --VALIDATION_STRINGENCY SILENT

        df -h .
        tree -h
    >>>

    output {
        File merged = "${merged_bam}"
        File merged_bai = basename(merged_bam, ".bam") + ".bai"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "20G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        #preemptible: 3
        maxRetries: 3
    }
}

task AlignmentStats {
    File bam_file
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    String bam_report
    String docker_image

    Int cpus = 1
    Int disk_size = ceil((size(bam_file, "GB") + size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(ref_dict, "GB")) * 1.1)

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        java -Xmx4g -jar /gatk.jar CollectAlignmentSummaryMetrics -R ${ref_fasta} -I ${bam_file} -O ${bam_report} --USE_JDK_DEFLATER --USE_JDK_INFLATER --VALIDATION_STRINGENCY SILENT

        df -h .
        tree -h
    >>>

    output {
        File report = "${bam_report}"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "20G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        preemptible: 3
        maxRetries: 3
    }
}

task Depth {
    File input_bam
    File input_bai
    File ref_fasta
    File ref_fasta_fai
    String docker_image

    Int cpus = 2
    Int disk_size = ceil((size(input_bam, "GB") + size(input_bai, "GB") + size(ref_fasta, "GB") + size(ref_fasta_fai, "GB")) * 1.1)

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        samtools view -H ${input_bam} | grep SQ | cut -f2,3 | sed 's/[SL]N://g' > chrs.txt
        bedtools makewindows -g chrs.txt -s 10000 -w 10000 > windows.bed
        bedtools coverage -a windows.bed -b ${input_bam} -mean > coverage.txt
        bedtools genomecov -ibam ${input_bam} > hist.txt
        bedtools nuc -fi ${ref_fasta} -bed windows.bed | cut -f1-3,5 | grep -v '^#' > gc.bed

        df -h .
        tree -h
    >>>

    output {
        File wgsmetrics = "coverage.txt"
        File gcbed = "gc.bed"
        File hist = "hist.txt"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "20G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        preemptible: 3
        maxRetries: 3
    }
}

task HaplotypeCaller {
    File input_bam
    File input_bai
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    String calls
    String docker_image

    Int cpus = 2
    Int disk_size = ceil((size(input_bam, "GB") + size(input_bai, "GB") + size(ref_fasta, "GB") + size(ref_dict, "GB") + size(ref_fasta_fai, "GB")) * 1.1)

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        java -Xmx16g -jar /gatk.jar HaplotypeCaller -R ${ref_fasta} -I ${input_bam} -O ${calls} --use-jdk-deflater --use-jdk-inflater -VS SILENT

        df -h .
        tree -h
    >>>

    output {
        File calls_vcf = "${calls}"
        File calls_idx = "${calls}.idx"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "20G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        preemptible: 3
        maxRetries: 3
    }
}

task PBSV {
    File input_bam
    File input_bai
    File ref_fasta
    File ref_fasta_fai
    String calls
    String docker_image

    Int cpus = 2
    Int disk_size = ceil((size(input_bam, "GB") + size(input_bai, "GB") + size(ref_fasta, "GB") + size(ref_fasta_fai, "GB")) * 1.1)

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        pbsv discover --log-level INFO ${input_bam} calls.svsig.gz
        pbsv call --log-level INFO --ccs -j ${cpus} ${ref_fasta} calls.svsig.gz ${calls}

        df -h .
        tree -h
    >>>

    output {
        File svsigfile = "calls.svsig.gz"
        File calls_vcf = "${calls}"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "20G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        preemptible: 3
        maxRetries: 3
    }
}

task Sniffles {
    File input_bam
    File input_bai
    String calls
    String docker_image

    Int cpus = 3
    Int disk_size = ceil((size(input_bam, "GB") + size(input_bai, "GB")) * 1.1)

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        sniffles -s 5 -t ${cpus} --genotype --cluster --report_seq --report_read_strands --ccs_reads -m ${input_bam} -v ${calls}

        df -h .
        tree -h
    >>>

    output {
        File calls_vcf = "${calls}"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "20G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        preemptible: 3
        maxRetries: 3
    }
}
