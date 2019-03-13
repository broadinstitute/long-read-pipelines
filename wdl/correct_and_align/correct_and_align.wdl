workflow CorrectAndAlignWorkflow {
    String input_bam
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    Int num_reads_per_split=200000

    String output_prefix="."
    String base_image="kgarimella/pbtools@sha256:60720374e7ade36d4b9aeab997f688c6b58513ffc31302d5069bbc1810bb8967"

    call ReadLengths {
        input:
            input_bam=input_bam,
            output_lengths=basename(input_bam, ".bam") + ".readlengths.txt",
            docker_image=base_image
    }

    call SplitSubreads {
        input:
            input_bam=input_bam,
            output_prefix=output_prefix,
            num_reads_per_split=num_reads_per_split,
            docker_image=base_image
    }

    scatter (subread_file in SplitSubreads.subread_files) {
        call Minimap2 as Minimap2Uncorrected {
            input:
                ref_fasta=ref_fasta,
                subread_file=subread_file,
                subread_aligned=basename(subread_file, ".bam") + ".aligned.bam",
                docker_image=base_image
        }

        call CCS {
            input:
                subread_file=subread_file,
                subread_ccs=basename(subread_file, ".bam") + ".ccs.bam",
                docker_image=base_image
        }

        call Minimap2 as Minimap2Corrected {
            input:
                ref_fasta=ref_fasta,
                subread_file=CCS.ccs,
                subread_aligned=basename(subread_file, ".bam") + ".ccs.aligned.bam",
                correct_reads=true,
                docker_image=base_image
        }
    }

    call MergeBams as MergeUncorrected {
        input:
            bam_outs=Minimap2Uncorrected.aligned,
            merged_bam=basename(input_bam, ".bam") + ".aligned.merged.bam",
            docker_image=base_image
    }

    call AlignmentStats as AlignmentStatsUncorrected {
        input:
            bam_file=MergeUncorrected.merged,
            ref_fasta=ref_fasta,
            ref_fasta_fai=ref_fasta_fai,
            ref_dict=ref_dict,
            bam_report=basename(MergeUncorrected.merged, ".bam") + ".alignment.report.txt",
            docker_image=base_image
    }

    call MergeBams as MergeCorrected {
        input:
            bam_outs=Minimap2Corrected.aligned,
            merged_bam=basename(input_bam, ".bam") + ".ccs.aligned.merged.bam",
            docker_image=base_image
    }

    call AlignmentStats as AlignmentStatsCorrected {
        input:
            bam_file=MergeCorrected.merged,
            ref_fasta=ref_fasta,
            ref_fasta_fai=ref_fasta_fai,
            ref_dict=ref_dict,
            bam_report=basename(MergeCorrected.merged, ".bam") + ".alignment.report.txt",
            docker_image=base_image
    }
}

task ReadLengths {
    String input_bam
    String output_lengths
    String docker_image

    Int cpus = 1
    Int disk_size = 20

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        samtools view ${input_bam} | awk -F"[\t/]" '{ print $2, $3, length($12) }' > ${output_lengths}

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
    }
}

task SplitSubreads {
    String input_bam
    String output_prefix
    String docker_image
    Int num_reads_per_split

    Int cpus = 1
    Int disk_size = 2*ceil(size(input_bam, "GB"))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        java -jar /gatk.jar SplitSubreadsByZmw -I ${input_bam} -O ${output_prefix} -nr ${num_reads_per_split}

        df -h .
        tree -h
    >>>

    output {
        Array[File] subread_files = glob("*.bam")
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "2G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
    }
}

task Minimap2 {
    File ref_fasta
    File subread_file
    String subread_aligned
    String docker_image

    Boolean? correct_reads
    Boolean correct = select_first([correct_reads, true])
    String correction_arg = if (correct) then " --preset CCS " else ""

    Int cpus = 1
    Int disk_size = ceil(size(ref_fasta, "GB")) + 3*ceil(size(subread_file, "GB"))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        pbmm2 align ${ref_fasta} ${subread_file} ${subread_aligned} --sort ${correction_arg}

        df -h .
        tree -h
    >>>

    output {
        File aligned = "${subread_aligned}"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        disks: "local-disk ${disk_size} SSD"
        memory: "20G"
        bootDiskSizeGb: 20
    }
}

task CCS {
    File subread_file
    String subread_ccs
    String docker_image

    Int cpus = 1
    Int disk_size = 2*ceil(size(subread_file, "GB"))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        ccs --minLength 10000 --maxLength 16000 ${subread_file} ${subread_ccs}

        df -h .
        tree -h
    >>>

    output {
        File ccs = "${subread_ccs}"
        File report = "ccs_report.txt"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "20G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
    }
}

task MergeBams {
    Array[File] bam_outs
    String merged_bam
    String docker_image

    Int cpus = 1
    Int disk_size = 3*ceil(length(bam_outs)*size(bam_outs[0], "GB"))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        java -Xmx4g -jar /gatk.jar MergeSamFiles -I ${sep=" -I " bam_outs} -O ${merged_bam} -AS --CREATE_INDEX

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
    Int disk_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(ref_dict, "GB"))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        java -Xmx4g -jar /gatk.jar CollectAlignmentSummaryMetrics -R ${ref_fasta} -I ${bam_file} -O ${bam_report}

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
    }
}
