workflow CorrectAndAlignWorkflow {
    String input_bam
    String sm
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    File trf
    Int num_reads_per_split=200000

    String output_prefix="."
    String base_image="kgarimella/pbtools@sha256:ae9451dd3e7857c999d6477893eb573da4cb4ab1c2f3e21a178d585c41a3aad7"

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
                sample_name=sm,
                docker_image=base_image
        }

        call CCS {
            input:
                subread_file=subread_file,
                subread_ccs=basename(subread_file, ".bam") + ".ccs.corrected.bam",
                docker_image=base_image
        }

        call RecoverUncorrectedReads {
            input:
                subread_file=subread_file,
                subread_ccs=CCS.ccs,
                subread_remaining=basename(subread_file, ".bam") + ".ccs.uncorrected.bam",
                docker_image=base_image
        }

        call Minimap2 as Minimap2Corrected {
            input:
                ref_fasta=ref_fasta,
                subread_file=CCS.ccs,
                subread_aligned=basename(subread_file, ".bam") + ".ccs.corrected.aligned.bam",
                sample_name=sm,
                correct_reads=true,
                docker_image=base_image
        }

        call Minimap2 as Minimap2Remaining {
            input:
                ref_fasta=ref_fasta,
                subread_file=RecoverUncorrectedReads.remaining,
                subread_aligned=basename(subread_file, ".bam") + ".ccs.uncorrected.aligned.bam",
                sample_name=sm,
                docker_image=base_image
        }
    }

    # uncorrected
    call MergeBams as MergeUncorrected {
        input:
            bam_outs=Minimap2Uncorrected.aligned,
            merged_bam=basename(input_bam, ".bam") + ".aligned.merged.bam",
            docker_image=base_image
    }

    call ReadLengths as ReadLengthsUncorrected {
        input:
            input_bam=MergeUncorrected.merged,
            output_lengths=basename(MergeUncorrected.merged, ".bam") + ".readlengths.txt",
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

    # corrected
    call MergeCCSReports {
        input:
            reports=CCS.report,
            docker_image=base_image
    }

    call MergeBams as MergeCorrected {
        input:
            bam_outs=Minimap2Corrected.aligned,
            merged_bam=basename(input_bam, ".bam") + ".ccs.aligned.merged.bam",
            docker_image=base_image
    }

    call ReadLengths as ReadLengthsCorrected {
        input:
            input_bam=MergeCorrected.merged,
            output_lengths=basename(MergeCorrected.merged, ".bam") + ".readlengths.txt",
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

    # remaining
    call MergeBams as MergeRemaining {
        input:
            bam_outs=Minimap2Remaining.aligned,
            merged_bam=basename(input_bam, ".bam") + ".ccs.uncorrected.aligned.merged.bam",
            disk_inflate_factor=10,
            docker_image=base_image
    }

    call ReadLengths as ReadLengthsRemaining {
        input:
            input_bam=MergeRemaining.merged,
            output_lengths=basename(MergeRemaining.merged, ".bam") + ".readlengths.txt",
            docker_image=base_image
    }

    call AlignmentStats as AlignmentStatsRemaining {
        input:
            bam_file=MergeRemaining.merged,
            ref_fasta=ref_fasta,
            ref_fasta_fai=ref_fasta_fai,
            ref_dict=ref_dict,
            bam_report=basename(MergeRemaining.merged, ".bam") + ".alignment.report.txt",
            docker_image=base_image
    }

    call Depth as DepthUncorrected {
        input:
            input_bam=MergeUncorrected.merged,
            input_bai=MergeUncorrected.merged_bai,
            docker_image=base_image
    }

    call Depth as DepthCorrected {
        input:
            input_bam=MergeCorrected.merged,
            input_bai=MergeCorrected.merged_bai,
            docker_image=base_image
    }

    call Depth as DepthRemaining {
        input:
            input_bam=MergeRemaining.merged,
            input_bai=MergeRemaining.merged_bai,
            docker_image=base_image
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
        preemptible: 1
    }
}

task SplitIntervalsByChr {
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    String docker_image

    Int cpus = 1
    Int disk_size = ceil(size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(ref_dict, "GB"))

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
        preemptible: 1
    }
}

task SplitSubreads {
    File input_bam
    String output_prefix
    String docker_image
    Int num_reads_per_split

    Int cpus = 1
    Int disk_size = 4*ceil(size(input_bam, "GB"))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        java -Dsamjdk.compression_level=0 -jar /gatk.jar SplitSubreadsByZmw -I ${input_bam} -O ${output_prefix} -nr ${num_reads_per_split} --use-jdk-deflater --use-jdk-inflater

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
        preemptible: 1
    }
}

task Minimap2 {
    File ref_fasta
    File subread_file
    String subread_aligned
    String sample_name
    String docker_image

    Boolean? correct_reads
    Boolean correct = select_first([correct_reads, false])
    String correction_arg = if (correct) then "asm10" else "map-pb"
    String read_type = if (correct) then "CCS" else "SUBREAD"

    Int cpus = 3
    Int disk_size = ceil(size(ref_fasta, "GB")) + 3*ceil(size(subread_file, "GB"))

    command <<<
        set -uxo pipefail
        df -h .
        tree -h

        ((samtools view -H ${subread_file} | grep '^@RG' | sed 's/\t/\n/g' | grep -v SM) && echo -n 'SM:${sample_name}') | tr "\n" "\t" | sed 's/\t/\\t/g' | sed 's/SUBREAD/${read_type}/' > rg.txt
        echo $?
        samtools view ${subread_file} | head -1 | sed 's/\t/\n/g' | grep -v '^RG' | grep -v ':B:' | grep '^.\+:.:' | cut -d':' -f1 | head -c -1 | tr '\n' ',' > tags.txt
        echo $?
        samtools fastq -T `cat tags.txt` ${subread_file} | minimap2 -R `cat rg.txt` -ax ${correction_arg} -t${cpus} -y ${ref_fasta} - | samtools sort -m4G -o ${subread_aligned} -
        echo $?

        df -h .
        tree -h
    >>>

    output {
        File aligned = "${subread_aligned}"
        File rg = "rg.txt"
        File tags = "tags.txt"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        disks: "local-disk ${disk_size} SSD"
        memory: "20G"
        bootDiskSizeGb: 20
        preemptible: 1
    }
}

task CCS {
    File subread_file
    String subread_ccs
    String docker_image

    Int max_length = 31000
    Int min_passes = 3
    Int cpus = 4
    Int disk_size = 2*ceil(size(subread_file, "GB"))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        ccs --maxLength ${max_length} --minPasses ${min_passes} -j ${cpus} --richQVs ${subread_file} ${subread_ccs}

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
        memory: "40G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        preemptible: 1
    }
}

task RecoverUncorrectedReads {
    File subread_file
    File subread_ccs
    String subread_remaining
    String docker_image

    Int cpus = 2
    Int disk_size = 2*ceil(size(subread_file, "GB"))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        java -Dsamjdk.compression_level=0 -Xmx4g -jar /gatk.jar RecoverUncorrectedReads -I ${subread_file} -C ${subread_ccs} -O ${subread_remaining} --use-jdk-deflater --use-jdk-inflater

        df -h .
        tree -h
    >>>

    output {
        File remaining = "${subread_remaining}"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "20G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        preemptible: 1
    }
}

task MergeCCSReports {
    Array[File] reports
    String docker_image

    Int cpus = 1
    Int disk_size = ceil(length(reports)*size(reports[0], "GB"))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        python3 /combine_ccs_reports.py ${sep=" " reports} > ccs_report.txt

        df -h .
        tree -h
    >>>

    output {
        File final_report = "ccs_report.txt"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "20G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        preemptible: 1
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
        preemptible: 1
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
    Int disk_size = ceil(size(bam_file, "GB") + size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(ref_dict, "GB"))

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
        preemptible: 1
    }
}

task Depth {
    File input_bam
    File input_bai
    String docker_image

    Int cpus = 2
    Int disk_size = ceil(size(input_bam, "GB") + size(input_bai, "GB"))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        samtools depth -r `samtools view -H ${input_bam} | head -2 | tail -1 | cut -f2 | sed 's/SN://'` -a ${input_bam} | awk '{ chr = $1; d = $3 - mean; mean += d/NR; M2 += d*($3 - mean); } END { print sprintf("%s %.2f %.2f", chr, mean, sqrt(M2/(NR - 1))); }' > coverage.txt

        df -h .
        tree -h
    >>>

    output {
        File wgsmetrics = "coverage.txt"
    }

    runtime {
        docker: "${docker_image}"
        cpu: "${cpus}"
        memory: "20G"
        bootDiskSizeGb: 20
        disks: "local-disk ${disk_size} SSD"
        preemptible: 1
    }
}