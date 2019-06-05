workflow PBSingleSampleWorkflow {
    String input_bam
    String sm
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    File trf
    Int num_reads_per_split=100000

    String output_prefix="."
    String base_image="kgarimella/pbtools@sha256:304d90d8130fc8a8704eef566a7de7b614e386e1c02964ec7821f9bd0c25a046"

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

    call SplitIntervalsByChr {
        input:
            ref_fasta=ref_fasta,
            ref_fasta_fai=ref_fasta_fai,
            ref_dict=ref_dict,
            docker_image=base_image
    }

    scatter (intervals_file in SplitIntervalsByChr.interval_files) {
        call Depth as DepthCorrected {
            input:
                input_bam=MergeCorrected.merged,
                input_bai=MergeCorrected.merged_bai,
                ref_fasta=ref_fasta,
                ref_fasta_fai=ref_fasta_fai,
                intervals_file=intervals_file,
                docker_image=base_image
            }
<<<<<<< HEAD

        call Depth as DepthRemaining {
            input:
                input_bam=MergeRemaining.merged,
                input_bai=MergeRemaining.merged_bai,
                ref_fasta=ref_fasta,
                ref_fasta_fai=ref_fasta_fai,
                intervals_file=intervals_file,
                docker_image=base_image
            }

        call Depth as DepthUncorrected {
            input:
                input_bam=MergeUncorrected.merged,
                input_bai=MergeUncorrected.merged_bai,
                ref_fasta=ref_fasta,
                ref_fasta_fai=ref_fasta_fai,
                intervals_file=intervals_file,
                docker_image=base_image
            }

        call HaplotypeCaller {
            input:
                input_bam=MergeCorrected.merged,
                input_bai=MergeCorrected.merged_bai,
                ref_fasta=ref_fasta,
                ref_fasta_fai=ref_fasta_fai,
                ref_dict=ref_dict,
                intervals_file=intervals_file,
                calls=basename(MergeCorrected.merged, ".bam") + ".hc.vcf",
                docker_image=base_image
        }
    }

    call MergeVCFs as MergeHaplotypeCallerVCFs {
        input:
            input_vcfs = HaplotypeCaller.calls_vcf,
            input_vcf_indices = HaplotypeCaller.calls_idx,
            output_vcf = basename(MergeCorrected.merged, ".bam") + ".hc.vcf",
            docker_image=base_image
    }
=======
    }

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
>>>>>>> 0c345eea93ffdef78057021d4ac9ce17a2cdb8ef

    call PBSV {
        input:
            input_bam=MergeCorrected.merged,
            input_bai=MergeCorrected.merged_bai,
            ref_fasta=ref_fasta,
            ref_fasta_fai=ref_fasta_fai,
            calls=basename(MergeCorrected.merged, ".bam") + ".pbsv.vcf",
            docker_image=base_image
    }

    call Sniffles {
        input:
            input_bam=MergeCorrected.merged,
            input_bai=MergeCorrected.merged_bai,
            calls=basename(MergeCorrected.merged, ".bam") + ".sniffles.vcf",
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

#    call Depth as DepthUncorrected {
#        input:
#            input_bam=MergeUncorrected.merged,
#            input_bai=MergeUncorrected.merged_bai,
#            ref_fasta=ref_fasta,
#            ref_fasta_fai=ref_fasta_fai,
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
#    call Depth as DepthRemaining {
#        input:
#            input_bam=MergeRemaining.merged,
#            input_bai=MergeRemaining.merged_bai,
#            ref_fasta=ref_fasta,
#            ref_fasta_fai=ref_fasta_fai,
#            docker_image=base_image
#    }
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
        preemptible: 3
        maxRetries: 3
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

    Int cpus = 3
    Int disk_size = ceil(size(ref_fasta, "GB")) + 4*ceil(size(subread_file, "GB"))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

<<<<<<< HEAD
        samtools fastq ${subread_file} | minimap2 -ayY --MD --eqx -x ${correction_arg} -t ${cpus} ${ref_fasta} - | samtools view -bS - > temp.aligned.unsorted.bam
=======
        samtools fastq ${subread_file} | minimap2 -ayY --MD -x ${correction_arg} -t ${cpus} ${ref_fasta} - | samtools view -bS - > temp.aligned.unsorted.bam
>>>>>>> 0c345eea93ffdef78057021d4ac9ce17a2cdb8ef
        java -Dsamjdk.compression_level=0 -Xmx4g -jar /gatk.jar RepairPacBioBam -I ${subread_file} -A temp.aligned.unsorted.bam -O temp.aligned.unsorted.repaired.bam -S ${sample_name} --use-jdk-deflater --use-jdk-inflater
        samtools sort -@${cpus} -m4G -o ${subread_aligned} temp.aligned.unsorted.repaired.bam

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
        preemptible: 3
        maxRetries: 3
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
        #preemptible: 3
        maxRetries: 3
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
        preemptible: 3
        maxRetries: 3
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
        preemptible: 3
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
    File intervals_file
    String docker_image

    Int cpus = 2
<<<<<<< HEAD
    Int disk_size = ceil((size(input_bam, "GB") + size(input_bai, "GB") + size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(intervals_file, "GB")) * 2.0)
=======
    Int disk_size = ceil((size(input_bam, "GB") + size(input_bai, "GB") + size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(intervals_file, "GB")) * 1.1)
>>>>>>> 0c345eea93ffdef78057021d4ac9ce17a2cdb8ef

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        tail -1 ${intervals_file} | cut -f1 > chr.txt
<<<<<<< HEAD
        samtools depth -aa -r `cat chr.txt` ${input_bam} | bgzip > depth.txt.gz
=======
        samtools depth -aa -r `cat chr.txt` ${input_bam} | bgzip > `cut -f1 chr.txt`.depth.txt.gz
>>>>>>> 0c345eea93ffdef78057021d4ac9ce17a2cdb8ef
        tabix -s 1 -b 2 -e 2 depth.txt.gz

        df -h .
        tree -h
    >>>

    output {
<<<<<<< HEAD
        File chr = "chr.txt"
        File depth = "depth.txt.gz"
        File depth_index = "depth.txt.gz.tbi"
=======
        Array[File] depth = glob("*.depth.txt.gz")
        Array[File] depth_index = glob("*.depth.txt.gz.tbi")
>>>>>>> 0c345eea93ffdef78057021d4ac9ce17a2cdb8ef
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
<<<<<<< HEAD
    File intervals_file
=======
>>>>>>> 0c345eea93ffdef78057021d4ac9ce17a2cdb8ef
    String calls
    String docker_image

    Int cpus = 2
<<<<<<< HEAD
    Int disk_size = ceil((size(input_bam, "GB") + size(input_bai, "GB") + size(ref_fasta, "GB") + size(ref_dict, "GB") + size(ref_fasta_fai, "GB") + size(intervals_file, "GB")) * 1.1)
=======
    Int disk_size = ceil((size(input_bam, "GB") + size(input_bai, "GB") + size(ref_fasta, "GB") + size(ref_dict, "GB") + size(ref_fasta_fai, "GB")) * 1.1)
>>>>>>> 0c345eea93ffdef78057021d4ac9ce17a2cdb8ef

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

<<<<<<< HEAD
        java -Xmx16g -jar /gatk.jar HaplotypeCaller -R ${ref_fasta} -I ${input_bam} -L ${intervals_file} -O ${calls} --use-jdk-deflater --use-jdk-inflater -VS SILENT
=======
        java -Xmx16g -jar /gatk.jar HaplotypeCaller -R ${ref_fasta} -I ${input_bam} -O ${calls} --use-jdk-deflater --use-jdk-inflater -VS SILENT
>>>>>>> 0c345eea93ffdef78057021d4ac9ce17a2cdb8ef

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

<<<<<<< HEAD
task MergeVCFs {
    Array[File] input_vcfs
    Array[File] input_vcf_indices
    String output_vcf
    String docker_image

    Int cpus = 2
    Int disk_size = ceil(3*(size(input_vcfs[0], "GB")*length(input_vcfs)))

    command <<<
        set -euxo pipefail
        df -h .
        tree -h

        java -Xmx16g -jar /gatk.jar MergeVcfs -I ${sep=' -I ' input_vcfs} -O ${output_vcf}

        df -h .
        tree -h
    >>>

    output {
        File calls_vcf = "${output_vcf}"
        File calls_idx = "${output_vcf}.idx"
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

=======
>>>>>>> 0c345eea93ffdef78057021d4ac9ce17a2cdb8ef
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
