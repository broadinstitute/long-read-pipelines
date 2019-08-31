version 1.0

# Copyright Broad Institute, 2019
#
# About:
#   This WDL pipeline processes long read data from a single sample (which may be split among multiple PacBio
#   SMRTCells or Oxford Nanopore flowcells).  We perform a variety of tasks (including CCS error correction,
#   alignment, flowcell merging, de novo assembly, SNP-to-SV variant discovery, variant filtration, methylation
#   calling, and automated QC. The results of pipeline are intended to be "analysis-ready".
#
#   This pipeline is capable of processing any combination of the following datasets:
#     - PacBio raw CCS data
#     - PacBio raw CLR data
#     - Oxford Nanopore raw data
#
#   Data may be presented as a PacBio run directory, Oxford Nanopore run directory, or just a collection of
#   BAM/fastq files.  Generally the data type and sample information are autodetected, but can also be manually
#   overridden in the input JSON file.
#
#
# Description of inputs:
#   Required:
#       Array[String] gcs_dirs          - The GCS directories wherein the data is stored.
#
#   Optional:
#       String sample_name              - The sample name to use in place of the autodetected name.
#
#
# Licensing:
#   This script is released under the WDL source code license (BSD-3) (see LICENSE in
#   https://github.com/broadinstitute/wdl). Note however that the programs it calls may be subject to different
#   licenses. Users are responsible for checking that they are authorized to run all programs before running
#   this script.

workflow LRWholeGenomeSingleSample {
    input {
        Array[String] gcs_dirs
        File ref_fasta
        File ref_fasta_fai
        File ref_dict
        String? sample_name
    }

    String docker_align = "kgarimella/lr-align:0.01.15"
    String docker_asm = "kgarimella/lr-asm:0.01.00"

    scatter (gcs_dir in gcs_dirs) {
        call DetectRunInfo {
            input:
                gcs_dir = gcs_dir,
                sample_name = sample_name,
                docker = docker_align
        }

        call PrepareRun {
            input:
                files = DetectRunInfo.files,
                docker = docker_align
        }

        call ShardLongReads {
            input:
                unmapped_bam=PrepareRun.unmapped_bam,
                docker = docker_align
        }

        scatter (unmapped_shard in ShardLongReads.unmapped_shards) {
            call CCS {
                input:
                    unmapped_shard=unmapped_shard,
                    platform=DetectRunInfo.run_info['PL'],
                    docker = docker_align
            }

            call Minimap2 as AlignCCS {
                input:
                    shard=CCS.ccs_shard,
                    ref_fasta=ref_fasta,
                    SM=DetectRunInfo.run_info['SM'],
                    ID=DetectRunInfo.run_info['ID'] + ".corrected",
                    PL=DetectRunInfo.run_info['PL'],
                    reads_are_corrected=true,
                    docker = docker_align
            }

            call RecoverCCSRemainingReads {
                input:
                    unmapped_shard=unmapped_shard,
                    ccs_shard=CCS.ccs_shard,
                    docker = docker_align
            }

            call Minimap2 as AlignRemaining {
                input:
                    shard=RecoverCCSRemainingReads.remaining_shard,
                    ref_fasta=ref_fasta,
                    SM=DetectRunInfo.run_info['SM'],
                    ID=DetectRunInfo.run_info['ID'] + ".remaining",
                    PL=DetectRunInfo.run_info['PL'],
                    reads_are_corrected=false,
                    docker = docker_align
            }
        }

        call MergeBams as MergeCorrected {
            input:
                aligned_shards=AlignCCS.aligned_shard,
                merged_name="corrected.bam",
                docker = docker_align
        }

        call MergeBams as MergeRemaining {
            input:
                aligned_shards=AlignRemaining.aligned_shard,
                merged_name="remaining.bam",
                docker = docker_align
        }
    }

    call MergeBams as MergeAllCorrected {
        input:
            aligned_shards=MergeCorrected.merged,
            merged_name="all.corrected.bam",
            docker = docker_align
    }

    call ValidateBam as ValidateAllCorrected {
        input:
            input_bam=MergeAllCorrected.merged,
            docker = docker_align
    }

    call MergeBams as MergeAllRemaining {
        input:
            aligned_shards=MergeRemaining.merged,
            merged_name="all.remaining.bam",
            docker = docker_align
    }

    call ValidateBam as ValidateAllRemaining {
        input:
            input_bam=MergeAllRemaining.merged,
            docker = docker_align
    }
}

task DetectRunInfo {
    input {
        String gcs_dir
        String? sample_name
        String docker
    }

    String SM = if defined(sample_name) then "--SM " + sample_name else ""

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        python /usr/local/bin/detect_run_info.py ~{SM} ~{gcs_dir} > run_info.txt
        gsutil ls ~{gcs_dir} | grep -v scraps | grep -e '\.bam$' -e '\.f\(ast\)\?q\(\.gz\)\?$' > files.txt
    >>>

    output {
        File run_info_file = "run_info.txt"
        File fofn = "files.txt"
        Map[String, String] run_info = read_map("run_info.txt")
        Array[String] files = read_lines("files.txt")
    }

    runtime {
        cpu: 1
        memory: "1GB"
        disks: "local-disk 1 SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}

task PrepareRun {
    input {
        Array[File] files
        String docker
    }

    Int disk_size = 2*ceil(size(files, "GB"))

    command <<<
        set -euxo pipefail

        python /usr/local/bin/prepare_run.py ~{sep=' ' files}
    >>>

    output {
        File unmapped_bam = "unmapped.bam"
    }

    runtime {
        cpu: 2
        memory: "2GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}

task ShardLongReads {
    input {
        File unmapped_bam
        String docker
        Int? num_reads_per_split
    }

    Int cpus = 1
    Int disk_size = 4*ceil(size(unmapped_bam, "GB"))
    Int nr = select_first([num_reads_per_split, 10000])

    command <<<
        set -euxo pipefail

        java -Dsamjdk.compression_level=0 -jar /usr/local/bin/gatk.jar ShardLongReads -I ~{unmapped_bam} -nr ~{nr} -O ./ -DF WellformedReadFilter --use-jdk-deflater --use-jdk-inflater
    >>>

    output {
        Array[File] unmapped_shards = glob("*.bam")
    }

    runtime {
        cpu: "~{cpus}"
        memory: "2G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}

# Note: this task changes the incoming read group name
task CCS {
    input {
        File unmapped_shard
        String platform
        String docker
    }

    String ccs_shard_name = basename(unmapped_shard, ".bam") + ".ccs.corrected.bam"
    Int max_length = 21000
    Int min_passes = 2
    Int cpus = 4
    Int disk_size = 2*ceil(size(unmapped_shard, "GB"))

    command <<<
        set -euxo pipefail

        PLATFORM="~{platform}"

        if [ $PLATFORM == "ONT" ]
        then
            samtools view -H ~{unmapped_shard} | samtools view -b > ~{ccs_shard_name}

            echo "ZMWs input          (A)  : 0
ZMWs generating CCS (B)  : 0 (100.00%)
ZMWs filtered       (C)  : 0 (0.00%)

Exclusive ZMW counts for (C):
No usable subreads       : 0 (0.00%)
Below SNR threshold      : 0 (0.00%)
Lacking full passes      : 0 (0.00%)
Heteroduplexes           : 0 (0.00%)
Min coverage violation   : 0 (0.00%)
Draft generation error   : 0 (0.00%)
Draft above --max-length : 0 (0.00%)
Draft below --min-length : 0 (0.00%)
Lacking usable subreads  : 0 (0.00%)
CCS did not converge     : 0 (0.00%)
CCS below minimum RQ     : 0 (0.00%)
Unknown error            : 0 (0.00%)" > ccs_report.txt
        else
            ccs --max-length ~{max_length} --min-passes ~{min_passes} -j ~{cpus} ~{unmapped_shard} ~{ccs_shard_name}
        fi
    >>>

    output {
        File ccs_shard = "~{ccs_shard_name}"
        File ccs_report = "ccs_report.txt"
    }

    runtime {
        cpu: "~{cpus}"
        memory: "40G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}

task RecoverCCSRemainingReads {
    input {
        File unmapped_shard
        File ccs_shard
        String docker
    }

    String remaining_shard_name = basename(unmapped_shard, ".bam") + ".ccs.remaining.bam"
    Int cpus = 2
    Int disk_size = 2*ceil(size(unmapped_shard, "GB"))

    command <<<
        set -euxo pipefail

        java -Dsamjdk.compression_level=0 -Xmx4g -jar /usr/local/bin/gatk.jar RecoverUncorrectedReads -I ~{unmapped_shard} -C ~{ccs_shard} -O ~{remaining_shard_name} -DF WellformedReadFilter --use-jdk-deflater --use-jdk-inflater
    >>>

    output {
        File remaining_shard = "~{remaining_shard_name}"
    }

    runtime {
        cpu: "~{cpus}"
        memory: "20G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}

task Minimap2 {
    input {
        File shard
        File ref_fasta
        String SM
        String ID
        String PL
        Boolean? reads_are_corrected
        String docker
    }

    Boolean correct = select_first([reads_are_corrected, false])
    String map_arg = if (PL == "ONT") then "map-ont" else "map-pb"
    String correction_arg = if (correct) then "asm20" else map_arg

    Int cpus = 4
    Int disk_size = ceil(size(ref_fasta, "GB")) + 4*ceil(size(shard, "GB"))

    String aligned_shard_name = basename(shard, ".bam") + ".aligned.bam"

    command <<<
        set -euxo pipefail

        RG=`python /usr/local/bin/merge_read_group_tags.py --ID ~{ID} --SM ~{SM} --PL ~{PL} ~{shard}`
        samtools fastq ~{shard} | minimap2 -ayY --MD --eqx -x ~{correction_arg} -R ${RG} -t ~{cpus} ~{ref_fasta} - | samtools view -b - > temp.aligned.unsorted.bam
        java -Dsamjdk.compression_level=0 -Xmx4g -jar /usr/local/bin/gatk.jar RepairLongReadBam -I ~{shard} -A temp.aligned.unsorted.bam -O temp.aligned.unsorted.repaired.bam -DF WellformedReadFilter --use-jdk-deflater --use-jdk-inflater
        samtools sort -@~{cpus} -m4G -o ~{aligned_shard_name} temp.aligned.unsorted.repaired.bam
    >>>

    output {
        File temp_shard = "temp.aligned.unsorted.bam"
        File temp_repaired = "temp.aligned.unsorted.repaired.bam"
        File aligned_shard = "~{aligned_shard_name}"
    }

    runtime {
        cpu: "~{cpus}"
        memory: "20G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}

task MergeBams {
    input {
        Array[File] aligned_shards
        String merged_name
        String docker
    }

    Int cpus = 2
    Int disk_size = 3*ceil(size(aligned_shards, "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx4g -jar /usr/local/bin/gatk.jar MergeSamFiles -I ~{sep=" -I " aligned_shards} -O ~{merged_name} -AS --CREATE_INDEX --USE_JDK_DEFLATER --USE_JDK_INFLATER --VALIDATION_STRINGENCY SILENT
    >>>

    output {
        File merged = "~{merged_name}"
        File merged_bai = basename(merged_name, ".bam") + ".bai"
    }

    runtime {
        cpu: "~{cpus}"
        memory: "20G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}

task ValidateBam {
    input {
        File input_bam
        String docker
    }

    Int cpus = 2
    Int disk_size = ceil(size(input_bam, "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx4g -jar /usr/local/bin/gatk.jar ValidateSamFile -I ~{input_bam} -O bam_validation_report.txt
    >>>

    output {
        File report = "bam_validation_report.txt"
    }

    runtime {
        cpu: "~{cpus}"
        memory: "8G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}

task AssembleMT {
    input {
        File input_bam
        File input_bai
        String ref_fasta
        String mt_chr_name
        String SM
        String ID
        String docker
    }

    Int cpus = 4
    Int disk_size = ceil(2*(size(input_bam, "GB") + size(input_bai, "GB")))

    String reads = "reads.fastq.gz"

    command <<<
        set -euxo pipefail

        # select MT
        samtools view -hb ~{input_bam} ~{mt_chr_name} | samtools fastq - | gzip -1 > ~{reads}

        # align
        minimap2 -x ava-pb -t ~{cpus} ~{reads} ~{reads} | gzip -1 > reads.paf.gz

        # layout
        miniasm -f ~{reads} reads.paf.gz > reads.gfa
        awk '$1 ~/S/ { print "\x3E" $2 "\n" $3 }' reads.gfa > reads.fasta

        # correct 1
        minimap2 -t ~{cpus} reads.fasta ~{reads} > reads.gfa1.paf
        racon -t ~{cpus} ~{reads} reads.gfa1.paf reads.fasta > reads.racon1.fasta

        # correct 2
        minimap2 -t ~{cpus} reads.racon1.fasta ~{reads} > reads.gfa2.paf
        racon -t ~{cpus} ~{reads} reads.gfa2.paf reads.racon1.fasta > mt.fasta

        # align to ref
        minimap2 -ayY --MD --eqx -x asm20 -R '@RG\tID:MT\tSM:MT' -t ~{cpus} ~{ref_fasta} mt.fasta | samtools sort -@~{cpus} -m4G -o mt.bam
        samtools index mt.bam

        # call variants
        bcftools mpileup -Ou -f ~{ref_fasta} mt.bam | bcftools call -mv -Ov --ploidy 1 -o mt.vcf
    >>>

    output {
        File mt_asm = "mt.fasta"
        File mt_bam = "mt.bam"
        File mt_vcf = "mt.vcf"
    }

    runtime {
        cpu: "~{cpus}"
        memory: "4G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}
