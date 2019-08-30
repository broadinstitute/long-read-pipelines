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

    scatter (gcs_dir in gcs_dirs) {
        call DetectRunInfo {
            input:
                gcs_dir = gcs_dir,
                sample_name = sample_name
        }

        call PrepareRun {
            input:
                files = DetectRunInfo.files
        }

        call ShardLongReads {
            input:
                unmapped_bam=PrepareRun.unmapped_bam
        }

        scatter (unmapped_shard in ShardLongReads.unmapped_shards) {
            call CCS {
                input:
                    unmapped_shard=unmapped_shard,
                    platform=DetectRunInfo.run_info['PL']
            }

            call Minimap2 as AlignCCS {
                input:
                    shard=CCS.ccs_shard,
                    ref_fasta=ref_fasta,
                    SM=DetectRunInfo.run_info['SM'],
                    ID=DetectRunInfo.run_info['ID'] + ".corrected",
                    PL=DetectRunInfo.run_info['PL'],
                    reads_are_corrected=true
            }

            call RecoverCCSRemainingReads {
                input:
                    unmapped_shard=unmapped_shard,
                    ccs_shard=CCS.ccs_shard
            }

            call Minimap2 as AlignRemaining {
                input:
                    shard=RecoverCCSRemainingReads.remaining_shard,
                    ref_fasta=ref_fasta,
                    SM=DetectRunInfo.run_info['SM'],
                    ID=DetectRunInfo.run_info['ID'] + ".remaining",
                    PL=DetectRunInfo.run_info['PL'],
                    reads_are_corrected=true
            }
        }

        call MergeBams as MergeCorrected {
            input:
                aligned_shards=AlignCCS.aligned_shard,
                merged_name="corrected.bam"
        }

        call MergeBams as MergeRemaining {
            input:
                aligned_shards=AlignRemaining.aligned_shard,
                merged_name="remaining.bam"
        }
    }

    call MergeBams as MergeAllCorrected {
        input:
            aligned_shards=MergeCorrected.merged,
            merged_name="all.corrected.bam"
    }

    call MergeBams as MergeAllRemaining {
        input:
            aligned_shards=MergeRemaining.merged,
            merged_name="all.remaining.bam"
    }
}

task DetectRunInfo {
    input {
        String gcs_dir
        String? sample_name
    }

    String SM = if defined(sample_name) then "--SM " + sample_name else ""

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        python /usr/local/bin/detect_run_info.py ~{SM} ~{gcs_dir} > run_info.txt
        gsutil ls ~{gcs_dir} | grep -v scraps | grep -e '\.bam$' -e '\.f\(ast\)\?q\(\.gz\)\?$' > files.txt
    >>>

    output {
        Map[String, String] run_info = read_map("run_info.txt")
        Array[String] files = read_lines("files.txt")
    }

    runtime {
        cpu: 1
        memory: "1GB"
        disks: "local-disk 1 SSD"
        preemptible: 1
        maxRetries: 0
        docker: "kgarimella/lr-test:0.01.13"
    }
}

task PrepareRun {
    input {
        Array[File] files
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
        docker: "kgarimella/lr-test:0.01.13"
    }
}

task ShardLongReads {
    input {
        File unmapped_bam
    }

    Int cpus = 1
    Int disk_size = 4*ceil(size(unmapped_bam, "GB"))

    command <<<
        set -euxo pipefail

        java -Dsamjdk.compression_level=0 -jar /usr/local/bin/gatk.jar ShardLongReads -I ~{unmapped_bam} -O ./ -DF WellformedReadFilter --use-jdk-deflater --use-jdk-inflater
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
        docker: "kgarimella/lr-test:0.01.13"
    }
}

task CCS {
    input {
        File unmapped_shard
        String platform
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
            touch ccs_report.txt
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
        docker: "kgarimella/lr-test:0.01.13"
    }
}

task RecoverCCSRemainingReads {
    input {
        File unmapped_shard
        File ccs_shard
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
        docker: "kgarimella/lr-test:0.01.13"
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
    }

    Boolean correct = select_first([reads_are_corrected, false])
    String map_arg = if (PL == "ONT") then "map-ont" else "map-pb"
    String correction_arg = if (correct) then "asm20" else map_arg

    String rg = "@RG\\tID:~{ID}\\tSM:~{SM}\\tPL:~{PL}"

    Int cpus = 4
    Int disk_size = ceil(size(ref_fasta, "GB")) + 4*ceil(size(shard, "GB"))

    String aligned_shard_name = basename(shard, ".bam") + ".aligned.bam"

    command <<<
        set -euxo pipefail

        samtools fastq ~{shard} | minimap2 -ayY --MD --eqx -x ~{correction_arg} -t ~{cpus} -R '~{rg}' ~{ref_fasta} - | samtools view -bS - > temp.aligned.unsorted.bam

        NUM_RECORDS=`samtools view temp.aligned.unsorted.bam | wc -l`
        if [ $NUM_RECORDS -eq 0 ]
        then
            cp temp.aligned.unsorted.bam temp.aligned.unsorted.repaired.bam
        else
            java -Dsamjdk.compression_level=0 -Xmx4g -jar /usr/local/bin/gatk.jar RepairLongReadBam -I ~{shard} -A temp.aligned.unsorted.bam -O temp.aligned.unsorted.repaired.bam -DF WellformedReadFilter --use-jdk-deflater --use-jdk-inflater
        fi

        samtools sort -@~{cpus} -m4G -o ~{aligned_shard_name} temp.aligned.unsorted.repaired.bam
    >>>

    output {
        File aligned_shard = "~{aligned_shard_name}"
    }

    runtime {
        cpu: "~{cpus}"
        memory: "20G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "kgarimella/lr-test:0.01.13"
    }
}

task MergeBams {
    input {
        Array[File] aligned_shards
        String merged_name
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
        docker: "kgarimella/lr-test:0.01.13"
    }
}
