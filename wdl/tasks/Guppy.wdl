version 1.0

import "Finalize.wdl" as FF

workflow Guppy {
    input {
        String gcs_fast5_dir
    }

    call ListFastqFiles {
        input:
            gcs_fast5_dir = gcs_fast5_dir
    }

    call Basecall {
        input:
           fast5_files = ListFastqFiles.fast5_files
    }

    output  {
        Array[File] output_files = Basecall.guppy_output_files
    }
}

task ListFastqFiles {
    input {
        String gcs_fast5_dir
    }

    String indir = sub(gcs_fast5_dir, "/$", "")

    command <<<
        gsutil ls ~{indir}/*.fast5> fast5_files.txt
    >>>

    output {
        Array[File] fast5_files = read_lines("fast5_files.txt")
    }

    runtime {
        cpu:                    1
        memory:                 "1 GiB"
        disks:                  "local-disk 1 HDD"
        bootDiskSizeGb:         10
        preemptible:            0
        maxRetries:             0
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
    }
}

task MergeFastq {
    input {
        Array[File] guppy_output_files
    }

    Int disk_size = 3 * ceil(size(guppy_output_files, "GB"))

    command <<<
        mkdir tmp
        mv -t tmp ~{sep=' ' guppy_output_files}

        cat tmp/*.fastq > merged.fastq
    >>>

    output {
        File merged_fastq = "merged.fastq"
    }

     runtime {
        cpu:                    1
        memory:                 "4 GiB"
        disks:                  "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:         10
        preemptible:            0
        maxRetries:             0
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.6"
    }
}

task Basecall {
    input {
        Array[File] fast5_files

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3 * ceil(size(fast5_files, "GB"))

    command <<<
        mkdir fast5
        mv -t fast5/ ~{sep=' ' fast5_files}

        guppy_basecaller -i fast5/ -s guppy_output/ -x "cuda:all" -c dna_r9.4.1_450bps_hac_prom.cfg
    >>>

    output {
        Array[File] guppy_output_files = glob("guppy_output/*")
    }

    runtime {
        cpu:                    4
        memory:                 "8 GiB"
        disks:                  "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:         30
        preemptible:            0
        maxRetries:             0
        gpuType:                "nvidia-tesla-p100"
        gpuCount:               1
        nvidiaDriverVersion:    "418.87.00"
        zones:                  ["us-east1-c"]
        cpuPlatform:            "Intel Haswell"
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-guppy:4.0.14"
    }

}
