version 1.0

import "tasks/Finalize.wdl" as FF

workflow GuppyBasecaller {
    input {
        String gcs_input_dir
        String gcs_output_dir
    }

    call ListFastqFiles {
        input:
            gcs_input_dir = gcs_input_dir
    }

    call Basecall {
        input:
           fast5_files = ListFastqFiles.fast5_files
    }

    call FF.FinalizeToDir {
        input:
            files = Basecall.guppy_output_files,
            outdir = gcs_output_dir
    }

}

task ListFastqFiles {
    input {
        String gcs_input_dir
    }

    String indir = sub(gcs_input_dir, "/$", "")

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

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       30,
        preemptible_tries:  0,
        max_retries:        0,
        gpuType:            "nvidia-tesla-p100",
        gpuCount:           1,
        nvidiaDriverVersion: "418.87.00",
        docker:             "quay.io/broad-long-read-pipelines/lr-guppy:0.4.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        gpuType:                "nvidia-tesla-p100"
        gpuCount:               1
        nvidiaDriverVersion:    "418.87.00"
        zones:                  ["us-central1-c"]
        cpuPlatform:            "Intel Haswell"
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

}
