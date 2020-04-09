version 1.0

import "tasks/Structs.wdl"

workflow GuppyBasecaller {
    input {
        File fast5
    }

    call Basecall { input: fast5 = fast5 }
}

task Basecall {
    input {
        File fast5

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        apt-get install -y pciutils
        lspci

        nvidia-smi

        mkdir fast5
        mv ~{fast5} fast5/

        guppy_basecaller --compress_fastq -i fast5/ -s basecall/ -x "cuda:0" -c dna_r9.4.1_450bps_hac_prom.cfg

        ls -ahl basecall/

        cat basecall/sequencing_summary.txt

    >>>
    # --gpu_runners_per_device 1 --num_callers 1 --chunks_per_runner 1

    output {
        String out = read_string(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             16,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        gpuType:            "nvidia-tesla-p4",
        gpuCount:           1,
        nvidiaDriverVersion: "418.87.00",
        docker:             "us.gcr.io/broad-dsp-lrma/lr-guppy:0.1.0"
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
