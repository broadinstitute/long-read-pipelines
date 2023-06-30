version 1.0

import "../../../tasks/Utility/Finalize.wdl" as FF

workflow PrepareReferenceMap {
    input {
        String ref_name
        File ref_fasta
        File ref_fai

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PrepareReferenceMap/"

    call GenerateSequenceDict { input: fasta = ref_fasta }

    call FF.FinalizeToFile as FinalizeRefDict { input: outdir = outdir, file = GenerateSequenceDict.dict }

    call CreateRefMap {
        input:
            ref_name = ref_name,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = FinalizeRefDict.gcs_path
    }

    call FF.FinalizeToFile as FinalizeRefMap  { input: outdir = outdir, file = CreateRefMap.ref_map }

    output {
        File ref_dict = FinalizeRefDict.gcs_path
        File ref_map = FinalizeRefMap.gcs_path
    }
}

task GenerateSequenceDict {
    input {
        File fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(fasta, "GB"))
    String basename = basename(fasta, '.fasta')

    command <<<
        set -euxo pipefail

        samtools dict "~{fasta}" > "~{basename}.dict"
    >>>

    output {
        File dict = "~{basename}.dict"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CreateRefMap {
    input {
        String ref_name
        String ref_fasta
        File ref_fai
        File ref_dict

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(ref_fasta, "GB"))
    String basename = basename(ref_fasta, '.fasta')

    command <<<
        set -euxo pipefail

        echo -e 'fasta\t{ref_fasta}' >> ~{ref_name}.txt
        echo -e 'fai\t{ref_fai}' >> ~{ref_name}.txt
        echo -e 'dict\t{ref_dict}' >> ~{ref_name}.txt
    >>>

    output {
        File ref_map = "~{ref_name}.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

