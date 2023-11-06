version 1.0

import "../../../tasks/Utility/Finalize.wdl" as FF

workflow RepairBam {

    meta {
        description: "Repair a corrupted BAM file."
    }
    parameter_meta {
        input_bam: "The BAM file to repair"
        prefix: "The name to give the repaired BAM file"
        gcs_out_root_dir: "The root directory in GCS to store the output files"
    }

    input {
        File input_bam
        String prefix

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/RepairBam/~{prefix}"

    call ValidateSamFile { input: input_bam = input_bam }
    call Repair { input: input_bam = input_bam }

    # Finalize
    call FF.FinalizeToFile as FinalizeValidationLog { input: outdir = outdir, file = ValidateSamFile.validation_log }
    call FF.FinalizeToFile as FinalizeRepairedBam { input: outdir = outdir, file = Repair.repaired_bam }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File validation_log = FinalizeValidationLog.gcs_path
        File repaired_bam = FinalizeRepairedBam.gcs_path
    }
}

task ValidateSamFile {
    input {
        File input_bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(size(input_bam, "GiB"))
    String log_name = basename(input_bam, ".bam") + ".log"

    command <<<
        set -x

        gatk --java-options -Xms6000m \
            ValidateSamFile \
            I=~{input_bam} \
            O=~{log_name} \
            IGNORE_WARNINGS=true \
            MODE=VERBOSE
    >>>

    output {
        File validation_log = log_name
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-gatk/gatk:4.4.0.0"
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

task Repair {
    input {
        File input_bam
        Int N = 100000

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(input_bam, "GiB"))
    String fixed_bam = basename(input_bam, ".bam") + ".fixed.bam"

    command <<<
        set -x

        ((samtools view -H ~{input_bam}) && (samtools view ~{input_bam} | head -n ~{N})) | samtools view -bS > ~{fixed_bam}
    >>>

    output {
        File repaired_bam = fixed_bam
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-gatk/gatk:4.4.0.0"
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
