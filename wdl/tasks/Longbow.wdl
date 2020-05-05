version 1.0

import "Structs.wdl"

task Annotate
{
    input {
        File reads
        String model = "mas15"
        String prefix = "longbow_annotated"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(reads, "GB"))

    String model_arg = if model == "mas15" then "" else "--m10"

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
#        longbow annotate -t8 --model ~{model} -v INFO ~{reads} -o ~{prefix}.bam
        longbow annotate -t8 ~{model_arg} -v INFO ~{reads} -o ~{prefix}.bam
    >>>

    output {
        File annotated_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
#        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.3.0"
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.2.2"
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

task Segment
{
    input {
        File annotated_reads
        String model = "mas15"
        String prefix = "longbow_segmented"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 15*ceil(size(annotated_reads, "GB")) + 20

    String model_arg = if model == "mas15" then "" else "--m10"

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
#        longbow segment --model ~{model} -v INFO -s ~{annotated_reads} -o ~{prefix}.bam
        longbow segment ~{model_arg}  -v INFO -s ~{annotated_reads} -o ~{prefix}.bam
    >>>

    output {
        File segmented_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.2.2"
#        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.3.0"
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

task ScSplit
{
    input {
        File reads_bam
        String model = "mas15"

        String force_option = ""

        Int umi_length = 10
        String cbc_dummy = "CTGCCTAACCTGATCC"

        String prefix = "scsplit"

        RuntimeAttr? runtime_attr_override
    }

    # On average bam compression is 10x and we have more data on output than that:
    Int disk_size = 15*ceil(size(reads_bam, "GB"))

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow scsplit -t4 -v INFO --model ~{model} -b ~{force_option} -o ~{prefix} -u ~{umi_length} -c ~{cbc_dummy} ~{reads_bam}
    >>>

    output {
        File mates1 = "~{prefix}_mates1.fastq"
        File mates2 = "~{prefix}_mates2.fastq"
        File annotated_bam = "~{prefix}.cbc_umi_annotated.bam"
        File whitelist = "~{prefix}_whitelist.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.3.0"
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

task Inspect
{
    input {
        File reads
        File reads_pbi

        File read_names

        String model = "mas15"
        String prefix = "longbow_inspected_reads"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(reads, "GB")) + size(reads_pbi, "GB") + size(read_names, "GB")

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow inspect --model ~{model} ~{reads} -p ~{reads_pbi} -r ~{read_names} -o ~{prefix}
        tar -zxf ~{prefix}.tar.gz ~{prefix}
    >>>

    output {
        File inspected_reads_tgz = "~{prefix}.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.3.0"
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

task Demultiplex
{
    input {
        File bam
        String prefix = "longbow_demultiplex"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow demultiplex -v INFO ~{bam} -o ~{prefix}

        # Create a list of models - one for each bam file created:
        # Do this safely (assume there can be spaces in the names even though this is generally bad form).
        # NOTE: the WDL glob() utility functions the same as bash globbing, so the order here should be the same:
        \ls ~{prefix}*.bam > tmp.txt
        while read file_name ; do
            echo "$file_name" | sed 's#^~{prefix}_\(.*\).bam$##g'
        done < tmp.txt >> file_model_list.txt
    >>>

    output {
        # TODO: Fix this to allow for an arbitrary number of models easily:
        File mas10_bam = "~{prefix}_mas10.bam"
        File mas15_bam = "~{prefix}_mas15.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
#        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.3.0"
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.2.2"
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

task CheckForAnnotatedArrayReads {
    input {
        String bam

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        gsutil cat ~{bam} | samtools view - | head -n1 | grep -q '[ \t]*SG:Z:'
        r=$?
        if [ $r -eq 1 ]; then
            echo "false" > bam_has_annotations.txt
        else
            echo "true" > bam_has_annotations.txt
        fi
    >>>

    output {
        # TODO: Fix this to allow for an arbitrary number of models easily:
        Boolean bam_has_annotations = read_boolean("bam_has_annotations.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             2,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
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

task Filter {
    input {
        File bam
        String model = "mas15"

        String prefix = "reads"


        File? bam_pbi

        RuntimeAttr? runtime_attr_override
    }

    String pbi_arg = if defined(bam_pbi) then " --pbi " else ""
    Int disk_size = ceil(4*ceil(size(bam, "GB")) + size(bam_pbi, "GB"))

    String model_arg = if model == "mas15" then "" else "--m10"

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
#        longbow filter \
#            -v INFO \
#            ~{pbi_arg}~{default="" sep=" --pbi " bam_pbi} \
#            --model ~{model} \
#            ~{bam} \
#            -o ~{prefix} 2> >(tee longbow_filter_log.txt >&2) # Get log data from stderr and reprint to stderr

        longbow filter \
            -v INFO \
            ~{pbi_arg}~{default="" sep=" --pbi " bam_pbi} \
            ~{model_arg} \
            ~{bam} \
            -o ~{prefix} 2> >(tee longbow_filter_log.txt >&2) # Get log data from stderr and reprint to stderr
    >>>

    output {
        File passed_reads = "~{prefix}_longbow_filter_passed.bam"
        File failed_reads = "~{prefix}_longbow_filter_failed.bam"
        File log = "longbow_filter_log.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.2.2"
#        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.3.0"
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


task Extract {
    input {
        File bam

        Int start_offset = 16+10
        Int base_padding = 2
        String leading_adapter = "10x_Adapter"
        String trailing_adapter = "Poly_A"

        String prefix = "reads"

        File? bam_pbi

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(4*ceil(size(bam, "GB")) + size(bam_pbi, "GB"))
    String pbi_arg = if defined(bam_pbi) then " --pbi " else ""

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow extract \
            -v INFO \
            --start-offset ~{start_offset} \
            --base-padding ~{base_padding} \
            --leading-adapter ~{leading_adapter} \
            --trailing-adapter ~{trailing_adapter} \
            ~{pbi_arg}~{default="" sep=" --pbi " bam_pbi} \
            ~{bam} \
            -o ~{prefix}.bam 2> >(tee longbow_extract_log.txt >&2) # Get log data from stderr and reprint to stderr
    >>>

    output {
        File extracted_reads = "~{prefix}.bam"
        File log = "longbow_extract_log.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.2.2"
#        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.3.0"
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
