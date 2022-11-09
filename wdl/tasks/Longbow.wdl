version 1.0

import "Structs.wdl"

task Annotate
{
    input {
        File reads
        String model = "mas_15+sc_10x5p"
        String prefix = "longbow_annotated"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow annotate -t8 --model ~{model} -v INFO ~{reads} -o ~{prefix}.bam
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.6.3"
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
        String model = "mas_15+sc_10x5p"
        String prefix = "longbow_segmented"

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 15*ceil(size(annotated_reads, "GB")) + 20

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow segment ~{extra_args} --model ~{model} -v INFO ~{annotated_reads} -o ~{prefix}.bam

        # Make sure this file exists:
        if [[ ! -e barcode_confidence_scores.txt ]] ; then
            touch barcode_confidence_scores.txt
        fi
    >>>

    output {
        File segmented_bam = "~{prefix}.bam"
        File barcode_conf_file = "barcode_confidence_scores.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.6.3"
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

task Sift
{
    input {
        File segmented_input_reads
        String model = "mas_15+sc_10x5p"

        String prefix = "longbow_sifted"

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 15*ceil(size(segmented_input_reads, "GB")) + 20

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate

        longbow sift ~{extra_args} \
            --model ~{model} \
            --stats ~{prefix}.stats.tsv \
            --summary-stats ~{prefix}.summary_stats.tsv \
            -o ~{prefix}.bam \
            -x ~{prefix}.sift_failed.bam \
            ~{segmented_input_reads}
    >>>

    output {
        File sifted_bam = "~{prefix}.bam"
        File sift_failed_bam = "~{prefix}.sift_failed.bam"
        File stats_tsv = "~{prefix}.stats.tsv"
        File summary_stats_tsv = "~{prefix}.summary_stats.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.6.3"
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

        String model = "mas_15+sc_10x5p"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.6.3"
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

        Array[String] models = ["mas10", "mas15"]

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        #set up memory logging daemon
        LOG_INTERVAL=5
        while true ; do
            echo "###################################"
            date
            cat /proc/meminfo
            sleep $LOG_INTERVAL
        done &
        pid=$!

        source /longbow/venv/bin/activate
        longbow demultiplex -v INFO --model ~{sep=' --model ' models} ~{bam} -o ~{prefix}

        # Create a list of models - one for each bam file created:
        # Do this safely (assume there can be spaces in the names even though this is generally bad form).
        # NOTE: the WDL glob() utility functions the same as bash globbing, so the order here should be the same:
        \ls ~{prefix}*.bam > tmp.txt
        while read file_name ; do
            echo "$file_name" | sed 's#^~{prefix}_\(.*\).bam$##g'
        done < tmp.txt >> file_model_list.txt

        kill -9 $pid
    >>>

    output {
        Array[File] demultiplexed_bams = glob("~{prefix}_*.bam")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.6.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
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
        String model = "mas_15+sc_10x5p"

        String prefix = "reads"

        File? bam_pbi

        RuntimeAttr? runtime_attr_override
    }

    String pbi_arg = if defined(bam_pbi) then " --pbi " else ""
    Int disk_size = ceil(4*ceil(size(bam, "GB")) + size(bam_pbi, "GB"))

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow filter \
            -v INFO \
            ~{pbi_arg}~{default="" sep=" --pbi " bam_pbi} \
            --model ~{model} \
            ~{bam} \
            -x ~{prefix}_longbow_filter_failed.bam \
            -o ~{prefix}_longbow_filter_passed.bam 2> >(tee longbow_filter_log.txt >&2) # Get log data from stderr and reprint to stderr
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.6.3"
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
        Int num_cpus = 2
        Int base_padding = 2

        Int? start_offset         # For mas15: 16+10
        String? leading_adapter   # For mas15: "10x_Adapter"
        String? trailing_adapter  # For mas15: "Poly_A"

        String prefix = "out"

        File? bam_pbi

        RuntimeAttr? runtime_attr_override
    }

    String pbi_arg = if defined(bam_pbi) then " --pbi " else ""
    String start_offset_arg = if defined(start_offset) then " --start-offset " else ""
    String leading_adapter_arg = if defined(leading_adapter) then " --leading-adapter " else ""
    String trailing_adapter_arg = if defined(trailing_adapter) then " --trailing-adapter " else ""

    Int disk_size = 10 * ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow extract \
            -v INFO \
            --base-padding ~{base_padding} \
            ~{start_offset_arg}~{default="" sep=" --start-offset " start_offset} \
            ~{leading_adapter_arg}~{default="" sep=" --leading-adapter " leading_adapter} \
            ~{trailing_adapter_arg}~{default="" sep=" --trailing-adapter " trailing_adapter} \
            ~{pbi_arg}~{default="" sep=" --pbi " bam_pbi} \
            ~{bam} \
            -o ~{prefix}.extracted.bam 2> >(tee longbow_extract_log.txt >&2) # Get log data from stderr and reprint to stderr
    >>>

    output {
        File extracted_bam = "~{prefix}.extracted.bam"
        File log = "longbow_extract_log.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             2*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.6.3"
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

task Pad
{
    input {
        File reads
        String model = "mas_15+sc_10x5p"
        String prefix = "longbow_padded"

        String tag_to_expand = "ZU"
        String new_tag_dest = "XM"
        Int padding = 2

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow pad --model ~{model} -v INFO --barcode-tag ~{tag_to_expand} -e ~{padding} -o tmp.bam -n ~{new_tag_dest} ~{reads}

        samtools sort tmp.bam -o ~{prefix}.bam
    >>>

    output {
        File padded_tag_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.6.3"
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

task Correct
{
    input {
        File reads
        File barcode_allow_list

        String model = "mas_15+sc_10x5p"
        String prefix = "longbow_correct"

        String raw_barcode_tag = "CR"
        String corrected_barcode_tag = "CB"

        Int ccs_lev_dist_threshold = 2
        Int clr_lev_dist_threshold = 3

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(reads, "GB")) + ceil(size(barcode_allow_list, "GB"))

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate
        longbow correct \
            -t 1 \
            --model ~{model} \
            --allow-list ~{barcode_allow_list} \
            -v INFO \
            --barcode-tag ~{raw_barcode_tag} \
            --corrected-tag ~{corrected_barcode_tag} \
            --max-hifi-dist ~{ccs_lev_dist_threshold} \
            --max-clr-dist ~{clr_lev_dist_threshold} \
            -o ~{prefix}_corrected_barcodes.bam \
            --barcode-uncorrectable-bam ~{prefix}_uncorrected_barcodes.bam \
            ~{reads} 2>&1 | tee longbow_correct.~{prefix}.log
    >>>

    output {
        File corrected_barcodes_bam = "~{prefix}_corrected_barcodes.bam"
        File uncorrected_barcodes_bam = "~{prefix}_uncorrected_barcodes.bam"
        File log = "longbow_correct.~{prefix}.log"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.6.3"
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

task AggregateCorrectLogStats
{
    input {
        Array[File] longbow_correct_log_files

        String out_name = "longbow_correct_stats.txt"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(longbow_correct_log_files, "GB"))

    # YES, this SHOULD be a proper tool, but right now it isn't.
    command <<<
python << CODE
import os

stats_dict = dict()
line_key = "STATS: "

for stats_file in ["~{sep='","' longbow_correct_log_files}"]:
    with open(stats_file, 'r') as f:
        for line in f:
            if line_key in line:
                line = line.strip()
                s = line[line.find(line_key) + len(line_key):]
                key, remainder = [t.strip() for t in s.split(":")]
                if "/" in remainder:
                    count = int(remainder[:remainder.find("/")])
                    tot = int(remainder[remainder.find("/")+1:remainder.find(" ")])
                else:
                    count = int(remainder)
                    tot = None

                try:
                    c, t = stats_dict[key]
                    if tot is not None:
                        tot += t
                    stats_dict[key] = (count + c, tot)
                except KeyError:
                    stats_dict[key] = (count, tot)

k_len = 0
for k in stats_dict.keys():
    if len(k) > k_len:
        k_len = len(k)

k_prefix = list(stats_dict.keys())[0]
k_prefix = k_prefix[:k_prefix.find(" ")]
with open("~{out_name}", 'w') as f:
    for k, v in stats_dict.items():

        if not k.startswith(k_prefix):
            f.write("\n")
            k_prefix = k[:k.find(" ")]

        k_spacing = k_len - len(k)

        count, tot = v
        if tot is None or tot == 0:
            f.write(f"{k}:{' '*k_spacing} {count}\n")
            print(f"WARNING: tot == {tot}")
        else:
            f.write(f"{k}:{' '*k_spacing} {count}/{tot} ({100.0*count/tot:2.4f}%)\n")

CODE
    >>>

    output {
        File stats = out_name
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.6.3"
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

task Stats
{
    input {
        File reads
        String model = "mas_15+sc_10x5p"
        String prefix = "longbow_stats"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        source /longbow/venv/bin/activate

        NUM_READS=`samtools view -c ~{reads}`
        if [ $NUM_READS -gt "0" ]; then
            longbow stats -s --model ~{model} -v INFO -o ~{prefix} ~{reads}
        else
            # This is a complete hack, meant to deal with simulated datasets that have zero reclaimed reads.
            # TODO: maybe fix longbow stats to just do the right thing with empty BAMs?

            touch ~{prefix}_summary_stats.txt
            touch ~{prefix}_00_MAS-seq_Array_Length_Counts_~{model}.png
            touch ~{prefix}_00_MAS-seq_Array_Length_Counts_~{model}.svg
            touch ~{prefix}_01_MAS-seq_Array_Element_Length_Counts_~{model}.png
            touch ~{prefix}_01_MAS-seq_Array_Element_Length_Counts_~{model}.svg
            touch ~{prefix}_02_MAS-seq_Ligations_~{model}_no_numbers.png
            touch ~{prefix}_02_MAS-seq_Ligations_~{model}_no_numbers.svg
            touch ~{prefix}_03_MAS-seq_Ligations_~{model}.png
            touch ~{prefix}_03_MAS-seq_Ligations_~{model}.svg
            touch ~{prefix}_04_MAS-seq_Ligations_~{model}_reduced_no_numbers.png
            touch ~{prefix}_04_MAS-seq_Ligations_~{model}_reduced_no_numbers.svg
            touch ~{prefix}_05_MAS-seq_Ligations_~{model}_reduced.png
            touch ~{prefix}_05_MAS-seq_Ligations_~{model}_reduced.svg
        fi
    >>>

    output {
        File summary_stats = "~{prefix}_summary_stats.txt"

        File array_length_counts_plot_png = "~{prefix}_00_MAS-seq_Array_Length_Counts_~{model}.png"
        File array_length_counts_plot_svg = "~{prefix}_00_MAS-seq_Array_Length_Counts_~{model}.svg"

        File array_element_length_counts_plot_png = "~{prefix}_01_MAS-seq_Array_Element_Length_Counts_~{model}.png"
        File array_element_length_counts_plot_svg = "~{prefix}_01_MAS-seq_Array_Element_Length_Counts_~{model}.svg"

        File ligation_heatmap_nn_png = "~{prefix}_02_MAS-seq_Ligations_~{model}_no_numbers.png"
        File ligation_heatmap_nn_svg = "~{prefix}_02_MAS-seq_Ligations_~{model}_no_numbers.svg"

        File ligation_heatmap_png = "~{prefix}_03_MAS-seq_Ligations_~{model}.png"
        File ligation_heatmap_svg = "~{prefix}_03_MAS-seq_Ligations_~{model}.svg"

        File ligation_heatmap_nn_reduced_png = "~{prefix}_04_MAS-seq_Ligations_~{model}_reduced_no_numbers.png"
        File ligation_heatmap_nn_reduced_svg = "~{prefix}_04_MAS-seq_Ligations_~{model}_reduced_no_numbers.svg"

        File ligation_heatmap_reduced_png = "~{prefix}_05_MAS-seq_Ligations_~{model}_reduced.png"
        File ligation_heatmap_reduced_svg = "~{prefix}_05_MAS-seq_Ligations_~{model}_reduced.svg"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.6.3"
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


task TagFix
{
    input {
        File bam

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        source /longbow/venv/bin/activate
        longbow tagfix -t${np} -v INFO -o tmp.bam ~{bam}

        samtools sort -@${np} tmp.bam -o ~{prefix}.alignment_tags_fixed.bam
    >>>

    output {
        File tag_fixed_bam = "~{prefix}.alignment_tags_fixed.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.6.3"
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

task Correct_UMI
{
    meta {
        description : "Corrects the UMIs in the given reads using a set cover algorithm."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }
    input {
        File bam

        String prefix = "out"
        Int umi_length = 10

        Int max_ccs_edit_dist = 2
        Int max_clr_edit_dist = 3
        Int max_ccs_length_diff = 50
        Int max_clr_length_diff = 150
        Float max_ccs_gc_diff = 0.05
        Float max_clr_gc_diff = 0.15
        Int max_ccs_umi_length_delta = 3
        Int max_clr_umi_length_delta = 4
        Int max_final_ccs_umi_length_delta = 3
        Int max_final_clr_umi_length_delta = 3
        Int min_back_seg_score = 10

        String umi_tag = "JX"
        String gene_tag = "XG"
        String eq_class_tag = "eq"
        String final_umi_tag = "BX"
        String umi_corrected_tag = "UX"

        Boolean pre_extracted = true

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 10*ceil(size(bam, "GB"))

    String pre_extracted_arg = if pre_extracted then " --pre-extracted " else ""

    command <<<
        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        source /longbow/venv/bin/activate
        longbow correct_umi \
            -v INFO \
            -l ~{umi_length} \
            ~{pre_extracted_arg} \
            --max-ccs-edit-dist ~{max_ccs_edit_dist} \
            --max-clr-edit-dist ~{max_clr_edit_dist} \
            --max-ccs-length-diff ~{max_ccs_length_diff} \
            --max-clr-length-diff ~{max_clr_length_diff} \
            --max-ccs-gc-diff ~{max_ccs_gc_diff} \
            --max-clr-gc-diff ~{max_clr_gc_diff} \
            --max-ccs-umi-length-delta ~{max_ccs_umi_length_delta} \
            --max-clr-umi-length-delta ~{max_clr_umi_length_delta} \
            --max-final-ccs-umi-length-delta ~{max_final_ccs_umi_length_delta} \
            --max-final-clr-umi-length-delta ~{max_final_clr_umi_length_delta} \
            --min-back-seg-score ~{min_back_seg_score} \
            --umi-tag ~{umi_tag} \
            --gene-tag ~{gene_tag} \
            --eq-class-tag ~{eq_class_tag} \
            --final-umi-tag ~{final_umi_tag} \
            --umi-corrected-tag ~{umi_corrected_tag} \
            -o ~{prefix}.corrected_umis.unsorted.bam \
            -x ~{prefix}.failed_umi_correction.bam \
            ~{bam}

        samtools sort -@${np} ~{prefix}.corrected_umis.unsorted.bam > ~{prefix}.corrected_umis.bam
        samtools index -@${np} ~{prefix}.corrected_umis.bam
    >>>

    output {
        File umi_corrected_bam = "~{prefix}.corrected_umis.bam"
        File umi_corrected_bam_index = "~{prefix}.corrected_umis.bam.bai"
        File failed_umi_correction_bam = "~{prefix}.failed_umi_correction.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.6.6"
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
