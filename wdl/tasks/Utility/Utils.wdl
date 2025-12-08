version 1.0

import "../../structs/Structs.wdl"

task ChunkManifest {

    meta {
        description: "Chunk a manifest file into smaller files"
    }

    parameter_meta {
        manifest: "The manifest file to chunk"
        manifest_lines_per_chunk: "The number of lines to include in each chunk"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File manifest
        Int manifest_lines_per_chunk

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        split -a 5 -d --additional-suffix=".txt" -l ~{manifest_lines_per_chunk} ~{manifest} chunk_
    >>>

    output {
        Array[File] manifest_chunks = glob("chunk_*")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        1,
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

task SortSam {
     meta {
        description: "Sort a BAM file by coordinate order"
    }

    parameter_meta {
        input_bam: "The BAM file to sort"
        prefix: "The basename for the output BAM file"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File input_bam
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 10*ceil(size(input_bam, "GB"))

    command <<<
        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        if [[ ${np} -gt 2 ]] ; then
            np=$((np-1))
        fi

        samtools sort -@ ${np} ~{input_bam} -o ~{prefix}.bam
        samtools index -@ ${np} ~{prefix}.bam
    >>>

    output {
        File output_bam = "~{prefix}.bam"
        File output_bam_index = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
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

task ChangeReadGroup {
     meta {
        description: "Change the read group of a BAM file."
    }

    parameter_meta {
        input_bam: "The BAM file to sort"
        ID: "The new ID for the read group"
        SM: "The new sample name for the read group"
        PL: "The new platform for the read group"
        LB: "The new library name for the read group"
        prefix: "The basename for the output BAM file"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File input_bam
        String prefix

        String ID
        String SM
        String PL
        String LB

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 10*ceil(size(input_bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools addreplacerg --no-PG \
            -r 'ID:~{ID}' \
            -r 'LB:~{LB}' \
            -r 'SM:~{SM}' \
            -r 'PL:~{PL}' \
            ~{input_bam} \
            -o ~{prefix}.bam
    >>>

    output {
        File output_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
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

task MakeChrIntervalList {

    meta {
        description: "Make a Picard-style list of intervals for each chromosome in the reference genome"
    }

    parameter_meta {
        ref_dict: "The reference dictionary"
        filter: "A list of strings to filter out of the reference dictionary"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File ref_dict
        Array[String] filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    command <<<
        set -euxo pipefail

        grep '^@SQ' ~{ref_dict} | \
            awk '{ print $2 "\t" 1 "\t" $3 }' | \
            sed 's/[SL]N://g' | \
            grep -v -e '^@HD' ~{true='-e' false='' length(filter) > 0} ~{sep=" -e " filter} | \
            tee chrs.txt

        cat chrs.txt | awk '{printf("%s:%d-%d\n", $1,$2,$3)}' > intervalList.intervals

        # Now make another output - a set of individual contig interval list files:
        while read line ; do
            contig=$(echo "${line}" | awk '{print $1}')
            echo "${line}" | awk '{printf("%s:%d-%d\n", $1,$2,$3)}' > contig.${contig}.intervals
        done < chrs.txt
    >>>

    output {
        Array[Array[String]] chrs = read_tsv("chrs.txt")
        File interval_list = "intervalList.intervals"
        Array[String] contig_interval_strings = read_lines("intervalList.intervals")
        Array[File] contig_interval_list_files = glob("contig.*.intervals")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
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

task ExtractIntervalNamesFromIntervalOrBamFile {

    meta {
        description: "Pulls the contig names and regions out of an interval list or bed file."
    }

    parameter_meta {
        interval_file: "Interval list or bed file from which to extract contig names and regions."
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File interval_file

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    String interval_tsv_filename = "intervals.tsv"

    command <<<
        set -euxo pipefail

        python3 <<CODE

        import re
        interval_re = re.compile('''(.*?):(\d+)-(\d+)''')

        def parse_bed_or_interval_file_for_names(interval_file):

            interval_names = []

            with open(interval_file, 'r') as f:
                if interval_file.endswith(".bed"):
                    for line in f:
                        fields = line.strip().split("\t")
                        contig = fields[0]
                        start = int(fields[1]) + 1  # Add 1 here so that the results will be valid *intervals* (bed files start from 0).
                        end = int(fields[2])

                        interval_names.append((contig, start, end))
                else:
                    for line in f:
                        match = interval_re.match(line.strip())
                        if not match:
                            raise RuntimeError(f"BAD NEWS!  DIDN'T MATCH OUR REGEX FOR INTERVAL LISTS!  ARE YOUR FILE EXTENSIONS CORRECT?  Line: {line}")
                        else:
                            contig = match.group(1)
                            start = int(match.group(2))
                            end = int(match.group(3))

                        interval_names.append((contig, start, end))
            return interval_names

        # Get our intervals:
        intervals_names = parse_bed_or_interval_file_for_names("~{interval_file}")

        # Print our interval names to stdout:
        print("Interval info:")
        for interval in intervals_names:
            print(f"{interval[0]}:{interval[1]}-{interval[2]}")
        print()

        # Write out the master interval list and individual interval files:
        with open("~{interval_tsv_filename}", 'w') as f:
            for interval in intervals_names:
                f.write(f"{interval[0]}\t{interval[1]}\t{interval[2]}\n")

        print(f"Wrote all interval info to: ~{interval_tsv_filename}")

        CODE
    >>>

    output {
        Array[Array[String]] interval_info = read_tsv(interval_tsv_filename)
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
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

task MakeIntervalListFromSequenceDictionary {

    meta {
        description: "Make a Picard-style list of intervals that covers the given reference genome dictionary, with intervals no larger than the given size limit."
    }

    parameter_meta {
        ref_dict: "The reference dictionary"
        ignore_contigs: "A list of strings to filter out of the reference dictionary"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File ref_dict
        Int max_interval_size = 10000
        Array[String] ignore_contigs = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    String out_master_interval_filename = "intervalList.intervals"
    String interval_tsv_filename = "intervals.tsv"

    command <<<
        set -euxo pipefail

        python3 <<CODE

        import re

        def create_intervals(sequence_dictionary_file: str, max_interval_size_bp: int, ignore_contigs: list = []) -> list:

            # First read in our sequence dictionary:
            sq_re = re.compile("@SQ[\t ]+SN:(\w+)[\t ]+LN:(\d+)[\t ].*")
            seq_dict = dict()
            with open(sequence_dictionary_file, 'r') as f:
                for line in f:
                    m = sq_re.match(line)
                    if m and m.group(1) not in ignore_contigs:
                        seq_dict[m.group(1)] = int(m.group(2))

            # Now create a list of intervals with size constraints:
            intervals = list()
            for contig, length in seq_dict.items():
                if length <= max_interval_size_bp:
                    intervals.append((contig, 1, length))
                else:
                    # We need to split the contig into parts:
                    for i in range(1, length+1, max_interval_size_bp):
                        if i+max_interval_size_bp-1 > length:
                            intervals.append((contig, i, length))
                        else:
                            intervals.append((contig, i, i+max_interval_size_bp-1))

                    if intervals[-1][-1] != length:
                        intervals.append((contig, i+max_interval_size_bp, length))

            return intervals

        # Get our intervals:
        intervals = create_intervals("~{ref_dict}", ~{max_interval_size}, ["~{sep='","' ignore_contigs}"])

        # Print our intervals to stdout:
        print(f"Generated {len(intervals)} intervals:")
        for interval in intervals:
            print(f"{interval[0]}:{interval[1]}-{interval[2]}")
        print()

        # Write out the master interval list and individual interval files:
        with open("~{out_master_interval_filename}", 'w') as f:
            with open("~{interval_tsv_filename}", 'w') as f2:
                for interval in intervals:
                    this_interval_filename = f"{interval[0]}.{interval[1]}_{interval[2]}.intervals"

                    f.write(f"{interval[0]}:{interval[1]}-{interval[2]}\n")
                    f2.write(f"{interval[0]}\t{interval[1]}\t{interval[2]}\n")

        print(f"Wrote all intervals to: ~{out_master_interval_filename}")
        print(f"Wrote individual intervals to separate named interval files.")

        CODE
    >>>

    output {
        File interval_list = out_master_interval_filename
        Array[Array[String]] interval_info = read_tsv(interval_tsv_filename)
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
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

task CreateIntervalListFileFromIntervalInfo {

    meta {
        description: "Make a Picard-style interval list file from the given interval info."
    }

    parameter_meta {
        contig: "Contig for the interval."
        start: "Start position for the interval."
        end: "End position for the interval."
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        String contig
        String start
        String end
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    String out_interval_list = "~{contig}.~{start}_~{end}.intervals"

    command <<<
        set -euxo pipefail

        echo "~{contig}:~{start}-~{end}" > ~{out_interval_list}
    >>>

    output {
        File interval_list = out_interval_list
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "ubuntu:22.04"
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

task CountBamRecords {

    meta {
        description: "Count the number of records in a bam file"
    }

    parameter_meta {
        bam: {
                 localization_optional: true,
                 description: "The bam file"
             }
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 100

    command <<<
        set -eux
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -c ~{bam} > "count.txt" 2>"error.log"
        if [[ -f "error.log" ]]; then
            if [[ -s "error.log" ]]; then echo "samtools has warn/error" && cat "error.log" && exit 1; fi
        fi
    >>>

    output {
        File? samools_error = "error.log"
        Int num_records = read_int("count.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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

task DownsampleSam {

    meta {
        description : "Downsample the given bam / sam file using Picard/GATK's DownsampleSam tool."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        bam:   "BAM file to be filtered."
        probability : "[Optional] Probability that a read will be emitted (Default: 0.01)."
        strategy : "[Optional] Strategy to use to downsample the given bam file (Default: HighAccuracy)."
        prefix : "[Optional] Prefix string to name the output file (Default: downsampled_reads)."
        extra_args : "[Optional] Extra arguments to pass into DownsampleSam."
    }

    input {
        File bam

        Float probability = 0.01
        String strategy = "HighAccuracy"
        String prefix = "downsampled_reads"

        Int random_seed = 1

        String extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(11 * size(bam, "GiB"))

    command <<<

        # Make sure we use all our proocesors:
        np=$(grep ^processor /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')

        gatk DownsampleSam --VALIDATION_STRINGENCY SILENT --RANDOM_SEED ~{random_seed} -I ~{bam} -O ~{prefix}.bam -S ~{strategy} -P ~{probability} ~{extra_args}
        samtools index -@$np ~{prefix}.bam
    >>>

    output {
        File output_bam = "~{prefix}.bam"
        File output_bam_index = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.2.0.0"
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

task Sum {

    meta {
        description : "Sum a list of integers."
    }

    parameter_meta {
        ints: "List of integers to be summed."
        prefix: "[Optional] Prefix string to name the output file (Default: sum)."
    }

    input {
        Array[Int] ints
        String prefix = "sum"

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        python -c "print(~{sep="+" ints})" > ~{prefix}.txt
    >>>

    output {
        Int sum = read_int("~{prefix}.txt")
        File sum_file = "~{prefix}.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
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

task Uniq {

    meta {
        description : "Find the unique elements in a list of strings."
    }

    parameter_meta {
        strings: "List of strings to be filtered."
    }

    input {
        Array[String] strings

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1

    command <<<
        set -euxo pipefail

        sort ~{write_lines(strings)} | uniq > uniq.txt
    >>>

    output {
        Array[String] unique_strings = read_lines("uniq.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        2,
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

task Timestamp {

    meta {
        description : "Get the current timestamp."
    }

    parameter_meta {
        dummy_dependencies: "List of dummy dependencies to force recomputation."
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        Array[String] dummy_dependencies

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        date --iso-8601=ns > timestamp.txt
    >>>

    output {
        String timestamp = read_string("timestamp.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       25,
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

task ConvertReads {

    meta {
        description : "Convert reads from one format to another."
    }

    parameter_meta {
        reads: "Reads to be converted."
        output_format: "Output format."
    }

    input {
        File reads
        String output_format
    }

    Int disk_size = 3 * ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        filename=~{reads}
        input_filetype=${filename##*.}
        output_filetype=~{output_format}

        if [[ ($input_filetype == "fastq" || $input_filetype == "fq") && $output_filetype == "fasta" ]]; then
            echo "Converting $input_filetype to $output_filetype"
            seqkit fq2fa $filename -o tmp.out
        elif [ $input_filetype == $output_filetype ]; then
            echo "Input filetype is the output filetype"
            mv $filename tmp.out
        else
            echo "ConvertReads does not know how to convert $input_filetype to $output_filetype"
            exit 1
        fi

        mv tmp.out converted_reads.$output_filetype
    >>>

    output {
        File converted_reads = "converted_reads.~{output_format}"
    }

    runtime {
        cpu:                    4
        memory:                 "8 GiB"
        disks:                  "local-disk " +  disk_size + " HDD"
        bootDiskSizeGb:         25
        preemptible:            2
        maxRetries:             0
        docker:                 "quay.io/broad-long-read-pipelines/lr-pacasus:0.3.0"
    }
}

task BamToBed {

    meta {
        description : "Convert a BAM file to a bed file."
    }

    parameter_meta {
        bam: "BAM file to be converted."
        prefix: "Prefix for the output bed file."
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail
        bedtools genomecov -ibam ~{bam} -bg > ~{prefix}.bed
    >>>

    output {
        File bed = "~{prefix}.bed"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
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

task BamToFastq {

    meta {
        description : "Convert a BAM file to a fastq file."
    }

    parameter_meta {
        bam: "BAM file to be converted."
        prefix: "Prefix for the output fastq file."
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools fastq ~{bam} | gzip > ~{prefix}.fq.gz
    >>>

    output {
        File reads_fq = "~{prefix}.fq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
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

task MergeFastqs {

    meta {
        description : "Merge fastq files."
    }

    parameter_meta {
        fastqs: "Fastq files to be merged."
        prefix: "Prefix for the output fastq file."
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        Array[File] fastqs

        String prefix = "merged"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3 * ceil(size(fastqs, "GB"))

    String disk_type = if disk_size < 375 then "LOCAL" else "HDD"

    Int memory = 8

    command <<<
        FILE="~{fastqs[0]}"
        if [[ "$FILE" =~ \.gz$ ]]; then
            cat ~{sep=' ' fastqs} > ~{prefix}.fq.gz
        else
            cat ~{sep=' ' fastqs} | gzip > ~{prefix}.fq.gz
        fi
    >>>

    output {
        File merged_fastq = "~{prefix}.fq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

# A utility to merge several input BAMs into a single BAM.
task MergeBams {

    meta {
        description : "Merge several input BAMs into a single BAM."
    }

    parameter_meta {
        bams: "Input array of BAMs to be merged."
        prefix: "Prefix for the output BAM."
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        Array[File] bams
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(bams, "GB"))

    command <<<
        set -euxo pipefail

        samtools merge \
            -p -c --no-PG \
            -@ 2 \
            --write-index \
            -o "~{prefix}.bam##idx##~{prefix}.bam.bai" \
            ~{sep=" " bams}
    >>>

    output {
        File merged_bam = "~{prefix}.bam"
        File merged_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Index {

    meta {
        description : "samtools index a BAM file."
    }

    parameter_meta {
        bam: "BAM file to be indexed."
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(bam, "GB"))

    String prefix = basename(bam)

    command <<<
        ################################
        # Standard Preamble

        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        if [[ ${np} -gt 2 ]] ; then
            let np=${np}-1
        fi

        tot_mem_mb=$(free -m | grep '^Mem' | awk '{print $2}')

        ################################

        mv ~{bam} ~{prefix}
        samtools index -@$((np-1)) ~{basename(prefix)}
    >>>

    output {
        File bai = "~{prefix}.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

# A utility to subset a BAM to specifed loci
task SubsetBam {

    meta {
        description : "Subset a BAM file to a specified locus."
    }

    parameter_meta {
        bam: {
            description: "bam to subset",
            localization_optional: true
        }
        bai:    "index for bam file"
        locus:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        File bai
        String locus
        String prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }



    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam
    >>>

    output {
        File subset_bam = "~{prefix}.bam"
        File subset_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

task ResilientSubsetBam {

    meta {
        description: "For subsetting a high-coverage BAM stored in GCS, without localizing (more resilient to auth. expiration)."
    }

    parameter_meta {
        bam: {
            localization_optional: true
        }
        interval_list_file:  "a Picard-style interval list file to subset reads with"
        interval_id:         "an ID string for representing the intervals in the interval list file"
        prefix: "prefix for output bam and bai file names"
    }

    input {
        File bam
        File bai

        File interval_list_file
        String interval_id
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Array[String] intervals = read_lines(interval_list_file)

    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    String subset_prefix = prefix + "." + interval_id

    command <<<

        # the way this works is the following:
        # 0) relying on the re-auth.sh script to export the credentials
        # 1) perform the remote sam-view subsetting in the background
        # 2) listen to the PID of the background process, while re-auth every 1200 seconds
        source /opt/re-auth.sh
        set -euxo pipefail

        # see man page for what '-M' means
        samtools view \
            -bhX \
            -M \
            -@ 1 \
            --verbosity=8 \
            --write-index \
            -o "~{subset_prefix}.bam##idx##~{subset_prefix}.bam.bai" \
            ~{bam} ~{bai} \
            ~{sep=" " intervals} && exit 0 || { echo "samtools seem to have failed"; exit 77; } &
        pid=$!

        set +e
        count=0
        while true; do
            sleep 1200 && date && source /opt/re-auth.sh
            count=$(( count+1 ))
            if [[ ${count} -gt 6 ]]; then exit 0; fi
            if ! pgrep -x -P $pid; then exit 0; fi
        done
    >>>

    output {
        File subset_bam = "~{subset_prefix}.bam"
        File subset_bai = "~{subset_prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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

task Bamtools {

    meta {
        description: "Runs a given bamtools command on a bam file"
    }

    parameter_meta {
        bamfile:    "bam file to run bamtools on"
        cmd:        "bamtools command to run"
        args:       "arguments to pass to bamtools"
    }

    input {
        File bamfile
        String cmd
        String args

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 1 + ceil(2 * size(bamfile, "GiB"))

    command <<<
        bamtools ~{cmd} -in ~{bamfile} -out ~{prefix}.bam ~{args}
    >>>

    output {
        File bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             2,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9.beta"
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


task DeduplicateBam {

    meta {
        description: "Utility to drop (occationally happening) duplicate records in input BAM"
    }

    parameter_meta {
        aligned_bam: "input BAM file"
        aligned_bai: "input BAM index file"
        same_name_as_input: "if true, output BAM will have the same name as input BAM, otherwise it will have the input basename with .dedup suffix"
        runtime_attr_override: "override default runtime attributes"
    }

    input {
        File aligned_bam
        File aligned_bai

        Boolean same_name_as_input = true

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3 * ceil(size(aligned_bam, "GB"))

    String base = basename(aligned_bam, ".bam")
    String prefix = if (same_name_as_input) then base else (base + ".dedup")

    command <<<
        echo "==========================================================="
        echo "collecting duplicate information"
        time \
            samtools view -@ 1 "~{aligned_bam}" | \
            awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5}' | \
            sort | uniq -d \
            > "~{aligned_bam}".duplicates.txt
        echo "==========================================================="
        echo "de-duplicating"
        time python3 /opt/remove_duplicate_ont_aln.py \
            --prefix "~{prefix}" \
            --annotations "~{aligned_bam}".duplicates.txt \
            "~{aligned_bam}"
        echo "==========================================================="
        echo "DONE"
        samtools index "~{prefix}.bam"
    >>>

    output {
        File corrected_bam = "~{prefix}.bam"
        File corrected_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.10"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Cat {

    meta {
        description: "Utility to concatenates a group of files into a single output file, with headers in the first line if has_header is true. If has_header is false, the script concatenates the files without headers."
    }

    parameter_meta {
        files:      "text files to combine"
        has_header: "files have a redundant header"
        out:        "[default-valued] output filename"
    }

    input {
        Array[File] files
        Boolean has_header = false
        String out = "out.txt"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(files, "GB"))

    command <<<
        set -euxo pipefail

        HAS_HEADER=~{true='1' false='0' has_header}

        if [ HAS_HEADER == 1 ]; then
            ((head -1 ~{files[0]}) && (cat ~{sep=' ' files} | xargs -n 1 tail -n +2)) > ~{out}
        else
            cat ~{sep=' ' files} > ~{out}
        fi
    >>>

    output {
        File combined = out
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

task ComputeGenomeLength {

    meta {
        description: "Utility to compute the length of a genome from a FASTA file"
    }

    parameter_meta {
        fasta:  "FASTA file"
    }

    input {
        File fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(fasta, "GB"))

    command <<<
        set -euxo pipefail

        samtools dict ~{fasta} | \
            grep '^@SQ' | \
            awk '{ print $3 }' | \
            sed 's/LN://' | \
            awk '{ sum += $1 } END { print sum }' > length.txt
    >>>

    output {
        Float length = read_float("length.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
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

task ListFilesOfType {

    meta {
        description: "Utility to list files of a given type in a directory"
    }

    parameter_meta {
        gcs_dir:  "input directory"
        suffixes: "suffix(es) for files"
        recurse:  "if true, recurse through subdirectories"
    }

    input {
        String gcs_dir
        Array[String] suffixes
        Boolean recurse = false

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1
    String in_dir = sub(gcs_dir, "/$", "")

    command <<<
        set -ex
        gsutil ls ~{true='-r' false='' recurse} ~{in_dir} > temp.txt
        grep -E '(~{sep="|" suffixes})$' temp.txt > files.txt || touch files.txt
        if [ ! -s files.txt ]; then echo "None found" && exit 1; fi
    >>>

    output {
        Array[String] files = read_lines("files.txt")
        File manifest = "files.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
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

task StopWorkflow {

    meta {
        description: "Utility to stop a workflow"
    }

    parameter_meta {
        reason: "reason for stopping"
    }

    input {
        String reason
    }
    command <<<
        echo -e "Workflow explicitly stopped because \n  ~{reason}." && exit 1
    >>>
    runtime {docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"}
}

task InferSampleName {
    meta {
        description: "Infer sample name encoded on the @RG line of the header section. Fails if multiple values found, or if SM ~= unnamedsample."
    }

    parameter_meta {
        bam: {
            localization_optional: true,
            description: "BAM file"
        }
    }

    input {
        File bam
        File bai
    }



    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} > header.txt
        if ! grep -q '^@RG' header.txt; then echo "No read group line found!" && exit 1; fi

        grep '^@RG' header.txt | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g' | sort | uniq > sample.names.txt
        if [[ $(wc -l sample.names.txt) -gt 1 ]]; then echo "Multiple sample names found!" && exit 1; fi
        if grep -iq "unnamedsample" sample.names.txt; then echo "Sample name found to be unnamedsample!" && exit 1; fi
    >>>

    output {
        String sample_name = read_string("sample.names.txt")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 25
        preemptible:    2
        maxRetries:     1
        docker:         "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task CheckOnSamplenames {

    meta {
        description: "Makes sure the provided sample names are all same, i.e. no mixture of sample names"
    }

    parameter_meta {
        sample_names: "sample names"
    }

    input {
        Array[String] sample_names
    }

    command <<<
        set -eux
        n_sm=$(sort ~{write_lines(sample_names)} | uniq | wc -l | awk '{print $1}')
        if [[ ${n_sm} -gt 1 ]]; then echo "Sample mixture!" && exit 1; fi
    >>>

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 25
        preemptible:    2
        maxRetries:     1
        docker:         "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

# todo: hook this in to all tasks using LOCAL ssds
task ComputeAllowedLocalSSD {
    # This exists because of the following error message
    #   Task PBFlowcell.ShardLongReads:NA:1 failed. The job was stopped before the command finished. PAPI error code 3.
    #   Execution failed: allocating: creating instance: inserting instance: Number of local SSDs for an instance of type custom-8-15360
    #   should be one of [0, 1, 2, 3, 4, 5, 6, 7, 8, 16, 24], while [9] is requested.
    meta {
        description: "Compute the number of LOCAL ssd's allowed by Google"
    }

    parameter_meta {
        intended_gb: "intended number of GB"
    }

    input {
        Int intended_gb
    }
        Int raw = intended_gb / 375
    command <<<
        if [[ ~{raw} -lt 1 ]]; then  ## we are pushing the boundary here a bit, based on the assumption that input is a convervative estimate
            echo "1" > "result.txt"
        elif [[ ~{raw} -lt 9 ]]; then
            echo ~{raw} > "result.txt"
        elif [[ ~{raw} -lt 16  ]]; then
            echo "16" > "result.txt"
        elif [[ ~{raw} -lt 24  ]]; then
            echo "24" > "result.txt"
        else
            echo "Would request ~{raw} local SSDs, more than possible (24)." && exit 1
        fi
    >>>

    output {
        Int numb_of_local_ssd = read_int("result.txt")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 25
        preemptible:    2
        maxRetries:     1
        docker:         "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task RandomZoneSpewer {

    meta {
        description: "Spews a random GCP zone"
        ## TODO: This is probably the right thing to do, but we need to test it:
        # volatile: true
    }

    parameter_meta {
        num_of_zones: "number of zones to spew"
    }

    input {
        Int num_of_zones

        String region = "us-central1"
    }

    command <<<
        set -eux

        # by no means a perfect solution, but that's not desired anyway

        # NOTE: as of May 19, 2025, requested zones must all be in the same region.
        #       The below code is a fix for that.

        rm -f zones.txt

        if [[ "~{region}" == "us-central1" ]]; then
            echo "us-central1-a" >> zones.txt
            echo "us-central1-b" >> zones.txt
            echo "us-central1-c" >> zones.txt
            echo "us-central1-f" >> zones.txt
        elif [[ "~{region}" == "us-east1" ]]; then
            echo "us-east1-b" >> zones.txt
            echo "us-east1-c" >> zones.txt
            echo "us-east1-d" >> zones.txt
            echo "us-east4-a" >> zones.txt
            echo "us-east4-b" >> zones.txt
            echo "us-east4-c" >> zones.txt
        elif [[ "~{region}" == "us-west1" ]]; then
            echo "us-west1-a" >> zones.txt
            echo "us-west1-b" >> zones.txt
            echo "us-west1-c" >> zones.txt
        elif [[ "~{region}" == "us-west2" ]]; then
            echo "us-west2-a" >> zones.txt
            echo "us-west2-b" >> zones.txt
            echo "us-west2-c" >> zones.txt
        elif [[ "~{region}" == "us-west3" ]]; then
            echo "us-west3-a" >> zones.txt
            echo "us-west3-b" >> zones.txt
            echo "us-west3-c" >> zones.txt
        elif [[ "~{region}" == "us-west4" ]]; then
            echo "us-west4-a" >> zones.txt
            echo "us-west4-b" >> zones.txt
            echo "us-west4-c" >> zones.txt
        else
            echo "Invalid region: ~{region}"
            echo "Defaulting to us-central1"

            echo "us-central1-a" >> zones.txt
            echo "us-central1-b" >> zones.txt
            echo "us-central1-c" >> zones.txt
            echo "us-central1-f" >> zones.txt
        fi
        
        shuf zones.txt | head -n ~{num_of_zones} > "chosen_zones.txt"

        cat chosen_zones.txt | tr '\n' ' ' | tr -d '\r' > "zone_string.txt"

        #########################################################

        echo "Zone list from which to randomly choose ~{num_of_zones}:"
        cat zones.txt

        echo 
        echo    

        echo "Final Zone List:"
        cat chosen_zones.txt

        echo 
        echo
        echo "Zone string:"
        cat zone_string.txt
    >>>

    output {
        Array[String] zones = read_lines("chosen_zones.txt")
        String zone_string = read_string("zone_string.txt")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 25
        preemptible:    2
        maxRetries:     1
        docker:         "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

# Get the current timestamp as a string.
# Levergaes the unix `date` command.
# You can enter your own date format string.
# The default date string is:
#     %Y%m%d_%H%M%S_%N
# which corresponds to a date of the following form:
# For August 10th, 2020 at 16:06:32.7091 EDT (20:06:32.7091 UTC):
#     20200810_200632_709100000
#
task GetCurrentTimestampString {

    meta {
        # The volatile keyword forces call caching to be turned off, which is
        # exactly what we want for this task.
        # For more info see: https://cromwell.readthedocs.io/en/stable/optimizations/VolatileTasks/
        volatile: true
        description: "Get the current timestamp as a string"
    }

    parameter_meta {
        date_format: "The date format string to use. See the unix `date` command for more info."
    }

    input {
        String date_format = "%Y%m%d_%H%M%S_%N"
    }

    String date_file = "the_date_file.txt"

    command {
        date +~{date_format} > ~{date_file}
        cat ~{date_file}
    }

    # ------------------------------------------------
    # Runtime settings:
     runtime {
         docker: "ubuntu:19.10"
         memory: "512 MB"
         disks: "local-disk 10 HDD"
         bootDiskSizeGb: 15
         preemptible: 0
         cpu: 1
     }

    output {
        String timestamp_string   = read_string(date_file)
    }
}


task GetRawReadGroup {

    meta {
        description: "Get the raw read group from a bam file (assumed to have 1 read group only)"
    }

    parameter_meta {
        gcs_bam_path: "path to bam file in GCS"
        runtime_attr_override: "override the runtime attributes"
    }

    input {
        String gcs_bam_path

        RuntimeAttr? runtime_attr_override
    }

    String out_file = "raw_read_group.txt"

    command {
        set -x

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

        # We need to escape the tabs and convert the spaces so that the read group will play nice with downstream processing:
        samtools view -H ~{gcs_bam_path} | grep -m1 '^@RG' | sed -e 's@\t@\\t@g' -e 's@ @_@g' > ~{out_file}

        echo "Raw Read Group:"
        cat ~{out_file}
    }

    output {
        String rg = read_string(out_file)
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            50,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.30"
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

task GetReadsInBedFileRegions {
    meta {
        desciption: "Get the reads from the given bam path which overlap the regions in the given bed file."
    }

    input {
        String gcs_bam_path
        File regions_bed

        String prefix = "reads"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        gcs_bam_path: "GCS URL to bam file from which to extract reads."
        regions_bed: "Bed file containing regions for which to extract reads."
        prefix:    "[default-valued] prefix for output BAM"
        runtime_attr_override: "Runtime attributes override struct."
    }

    Int disk_size = 2 * ceil(size([gcs_bam_path, regions_bed], "GB"))

    command <<<
        set -x
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

        # Make sure we use all our proocesors:
        np=$(grep ^processor /proc/cpuinfo | tail -n1 | awk '{print $NF+1}')

        samtools view -@${np} -b -h -L ~{regions_bed} ~{gcs_bam_path} | samtools sort - > ~{prefix}.bam
        samtools index -@${np} ~{prefix}.bam
    >>>

    output {
        File bam = "~{prefix}.bam"
        File bai = "~{prefix}.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.30"
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

task MapToTsv {

    meta {
        description: "Convert a map to a tsv file"
    }

    parameter_meta {
        my_map: "The map to convert"
        name_of_file: "The name of the file to write to"
    }

    input {
        Map[String, Float] my_map
        String name_of_file
    }

    command <<<
        cp ~{write_map(my_map)} ~{name_of_file}
    >>>

    output {
        File result = "~{name_of_file}"
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task CreateIGVSession {
    meta {
        description: "Create an IGV session given a list of IGV compatible file paths.  Adapted / borrowed from https://github.com/broadinstitute/palantir-workflows/blob/mg_benchmark_compare/BenchmarkVCFs ."
    }
    input {
        Array[String] input_bams
        Array[String] input_vcfs
        String reference_short_name
        String output_name

        RuntimeAttr? runtime_attr_override
    }

    Array[String] input_files = flatten([input_bams, input_vcfs])

    command {
        bash /usr/writeIGV.sh ~{reference_short_name} ~{sep=" " input_files} > "~{output_name}.xml"
    }

    output {
        File igv_session = "${output_name}.xml"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            50,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "quay.io/mduran/generate-igv-session_2:v1.0"
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

task SplitContigToIntervals {
    meta {
        author: "Jonn Smith"
        notes: "Splits the given contig into intervals of the given size."
    }

    input {
        File ref_dict
        String contig
        Int size = 200000

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2

    command <<<
        set -euxo pipefail

        awk '{print $2,$3}' ~{ref_dict} | grep '^SN' | sed -e 's@SN:@@' -e 's@LN:@@' | tr ' ' '\t' > genome.txt
        grep "~{contig}" genome.txt > genome.contig.txt

        bedtools makewindows -g genome.contig.txt -w ~{size} > ~{contig}.~{size}bp_intervals.bed

        max_pos=$(tail -n1 ~{contig}.~{size}bp_intervals.bed | awk '{print $3}')

        # Make individual bed files from each line:
        # NOTE: We need to add leading zeros here for sorting purposes.
        while read -r line ; do
            start=$(echo "${line}" | cut -d $'\t' -f 2)
            end=$(echo "${line}" | cut -d $'\t' -f 3)
            new_fn=$(printf "%s.%0${#max_pos}d-%0${#max_pos}d.single_interval.bed" ~{contig} "${start}" "${end}")
            echo "${line}" > "${new_fn}"
        done < ~{contig}.~{size}bp_intervals.bed
    >>>

    output {
        File full_bed_file = "~{contig}.~{size}bp_intervals.bed"
        Array[File] individual_bed_files = glob("*.single_interval.bed")
    }

    #########################
        RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
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

task ResolveMapKeysInPriorityOrder {
    meta {
        description: "Gets the first key in the map that exists.  If no keys exist, returns an empty string."
    }

    parameter_meta {
        map: "Map[String, String] to resolve." 
        keys: "Array[String] of keys to check in order of priority"
    }

    input {
        Map[String, String] map
        Array[String] keys
    }

    String out_file = "key.txt"

    command <<<
        touch ~{out_file}

        f=~{write_map(map)}
        awk '{print $1}' < ${f} > keys_in_map.txt
        while read key ; do 
            grep -q "^${key}$" keys_in_map.txt
            if [[ $? -eq 0 ]] ; then
                echo "${key}" > ~{out_file}
                exit 0
            fi
        done < ~{write_lines(keys)}
    >>>

    output {
        String key = read_string(out_file)
    }

    ###################
    runtime {
        cpu: 1
        memory:  "512 MiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 25
        preemptible:     3
        maxRetries:      2
        docker:"gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}