version 1.0

import "Structs.wdl"

# A wrapper to minimap2 for mapping & aligning (groups of) sequences to a reference
task Minimap2 {
    input {
        Array[File] reads
        File ref_fasta

        String RG
        String map_preset

        String? library

        Array[String] tags_to_preserve = []

        String prefix = "out"
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:            "query sequences to be mapped and aligned"
        ref_fasta:        "reference fasta"
        RG:               "read group information to be supplied to parameter '-R' (note that tabs should be input as '\t')"
        map_preset:       "preset to be used for minimap2 parameter '-x'"
        tags_to_preserve: "sam tags to carry over to aligned bam file"
        prefix:           "[default-valued] prefix for output BAM"
    }

    Boolean fix_library_entry = if defined(library) then true else false

    Int disk_size = 1 + 10*2*2*ceil(size(reads, "GB") + size(ref_fasta, "GB"))

    Boolean do_preserve_tags = if length(tags_to_preserve) != 0 then true else false

    Int cpus = 4
    Int mem = 30

    command <<<
        set -euxo pipefail

        NUM_CPUS=$( cat /proc/cpuinfo | grep '^processor' | tail -n1 | awk '{print $NF+1}' )
        RAM_IN_GB=$( free -g | grep "^Mem" | awk '{print $2}' )
        MEM_FOR_SORT=$( echo "" | awk "{print int(($RAM_IN_GB - 1)/$NUM_CPUS)}" )

        rg_len=$(echo -n '~{RG}' | wc -c | awk '{print $NF}')
        if [[ $rg_len -ne 0 ]] ; then
            MAP_PARAMS="-ayYL --MD --eqx -x ~{map_preset} -R ~{RG} -t ${NUM_CPUS} ~{ref_fasta}"
        else
            MAP_PARAMS="-ayYL --MD --eqx -x ~{map_preset} -t ${NUM_CPUS} ~{ref_fasta}"
        fi

        SORT_PARAMS="-@${NUM_CPUS} -m${MEM_FOR_SORT}G --no-PG -o ~{prefix}.pre.bam"
        FILE="~{reads[0]}"
        FILES="~{sep=' ' reads}"

        # We write to a SAM file before sorting and indexing because rarely, doing everything
        # in a single one-liner leads to a truncated file error and segmentation fault of unknown
        # origin.  Separating these commands requires more resources, but is more reliable overall.

        if [[ "$FILE" =~ \.fastq$ ]] || [[ "$FILE" =~ \.fq$ ]]; then
            cat $FILES | minimap2 $MAP_PARAMS - > tmp.sam
        elif [[ "$FILE" =~ \.fastq.gz$ ]] || [[ "$FILE" =~ \.fq.gz$ ]]; then
            zcat $FILES | minimap2 $MAP_PARAMS - > tmp.sam
        elif [[ "$FILE" =~ \.fasta$ ]] || [[ "$FILE" =~ \.fa$ ]]; then
            cat $FILES | python3 /usr/local/bin/cat_as_fastq.py | minimap2 $MAP_PARAMS - > tmp.sam
        elif [[ "$FILE" =~ \.fasta.gz$ ]] || [[ "$FILE" =~ \.fa.gz$ ]]; then
            zcat $FILES | python3 /usr/local/bin/cat_as_fastq.py | minimap2 $MAP_PARAMS - > tmp.sam
        elif [[ "$FILE" =~ \.bam$ ]]; then

            # samtools fastq takes only 1 file at a time so we need to merge them together:
            for f in "~{sep=' ' reads}" ; do
                if ~{do_preserve_tags} ; then
                    samtools fastq -T  ~{sep=',' tags_to_preserve} "$f"
                else
                    samtools fastq "$f"
                fi
            done > tmp.fastq

            echo "Memory info:"
            cat /proc/meminfo
            echo ""

            if ~{do_preserve_tags} ; then
                minimap2 ${MAP_PARAMS} -y tmp.fastq > tmp.sam
            else
                minimap2 ${MAP_PARAMS} tmp.fastq > tmp.sam
            fi
        else
            echo "Did not understand file format for '$FILE'"
            exit 1
        fi

        samtools sort ${SORT_PARAMS} tmp.sam

        if ~{fix_library_entry}; then
            mv ~{prefix}.pre.bam ~{prefix}.pre.tmp.bam
            samtools view --no-PG -H ~{prefix}.pre.tmp.bam > header.txt
            awk '$1 ~ /^@RG/' header.txt > rg_line.txt
            awk -v lib="~{library}" 'BEGIN {OFS="\t"} { for (i=1; i<=NF; ++i) { if ($i ~ "LB:") $i="LB:"lib } print}' \
                rg_line.txt \
                > fixed_rg_line.txt
            sed -n '/@RG/q;p' header.txt > first_half.txt
            sed -n '/@RG/,$p' header.txt | sed '1d' > second_half.txt

            cat first_half.txt fixed_rg_line.txt second_half.txt > fixed_header.txt

            date
            samtools reheader fixed_header.txt ~{prefix}.pre.tmp.bam > ~{prefix}.pre.bam
            rm ~{prefix}.pre.tmp.bam
            date
        fi

        samtools calmd -b --no-PG ~{prefix}.pre.bam ~{ref_fasta} > ~{prefix}.bam
        samtools index -@${NUM_CPUS} ~{prefix}.bam
    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             mem,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
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

# A simple task to covert SAM-formatted alignment to PAF format
task SAMtoPAF {
    input {
        File sam_formatted_file
        File? index

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        sam_formatted_file: "SAM-formated input file to be converted to PAF (note currently we only support SAM or BAM, not CRAM)"
        index:              "[optional] index for sam_formatted_file"
    }

    String prefix = basename(basename(sam_formatted_file, ".bam"), ".sam") # we have hack like this because WDL stdlib doesn't provide endsWith stuff

    Int disk_size = 2*ceil(size(sam_formatted_file, "GB"))

    command <<<
        set -eu

        MM2_VERSION="2.24"

        filename=$(basename -- ~{sam_formatted_file})
        extension="${filename##*.}"
        if [[ "$extension" == "sam" ]]; then
            /minimap2-${MM2_VERSION}_x64-linux/k8 \
                /minimap2-${MM2_VERSION}_x64-linux/paftools.js \
                sam2paf \
                -L \
                ~{sam_formatted_file} \
                > ~{prefix}".paf"
        elif [[ "$extension" == "bam" ]]; then
            samtools view -h ~{sam_formatted_file} | \
                /minimap2-${MM2_VERSION}_x64-linux/k8 \
                /minimap2-${MM2_VERSION}_x64-linux/paftools.js \
                sam2paf \
                -L \
                - \
                > ~{prefix}".paf"
        else
            echo "Currently we only support SAM or BAM (not CRAM)." && exit 1;
        fi
    >>>

    output {
        File pat_formatted_file = "~{prefix}.paf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
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
