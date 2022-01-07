version 1.0

import "Structs.wdl"

workflow AlignReads {
    meta {
        description : "This workflow aligns reads using minimap2."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }
    input {
        Array[File] reads
        File ref_fasta

        String map_preset

        String RG = ""

        String prefix = "out"
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:      "query sequences to be mapped and aligned"
        ref_fasta:  "reference fasta"
        RG:         "[optional] read group information to be supplied to parameter '-R' (note that tabs should be input as '\t')"
        map_preset: "preset to be used for minimap2 parameter '-x'"
        prefix:     "[default-valued] prefix for output BAM"
    }

    # Call our alignment task:
    call Minimap2 {
        input:
            reads = reads,
            ref_fasta = ref_fasta,
            map_preset = map_preset,
            RG = RG,
            prefix = prefix,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File aligned_bam = Minimap2.aligned_bam
        File aligned_bai = Minimap2.aligned_bai
    }
}

# A wrapper to minimap2 for mapping & aligning (groups of) sequences to a reference
task Minimap2 {
    input {
        Array[File] reads
        File ref_fasta

        String map_preset

        String RG = ""

        String prefix = "out"
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:      "query sequences to be mapped and aligned"
        ref_fasta:  "reference fasta"
        RG:         "[optional] read group information to be supplied to parameter '-R' (note that tabs should be input as '\t')"
        map_preset: "preset to be used for minimap2 parameter '-x'"
        prefix:     "[default-valued] prefix for output BAM"
    }

    # 10x for the decompressed file size
    # 2x for potential for FASTQ and SAM files (from file conversion).
    # 2x for extra "just in case" space.
    # +1 to handle small files
    Int disk_size = 1 + 10*2*2*ceil(size(reads, "GB") + size(ref_fasta, "GB"))

    Int cpus = 4
    Int mem = 30

    # This is a hack to fix the WDL parsing of ${} variables:
    String DOLLAR = "$"
    command <<<
        set -euxo pipefail

        rg_len=$(echo -n '~{RG}' | wc -c | awk '{print $NF}')
        if [[ $rg_len -ne 0 ]] ; then
            # Sometimes we have to sanitize our read groups:
            sanitized_read_group=$( echo "~{RG}" | sed -e 's# .*##g' | sed 's#\t.*##g' )

            echo "Original Read Group: ~{RG}"
            echo "Sanitized Read Group: $sanitized_read_group"

            MAP_PARAMS="-ayYL --MD --eqx -x ~{map_preset} -R $sanitized_read_group -t ~{cpus} ~{ref_fasta}"
        else
            MAP_PARAMS="-ayYL --MD --eqx -x ~{map_preset} -t ~{cpus} ~{ref_fasta}"
        fi
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
                samtools fastq "$f"
            done > tmp.fastq

            echo "Memory info:"
            cat /proc/meminfo
            echo ""

            minimap2 $MAP_PARAMS tmp.fastq > tmp.sam
        else
            echo "Did not understand file format for '$FILE'"
            exit 1
        fi

        samtools sort -@~{cpus} -m~{mem}G --no-PG -o ~{prefix}.bam tmp.sam
        samtools index ~{prefix}.bam
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

        filename=$(basename -- ~{sam_formatted_file})
        extension="${filename##*.}"
        if [[ "$extension" == "sam" ]]; then
            /minimap2-2.24_x64-linux/k8 \
                /minimap2-2.24_x64-linux/paftools.js \
                sam2paf \
                -L \
                ~{sam_formatted_file} \
                > ~{prefix}".paf"
        elif [[ "$extension" == "bam" ]]; then
            samtools view -h ~{sam_formatted_file} | \
                /minimap2-2.24_x64-linux/k8 \
                /minimap2-2.24_x64-linux/paftools.js \
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
