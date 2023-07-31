version 1.0

import "../../structs/Structs.wdl"

task Minimap2 {
    input {
        Array[File] reads
        Array[String] reads_file_basenames
        File ref_fasta

        String map_preset

        String? RG
        String? library

        Array[String] tags_to_preserve = []

        String prefix = "out"
        String disk_type = "SSD"
        RuntimeAttr? runtime_attr_override
    }
    meta {
        descrpiton: "A wrapper to minimap2 for mapping & aligning (groups of) sequences to a reference. Note that this only works for reads belonging to a single readgroup."
    }
    parameter_meta {
        reads: {
            desciption: "query sequences to be mapped and aligned",
            localization_optional: true
        }
        reads_file_basenames: "basenames of the BAM files, note this includes the extention (e.g. .bam, .fasta.gz, etc) of files"
        ref_fasta:        "reference fasta"
        RG:               "read group information to be supplied to parameter '-R' (note that tabs should be input as '\t'); if the input uBAM files contain the RG information, this overrides that; we recommend providing this when input isn't a uBAM with RG information."
        map_preset:       "preset to be used for minimap2 parameter '-x'"
        tags_to_preserve: "SAM tags to carry over to aligned bam file; ignored when inputs are not SAM/BAM/CRAM files"
        prefix:           "[default-valued] prefix for output BAM"
    }

    Boolean is_custom_rg = defined(RG)

    Boolean fix_library_entry = if defined(library) then true else false

    Int pd_disk_size = 1 + 10*ceil(size(reads, "GiB") + ceil(size(ref_fasta, "GiB")))
    Int local_disk_size = if(size(reads, "GiB")>150) then 750 else 375
    Int disk_size = if('LOCAL'==disk_type) then local_disk_size else pd_disk_size

    Boolean do_preserve_tags = if length(tags_to_preserve) != 0 then true else false

    # we limit the number of CPUs here because the number of input files could be huge (e.g. ONT fastqs)
    # also, the task is mostly CPU-bound (i.e. minimap2)
    Int max_cpus = 96
    Int desired_cpus = (if('LOCAL'==disk_type) then 32 else 24) * length(reads)
    Int cpus = if max_cpus < desired_cpus then max_cpus else desired_cpus  # WDL 1.0 doesn't have a max(,)....
    Int mem = cpus * (if('LOCAL'==disk_type) then 6 else 5)

    Int mm2_threads = cpus - 2

    command <<<
        set -euxo pipefail

        FILE="~{reads[0]}"

        ############
        # localize data (much faster than Cromwell)
        mkdir -p reads_dir
        time gcloud storage cp ~{sep=' ' reads} /cromwell_root/reads_dir/
        ls reads_dir

        cd reads_dir

        ############
        # parameter setting
        NUM_CPUS=$( cat /proc/cpuinfo | grep '^processor' | tail -n1 | awk '{print $NF+1}' )
        RAM_IN_GB=$( free -g | grep "^Mem" | awk '{print $2}' )
        echo "Memory info:"
        cat /proc/meminfo
        echo ""

        # mimimap2
        # todo: if input is SAM, take RG from the SAM file when custom RG isn't specified from the WDL input
        # rg_len=$(echo -n '~{RG}' | wc -c | awk '{print $NF}')
        # if [[ $rg_len -ne 0 ]] ; then
        if ~{is_custom_rg}; then
            MAP_PARAMS="-ayYL --MD --eqx --cs -x ~{map_preset} -t ~{mm2_threads} -K4G -R ~{RG}      ~{ref_fasta}"
        elif [[ "$FILE" =~ \.bam$ ]]; then
            ubam_rg=$(samtools view -H "$FILE" | grep "^@RG")
            MAP_PARAMS="-ayYL --MD --eqx --cs -x ~{map_preset} -t ~{mm2_threads} -K4G -R ${ubam_rg} ~{ref_fasta}"
        else
            echo "Warning! Input isn't an uBAM, but readgroup line isn't provided. This isn't recommended. We'll continue, though."
            MAP_PARAMS="-ayYL --MD --eqx --cs -x ~{map_preset} -t ~{mm2_threads} -K4G               ~{ref_fasta}"
        fi

        # samtools sort
        SORT_PARAMS="-@4 -m4G --no-PG -o ~{prefix}.pre.bam"

        ############
        # minimap2
        # We write to a SAM file before sorting and indexing because rarely, doing everything
        # in a single one-liner leads to a truncated file error and segmentation fault of unknown
        # origin.  Separating these commands requires more resources, but is more reliable overall.

        if [[ "$FILE" =~ \.bam$ ]]; then

            date
            idx=0
            for f in "~{sep=' ' reads_file_basenames}" ; do
                samtools fastq ~{true=' ' false='-t' is_custom_rg} \
                    ~{true='-T ' false=' ' do_preserve_tags } ~{sep=',' tags_to_preserve} \
                    "$f" \
                | minimap2 ${MAP_PARAMS} - \
                | samtools sort -@"~{mm2_threads}" -m4G --no-PG -o "to-merge.${idx}.bam" -
                idx=$(( idx + 1 ))
            done
            date

            if [[ ${idx} == 1 ]]; then
                mv to-merge.0.bam ~{prefix}.pre.bam
            else
                time \
                samtools merge -c -@1 --no-PG -o ~{prefix}.pre.bam to-merge.*.bam;
            fi
        elif [[ "$FILE" =~ \.fastq$ ]] || [[ "$FILE" =~ \.fq$ ]]; then
            cat ~{sep=' ' reads_file_basenames} | minimap2 $MAP_PARAMS - > tmp.sam
        elif [[ "$FILE" =~ \.fastq.gz$ ]] || [[ "$FILE" =~ \.fq.gz$ ]]; then
            zcat ~{sep=' ' reads_file_basenames} | minimap2 $MAP_PARAMS - > tmp.sam
        elif [[ "$FILE" =~ \.fasta$ ]] || [[ "$FILE" =~ \.fa$ ]]; then
            cat ~{sep=' ' reads_file_basenames} | python3 /usr/local/bin/cat_as_fastq.py | minimap2 $MAP_PARAMS - > tmp.sam
        elif [[ "$FILE" =~ \.fasta.gz$ ]] || [[ "$FILE" =~ \.fa.gz$ ]]; then
            zcat ~{sep=' ' reads_file_basenames} | python3 /usr/local/bin/cat_as_fastq.py | minimap2 $MAP_PARAMS - > tmp.sam
        else
            echo "Did not understand file format for '$FILE'"
            exit 1
        fi

        ############
        # sort, necessary only when input format isn't BAM
        if [[ -f tmp.sam ]]; then
            time \
            samtools sort "${SORT_PARAMS}" tmp.sam
            rm tmp.sam
        fi

        ############
        # post-fix
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

            time \
            samtools reheader -@"${NUM_CPUS}" fixed_header.txt ~{prefix}.pre.tmp.bam > ~{prefix}.pre.bam
            rm ~{prefix}.pre.tmp.bam
        fi

        ############
        # MD tag and index
        time \
        samtools calmd -@"${NUM_CPUS}" -b ~{prefix}.pre.bam ~{ref_fasta} > ~{prefix}.bam
        time \
        samtools index -@"${NUM_CPUS}" ~{prefix}.bam

        ############
        # move results up
        mv ~{prefix}.bam ~{prefix}.bam.bai /cromwell_root
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
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-minimap2:2.26-gcloud"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SAMtoPAF {
    input {
        File sam_formatted_file
        File? index

        RuntimeAttr? runtime_attr_override
    }
    meta {
        description: "Convert SAM-formatted alignment to PAF format using minimap2's paftools.js"
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
