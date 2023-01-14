version 1.0

import "../Structs.wdl"

task GetReadGroupInfo {
    meta {
        desciption:
        "Get some read group information given a single-readgroup BAM. If the requested keys are absent, a null value is assigned in the returned entry. If the BAM contains multiple read groups, results are undetermined."
    }

    input {
        String uBAM  # not using file as call-caching brings not much benefit

        Array[String] keys
        String null_value_representation = "None"
    }

    parameter_meta {
        keys: "A list of requested fields in the RG line, e.g. ID, SM, LB."
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{uBAM} | grep "^@RG" > one_rg_per_line.txt
        num_rgs=$(wc -l one_rg_per_line.txt | awk '{pritn $1}')
        if [[ num_rgs -gt 1 ]]; then exit 1; fi

        cat one_rg_per_line.txt | tr '\t' '\n' > rh_header.txt

        for attribute in ~{sep=' ' keys}; do
            if grep -q "^${attribute}" rh_header.txt; then
                value=$(grep "^${attribute}" rh_header.txt | awk -F ':' '{print $2}')
            else
                value="~{null_value_representation}"
            fi
            echo -e "${attribute}\t${value}" >> "result.txt"
        done
    >>>

    output {
        Map[String, String] read_group_info = read_map("result.txt")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task SplitByRG {
    meta {
        description: "Split a BAM file that was aggregated, for the same sample, into pieces by read group."
    }
    input {
        File bam

        String out_prefix

        Int? num_ssds

        Boolean retain_rgless_records = false
        Boolean sort_and_index = false
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: "BAM to be split"
        out_prefix: "prefix for output bam and bai file names"
        sort_and_index: "if the user wants to (pos-)sort and index the resulting BAMs; this indicates the input BAM is mapped"

        split_bam: "the resuling BAMs, each having reads only in a single read group"
        split_bai: "the accompanying BAIs, if possible and explicit requested"
    }

    Int disk_size = if defined(num_ssds) then 375*select_first([num_ssds]) else 1+3*ceil(size([bam], "GB"))

    Array[String] extra_args = if (retain_rgless_records) then ["-u", "~{out_prefix}_noRG.bam"] else [""]
    command <<<
        set -eux

        samtools view -H ~{bam} | grep "^@RG" > "read_groups_header.txt"
        cat "read_groups_header.txt" | tr '\t' '\n' | grep "^ID:"  | awk -F ':' '{print $2}' > "RG_ids.txt"

        samtools split -@3 \
            -f "~{out_prefix}_%#.bam" \
            ~{sep=" " extra_args} \
            ~{bam}
        if ~{sort_and_index} ;
        then
            # cleanup space for the sorting
            rm ~{bam}
            for split_bam in "~{out_prefix}_"*.bam;
            do
                mv "${split_bam}" temp.bam
                samtools sort \
                    --write-index \
                    -o "${split_bam}##idx##${split_bam}.bai" \
                    temp.bam
            done
        fi
    >>>

    output {
        File read_group_header = "read_groups_header.txt"
        Array[String] rg_ids   = read_lines("RG_ids.txt")
        Array[File]  split_bam = glob("*.bam")
        Array[File?] split_bai = glob("*.bai")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task GatherBamMetadata {
    meta {
        description: "Check several metadata of an input BAM (aliged? sort order? etc)"
    }
    input {
        File bam
    }
    parameter_meta {
        bam: {
            desciption: "BAM to be checked",
            localization_optional: true
        }
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} > header.txt

        grep -F "@HD" header.txt | tr '\t' '\n' > hd.lines.txt
        if grep -q "SO" hd.lines.txt;
        then
            echo "true" > "is_sorted.txt"
            grep "SO" hd.lines.txt | awk -F ':' '{print $2}' > "sort_order.txt"
        else
            echo "false" > "is_sorted.txt"
            echo "NA" > "sort_order.txt"
        fi

        # we use two conditions: @SQ lines in header, and at least some mapped reads
        if grep -q "@SQ" "header.txt";
        then
            export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
            if [[ $(samtools view -F 4 ~{bam} | head | wc -l | awk '{print $1}') -ge 1 ]];
            then
                echo "true" > "is_mapped.txt"
            else
                echo "unsure" > "is_mapped.txt" # this will trigger error in WDL later, but that's intentional because we cannot be sure
            fi
        else
            echo "false" > "is_mapped.txt"
        fi
    >>>

    output {
        Boolean is_aligned = read_boolean("is_mapped.txt")

        Boolean is_sorted = read_boolean("is_sorted.txt")
        String sort_order = read_string("sort_order.txt")
    }

    runtime {
        cpu:    2
        memory: "8 GiB"
        disks:  "local-disk 20 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
    }
}

task UnAlignBam {
    meta {
        desciption:
        "Go from an aligned BAM to unaligned BAM"
    }
    input {
        File bam

        String out_prefix

        Int? num_ssds
        RuntimeAttr? runtime_attr_override
    }
    Int disk_size = if defined(num_ssds) then 375*select_first([num_ssds]) else 1+3*ceil(size([bam], "GB"))

    String bam_emptyness = "bam_is_empty.txt"

    command <<<
        set -eux

        samtools fastq ~{bam} | gzip > ~{out_prefix}.fq.gz

        # need to reheader to keep history accumulated in the original bam, except those reference contigs
        samtools view -H ~{bam} | grep -v "@SQ" > original_header.txt
        cat original_header.txt

        cnt=$(samtools view -c ~{bam})
        if [[ ${cnt} -eq 0 ]]; then # if input bam is empty, then just output header and signal in the output
            touch "~{out_prefix}.unaligned.bam"
            samtools view -h -o "~{out_prefix}.unaligned.bam" original_header.txt
            echo "true" > ~{bam_emptyness}
            exit 0
        fi
        time \
        samtools import \
            -0 ~{out_prefix}.fq.gz \
            -o ~{out_prefix}.unaligned.tbrh.bam

        time \
        samtools reheader \
            original_header.txt \
            ~{out_prefix}.unaligned.tbrh.bam \
            > ~{out_prefix}.unaligned.bam
        echo "false" > ~{bam_emptyness}
    >>>

    output {
        File fq   = "~{out_prefix}.fq.gz"
        File uBAM = "~{out_prefix}.unaligned.bam"
        Boolean is_bam_empty = read_boolean(bam_emptyness)
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SamtoolsReset {
    meta {
        description: "Use samtools reset to drop alignment information from the input bam."
    }

    parameter_meta {
        bam: "aligned BAM to operate on"
        addtional_tags_to_drop: "tags in each alignment record to be dropped; usually these are tags produced by the mapper/aligner that generated the original alignment"
    }
    input {
        File bam
        # these are known mm2 tags and pbmm2 tags
        Array[String] addtional_tags_to_drop = ['cg', 'cm', 'cs',
                                                'de', 'dv',
                                                'ms',
                                                'nn',
                                                'rl',
                                                's1', 's2',
                                                'ts',
                                                'mc', 'mg', 'mi', 'rm']

        Int? num_ssds
        RuntimeAttr? runtime_attr_override
    }

    output {
        File res = "~{prefix}.unaligned.bam"
        File original_sam_flag_stats = "orignal.SAM-flag.stats.txt"
    }

    Array[String] std_tags_to_drop = ['MD', 'NM', 'AS', 'AS', 'XA']
    Array[String] tags_to_drop = flatten([std_tags_to_drop, addtional_tags_to_drop])

    String prefix = basename(bam, ".bam")

    Int disk_size = if defined(num_ssds) then 375*select_first([num_ssds]) else 1+10*ceil(size([bam], "GB"))
    String disk_type = if defined(num_ssds) then " LOCAL" else " SSD"
    command <<<
        set -eux

        samtools view -@1 ~{bam} | grep -v "^@" | awk -F '\t' '{print $2}' | sort | uniq -c > orignal.SAM-flag.stats.txt &

        samtools reset -@3 \
            --remove-tag ~{sep=',' tags_to_drop} \
            -o ~{prefix}.unaligned.bam \
            ~{bam}
        wait
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-bam-utils:0.1.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + disk_type
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Drop2304Alignments {
    meta {
        description: "Drop secondary and supplementary alignments from a BAM"
    }

    parameter_meta {
        bam: "input BAM"
        filtered_bam: "output BAM without secondary and supplementary alignment records"
        num_ssds: "Number of local SSDs to use; if not provided, will use SSD persistent disks"
    }
    input {
        File bam
        File bai
        Int? num_ssds
        RuntimeAttr? runtime_attr_override
    }

    output {
        File filtered_bam = "~{prefix}.no_SAM_2304.bam"
        File filtered_bai = "~{prefix}.no_SAM_2304.bam.bai"
    }

    String prefix = basename(bam, ".bam")

    Int disk_size = if defined(num_ssds) then 375*select_first([num_ssds]) else 1+3*ceil(size([bam], "GB"))
    String disk_type = if defined(num_ssds) then " LOCAL" else " SSD"

    command <<<
        set -eux

        samtools view -h -@3 \
            -F 2304 \
            --write-index \
            -o "~{prefix}.no_SAM_2304.bam##idx##~{prefix}.no_SAM_2304.bam.bai" \
            ~{bam}
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + disk_type
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SortBamByQName {
    meta {
        description: "QueryName-sort a BAM"
    }

    parameter_meta {
        bam: "input BAM"
        qnsort_bam: "output BAM sorted by query name"
        num_ssds: "Number of local SSDs to use; if not provided, will use SSD persistent disks"

        multi_record_queries: "File holding names of queries that has multiple records in the output"
    }
    input {
        File bam
        Int? num_ssds
        RuntimeAttr? runtime_attr_override
    }

    output {
        File qnsort_bam = "~{prefix}.qname-sorted.bam"
        File multi_record_queries = "multi_record_queries.txt"
    }

    String prefix = basename(bam, ".bam")

    Int disk_size = if defined(num_ssds) then 375*select_first([num_ssds]) else 1+4*ceil(size([bam], "GB"))
    String disk_type = if defined(num_ssds) then " LOCAL" else " SSD"

    command <<<
        set -eux

        samtools view -H ~{bam} | grep "^@HD" | tr '\t' '\n' > hd.line.txt
        if grep -q 'SO:queryname' hd.line.txt;
        then
            echo "already sorted"
            mv ~{bam} "~{prefix}.qname-sorted.bam"
            exit 0
        fi
        samtools sort -@3 -m2G \
            -n \
            -O BAM \
            ~{bam} \
        > "~{prefix}.qname-sorted.bam"

        touch multi_record_queries.txt
        samtools view "~{prefix}.qname-sorted.bam" | awk '{print $1}' | uniq -d > multi_record_queries.txt
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + disk_type
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MergeBams {
    input {
        Array[File] bams
        Boolean input_is_aligned
        String prefix = "out"

        Int? num_ssds
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bams:   "input array of BAMs to be merged"
        prefix: "[default-valued] prefix for output BAM"
    }

    Int disk_size = if defined(num_ssds) then 375*select_first([num_ssds]) else 1+3*ceil(size(bams, "GB"))
    String disk_type = if defined(num_ssds) then " LOCAL" else " SSD"

    command <<<
        set -euxo pipefail

        if [[ ~{input_is_aligned} ]]; then
            samtools merge \
                -p -c --no-PG \
                -@ 2 \
                --write-index \
                -o "~{prefix}.bam##idx##~{prefix}.bam.bai" \
                ~{sep=" " bams}
        else
            samtools merge \
                -p -c --no-PG \
                -@ 2 \
                -o "~{prefix}.bam" \
                ~{sep=" " bams}
        fi
    >>>

    output {
        File  merged_bam = "~{prefix}.bam"
        File? merged_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + disk_type
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ValidateSamFile {
    meta {
        desciption: "Call GATK/Picard ValidateSamFile to validate input BAM: https://bit.ly/3JMutxp."
    }
    parameter_meta {
        validation_mode: "Desired valiation mode"
        disk_type: "Type of disk to use for the computation."
    }

    input {
        File bam
        String validation_mode = "SUMMARY"

        Array[String] validation_errs_to_ignore = ["INVALID_TAG_NM",  # for the purpose we currently have, NM and CIGAR don't matter, and longreads have no mates
                                                    "MISSING_TAG_NM",
                                                    "INVALID_CIGAR",
                                                    "ADJACENT_INDEL_IN_CIGAR",
                                                    "CIGAR_MAPS_OFF_REFERENCE",
                                                    "MISMATCH_MATE_CIGAR_STRING",
                                                    "MATE_CIGAR_STRING_INVALID_PRESENCE",
                                                    "MATE_NOT_FOUND",
                                                    "INVALID_MAPPING_QUALITY",
                                                    "INVALID_FLAG_MATE_UNMAPPED",
                                                    "MISMATCH_FLAG_MATE_UNMAPPED",
                                                    "INVALID_FLAG_MATE_NEG_STRAND",
                                                    "MISMATCH_FLAG_MATE_NEG_STRAND",
                                                    "INVALID_MATE_REF_INDEX",
                                                    "MISMATCH_MATE_REF_INDEX",
                                                    "MISMATCH_MATE_ALIGNMENT_START",
                                                    "MATE_FIELD_MISMATCH",
                                                    "PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND"
                                                   ]

        String disk_type = "LOCAL"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(bam, "GiB")) + 50
    String output_basename = basename(basename(bam, ".bam"), ".cram")
    String output_name = "${output_basename}_${validation_mode}.txt"

    command <<<
        set -eux

        gatk ValidateSamFile \
            --INPUT ~{bam} \
            --OUTPUT ~{output_name} \
            --MODE ~{validation_mode} \
            ~{true="--IGNORE " false="" 0<length(validation_errs_to_ignore)} \
            ~{sep=" --IGNORE " validation_errs_to_ignore}

    >>>

    output {
        File validation_report = "${output_name}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.4.0.0"
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

task CountReadGroups {
    meta {
        desciption: "Count the number of RG lines in the header of the BAM file."
    }
    input {
        String bam  # not using file as call-caching brings not much benefit
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} | grep -c "^@RG" > "rg_cnt.txt"
    >>>

    output {
        Int num_rg = read_int("rg_cnt.txt")
    }
    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}



task SplitQNameSortedNo2304Bam {
    meta {
        description: "Split a BAM into pieces. The input BAM is assumed to have no supplementary/secondary alignments, and query-name sorted."
    }

    parameter_meta {
        bam: "input BAM that's assumed to be query-name sorted and without supplementary and secondary alignments"
        num_lines_per_chunk: "number to lines in each shard"
        split_bams: "output BAM without secondary and supplementary alignment records"

        num_ssds: "Number of local SSDs to use; if not provided, will use SSD persistent disks"
    }
    input {
        File bam
        Int num_lines_per_chunk
        Int? num_ssds
        RuntimeAttr? runtime_attr_override
    }

    output {
        Array[File] split_bams = glob("~{prefix}_split-*.bam")
    }

    String prefix = basename(bam, ".bam")

    Int disk_size = if defined(num_ssds) then 375*select_first([num_ssds]) else 1+10*ceil(size([bam], "GB"))
    String disk_type = if defined(num_ssds) then " LOCAL" else " SSD"


    # remove commented out codes, use --no-PG, add timestamps
    command <<<
        set -eux

        samtools view -H ~{bam} > original.sam.header

        #     > tmp.sam.txt.gz
        # ls -lh tmp.sam.txt.gz && rm ~{bam}

        # zcat tmp.sam.txt.gz \
        samtools view -O SAM ~{bam} \
            | grep -v "^@" \
        | split \
            -d \
            -l "~{num_lines_per_chunk}" \
            - \
            "~{prefix}_split-"
        # rm tmp.sam.txt.gz

        for ff in `ls "~{prefix}_split-"*`;
        do
            cat \
                original.sam.header \
                "${ff}" \
            | samtools view --no-PG -bh \
                -o "${ff}.bam" -
            rm "${ff}"
        done
        wait
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + disk_type
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task DropAlignmentInfo {
    meta {
        description: ""
    }

    parameter_meta {
        bam: "aligned BAM to operate on"
        addtional_tags_to_drop: "tags in each alignment record to be dropped; usually these are tags produced by the mapper/aligner that generated the original alignment"
    }
    input {
        File bam
        # these are known mm2 tags and pbmm2 tags
        Array[String] addtional_tags_to_drop = ['cg', 'cm', 'cs',
                                                'de', 'dv',
                                                'ms',
                                                'nn',
                                                'rl',
                                                's1', 's2',
                                                'ts',
                                                'mc', 'mg', 'mi', 'rm']

        Int? num_ssds
        RuntimeAttr? runtime_attr_override
    }

    output {
        File res = "~{prefix}.unaligned.bam"
        File original_sam_flag_stats = "orignal.SAM-flag.stats.txt"
    }

    String docker_version = "0.1.0"

    Array[String] std_tags_to_drop = ['MD', 'NM', 'AS', 'AS', 'XA']
    Array[String] tags_to_drop = flatten([std_tags_to_drop, addtional_tags_to_drop])

    String prefix = basename(bam, ".bam")

    Int disk_size = if defined(num_ssds) then 375*select_first([num_ssds]) else 1+10*ceil(size([bam], "GB"))
    String disk_type = if defined(num_ssds) then " LOCAL" else " SSD"
    command <<<
        set -eux

        # must save header first
        du -sh ~{bam}
        samtools view --no-PG -H ~{bam} | grep "^@SQ" > original.header.without.SQ.txt
        samtools view -@1 ~{bam} | grep -v "^@" | awk -F '\t' '{print $2}' | sort | uniq -c > orignal.SAM-flag.stats.txt &

        #####################
        # first, reset columns 3-9.
        # col-1 is QName, col-2 is SAM-Flag, col-10 is the seq, col-11 is the baseQ, on-wards are tags. Refer to the SAM-spec on the other cols.

        # make sure we record this and the next operation for documentation
        cat \
            original.header.without.SQ.txt \
            <(echo -e "@PG\tID:lrma-pipe-DropAlignmentInfo\tPN:DropAlignmentInfo\tVN:~{docker_version}\tDS:drop alignment information by resetting cols 3-9, and restoring original bases and qualities") \
            | gzip \
        > header.to.use.txt.gz

        # refer to the SAM-spec on cols-3~9.
        samtools view --no-PG ~{bam} \
            | grep -v "^@" \
            | awk -F '\t' 'BEGIN{OFS="\t"} {$3="*"; $4=0; $5=255; $6="*"; $7="*"; $8=0; $9=0; print}' \
            | gzip \
        > col-3-9_reset.txt.gz

        # wait before deleting
        wait && rm ~{bam}
        du -sh col-3-9_reset.txt.gz &

        zcat \
            header.to.use.txt.gz \
            col-3-9_reset.txt.gz \
        > col-3-9_reset.sam

        # wait before deleting
        wait && rm col-3-9_reset.txt.gz
        #####################
        # second, for records where SAM flag 16 is turned on (reverse strand), RC col-10, and reverse col-11, then reset col2 to 4 (unmapped)
        # this is critical to avoid potentially giving strand-bias induced by the previous mapping step
        # here the assumption is made that only primary alignment is kept (hence no hard-clipping)
        python3 /opt/reverse_seq_and_qual.py \
            -i col-3-9_reset.sam \
            -o col-3-9_reset_orig-baseQ-restored.bam
        rm col-3-9_reset.sam

        #####################
        # lastly, drop tags
        samtools view -bh -@3 \
            --remove-tag ~{sep=',' tags_to_drop} \
            -o ~{prefix}.unaligned.bam \
            col-3-9_reset_orig-baseQ-restored.bam
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-bam-utils:" + docker_version
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + disk_type
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}