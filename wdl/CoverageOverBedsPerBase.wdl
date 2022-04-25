version 1.0

import "tasks/Finalize.wdl" as FF

workflow CoverageOverBedsPerBase {
    input {
        File bam
        File bai
    	File four_col_bed
        String participant_name

        Int? min_base_q
        Int? min_map_q

    	Boolean fail_on_non_4_col_bed

        String workflow_root_name
        String gcs_out_root_dir
    }

    String dir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_root_name}/~{participant_name}/alignments"

    # call PadAndSplit {input: four_col_bed = four_col_bed, fail_on_non_4_col_bed = fail_on_non_4_col_bed, split_line_cnt = 500}

    # scatter (bed in PadAndSplit.splitted_padded_beds) {
    #     call ViewDepth {input: bam = bam, bai = bai, bed = bed, min_base_q = min_base_q, min_map_q = min_map_q}
    # }

    # call MergeDepth {input: perbase_depth_split = ViewDepth.subset_bed_depth, out_prefix = basename(bam, ".bam") + '.' + basename(four_col_bed, ".bed")}
    # call FF.FinalizeToFile { input: outdir = dir, file = MergeDepth.merged_sorted }

    call CheckAndSplit {input: four_col_bed = four_col_bed, fail_on_non_4_col_bed = fail_on_non_4_col_bed, split_line_cnt = 500}
    scatter (bed in CheckAndSplit.splitted_beds) {
        call PaddedViewExactDepth {input: bam = bam, bai = bai, bed = bed, min_base_q = min_base_q, min_map_q = min_map_q}
    }

    call MergeDepth {input: perbase_depth_split = PaddedViewExactDepth.subset_bed_depth, out_prefix = basename(bam, ".bam") + '.' + basename(four_col_bed, ".bed")}
    call FF.FinalizeToFile { input: outdir = dir, file = MergeDepth.merged_sorted }

    output {
        File depth = FinalizeToFile.gcs_path
    }
}

###############################################
task SplitCollectMerge {
    input {
    	File bam
        File bai
    	File four_col_bed
        Int? min_base_q
        Int? min_map_q
    	Boolean fail_on_non_4_col_bed
    }

    parameter_meta {
        bam: {
            localization_optional: true
        }
    }

    String bq_filter = if defined(min_base_q) then " -q ~{min_base_q} " else " "
    String mq_filter = if defined(min_map_q)  then " -Q ~{min_map_q} " else " "

    String prefix = basename(bam, ".bam") + '.' + basename(four_col_bed, ".bed")

    command <<<
        set -euxo pipefail


        #######################################################################
        if awk -F '\t' 'NF!=4{exit 1}' ~{four_col_bed} ; then
            echo;
        else
            if ~{fail_on_non_4_col_bed}; then
                echo "input BED must be exactly 4 columns" && exit 1;
            else
                echo "WARNING!!! Input BED has more than 4 columns, information held in columns beyond the 4th will be dropped."
            fi
        fi
        cut -f -4 ~{four_col_bed} > slim.bed

        # note --additional-suffix= is supported only after ≥ 8.16 of gnu split
        split -d -l 500 slim.bed --additional-suffix=".bed" splitted_padded_

        for chunk in splitted_padded_*.bed; do
            idx=$(echo "${chunk}" | awk -F '.' '{print $1}' | awk -F '_' '{print $3}')
            touch "temp.${idx}.tsv"

            set +x
            echo "${idx}"
            export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
            # Note it's done this way due to https://github.com/samtools/samtools/issues/1581
            while IFS=$'\t' read -r chr column2 end bed_id; do
                echo "${bed_id}"
                beg=$((column2+1))
                arr=$(samtools depth -aa -J -r "${chr}:${beg}-${end}" ~{bq_filter} ~{mq_filter} ~{bam} | awk '{print $3}' | paste -s -d, -)
                echo -e "${chr}\t${column2}\t${end}\t${bed_id}\t${arr}" >> "temp.${idx}.tsv"
            done < "${chunk}"
        done
        set -x
        cat temp.*.tsv > temp.merged.tsv
        out="~{prefix}.samtools.depth.txt"
        echo -e "chr\tBED_beg\tBED_end\tid\tdepth_array_aa.J" > "${out}"
        cat temp.merged.tsv >> "${out}"

        gzip -k "${out}"
    >>>

    output {
        File depth = "~{prefix}.samtools.depth.txt.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            375,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
    }
    RuntimeAttr runtime_attr = select_first([default_attr])
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

task ViewDepth {

    meta {
        description: "For subsetting a high-coverage BAM stored in GCS, without localizing (more resilient to auth. expiration)."
    }

    input {
        File bam
        File bai

        File bed
        Int? min_base_q
        Int? min_map_q

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: {
            localization_optional: true
        }
    }

    Int disk_size = ceil(size([bam, bai], "GB"))

    String subset_prefix = basename(bam, ".bam") + basename(bed, ".bed")

    String bq_filter = if defined(min_base_q) then " -q ~{min_base_q} " else " "
    String mq_filter = if defined(min_map_q)  then " -Q ~{min_map_q} " else " "

    command <<<

        # Note it's done this way due to https://github.com/samtools/samtools/issues/1581
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        set -euxo pipefail
        samtools view \
            -bhX -@ 1 \
            --verbosity=8 \
            --regions-file ~{bed} \
            --write-index \
            -o "~{subset_prefix}.bam##idx##~{subset_prefix}.bam.bai" \
            "~{bam}" "~{bai}" &
        set +x
        sleep 1200 && export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        set -x
        wait

        while IFS=$'\t' read -r chr column2 end bed_id; do
            echo "${bed_id}"
            beg=$((column2+1))
            arr=$(samtools depth -aa -J -r "${chr}:${beg}-${end}" ~{bq_filter} ~{mq_filter} "~{subset_prefix}.bam" | awk '{print $3}' | paste -s -d, -)
            echo -e "${chr}\t${column2}\t${end}\t${bed_id}\t${arr}" >> "~{subset_prefix}.tsv"
        done < "~{bed}"

        gzip -k "~{subset_prefix}.tsv"
    >>>

    output {
        File subset_bed_depth = "~{subset_prefix}.tsv.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
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

###############################################

task CheckAndSplit {
    input {
        File four_col_bed
        Boolean fail_on_non_4_col_bed
        Int split_line_cnt
    }

    command <<<

        set -eu
        if awk -F '\t' 'NF!=4{exit 1}' ~{four_col_bed} ; then
            echo;
        else
            if ~{fail_on_non_4_col_bed}; then
                echo "input BED must be exactly 4 columns" && exit 1;
            else
                echo "WARNING!!! Input BED has more than 4 columns, information held in columns beyond the 4th will be dropped."
            fi
        fi

        split -d -l ~{split_line_cnt} ~{four_col_bed} --additional-suffix=".bed" splitted_
    >>>

    output {
        Array[File] splitted_beds = glob("splitted_*.bed")
    }
    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker:         "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task PaddedViewExactDepth {

    meta {
        description: "For subsetting a high-coverage BAM stored in GCS, without localizing (more resilient to auth. expiration)."
    }

    input {
        File bam
        File bai

        File bed
        Int? min_base_q
        Int? min_map_q

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: {
            localization_optional: true
        }
    }

    Int disk_size = ceil(size([bam, bai], "GB"))

    String subset_prefix = basename(bam, ".bam") + basename(bed, ".bed")

    String bq_filter = if defined(min_base_q) then " -q ~{min_base_q} " else " "
    String mq_filter = if defined(min_map_q)  then " -Q ~{min_map_q} " else " "

    command <<<

        # Note it's done this way due to https://github.com/samtools/samtools/issues/1581
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        set -euxo pipefail

        awk 'BEGIN{OFS="\t"} {a=$2-100; b=$3+100; if(a<0) a=0; print $1, a, b, $4}' ~{bed} > padded.bed
        samtools view \
            -bhX -@ 1 \
            --verbosity=8 \
            --regions-file padded.bed \
            --write-index \
            -o "~{subset_prefix}.bam##idx##~{subset_prefix}.bam.bai" \
            "~{bam}" "~{bai}" &
        set +x
        sleep 1200 && export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        set -x
        wait

        while IFS=$'\t' read -r chr column2 end bed_id; do
            echo "${bed_id}"
            beg=$((column2+1))
            arr=$(samtools depth -aa -J -r "${chr}:${beg}-${end}" ~{bq_filter} ~{mq_filter} "~{subset_prefix}.bam" | awk '{print $3}' | paste -s -d, -)
            echo -e "${chr}\t${column2}\t${end}\t${bed_id}\t${arr}" >> "~{subset_prefix}.tsv"
        done < "~{bed}"

        gzip -k "~{subset_prefix}.tsv"
    >>>

    output {
        File subset_bed_depth = "~{subset_prefix}.tsv.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
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

task MergeDepth {
    input {
        Array[File] perbase_depth_split
        String out_prefix
    }

    command <<<
        set -euxo pipefail

        cat ~{sep=" " perbase_depth_split} | gunzip > merged.tsv
        wc -l merged.tsv

        sort -V merged.tsv | gzip > "~{out_prefix}.tsv.gz"
    >>>

    output {
        File merged_sorted = "~{out_prefix}.tsv.gz"
    }

    runtime {
        cpu:            2
        memory:         "8 GiB"
        disks:          "local-disk 500 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker:         "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
