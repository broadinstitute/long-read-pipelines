version 1.0

import "tasks/Finalize.wdl" as FF

workflow CoverageOverBedsPerBase {
    input {
        File bam
    	File four_col_bed
        String participant_name

        Int? min_base_q
        Int? min_map_q

    	Boolean fail_on_non_4_col_bed

        String workflow_root_name
        String gcs_out_root_dir
    }

    String dir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_root_name}/~{participant_name}/alignments"

    call Collect {input:
        bam = bam, four_col_bed = four_col_bed,
        min_base_q = min_base_q, min_map_q = min_map_q,
        fail_on_non_4_col_bed = fail_on_non_4_col_bed
    }

    call FF.FinalizeToFile { input: outdir = dir, file = Collect.depth }

    output {
        File depth = FinalizeToFile.gcs_path
    }
}

task Collect {
    input {
    	File bam
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

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

        out="~{prefix}.samtools.depth.txt"
        echo -e "chr\tBED_beg\tBED_end\tid\tdepth_array_aa.J" > "${out}"

        # Note it's done this way due to https://github.com/samtools/samtools/issues/1581
        while IFS=$'\t' read -r chr column2 end bed_id; do
            echo "${bed_id}"
            beg=$((column2+1))
            arr=$(samtools depth -aa -J -r "${chr}:${beg}-${end}" ~{bq_filter} ~{mq_filter} "~{bam}" | awk '{print $3}' | paste -s -d, -)
            echo -e "${chr}\t${column2}\t${end}\t${bed_id}\t${arr}" >> "${out}"
        done < "~{four_col_bed}"

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
        max_retries:        0,
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