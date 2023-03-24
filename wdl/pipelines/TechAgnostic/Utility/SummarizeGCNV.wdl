version 1.0

import "../../../structs/Structs.wdl"

workflow SummarizeGCNV {
    input {
        Array[File] gcnv_cohort_geontyped_segments_vcfs
        File ref_map_file
    }

    scatter (vcf in gcnv_cohort_geontyped_segments_vcfs) {
        call ConvertToBed { input: input_vcf = vcf}
    }

    Map[String, String] ref_map = read_map(ref_map_file)
    call MergeAndSortBed { input: bed_files = ConvertToBed.bed,  fai_idx = ref_map['fai'] }

    output {
        File merged_bed = MergeAndSortBed.bed
    }
}

task ConvertToBed {
    input {
        File input_vcf
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eux

        smid=$(zgrep "^#CHROM" "~{input_vcf}" | awk -F '\t' '{print $10}')
        zgrep -v "^#" "~{input_vcf}" \
            | awk -F '\t' 'BEGIN{OFS="\t"} {if($7=="PASS" || $7==".") print}' \
            | awk -F '\t' 'BEGIN{OFS="\t"} $1 ~ /^chr[0-9]{1,2}$/ {print}' \
            | awk -v sm="$smid" -F '\t' 'BEGIN{OFS="\t"}
                {
                    match($8, /END=[0-9]+/); {
                        end=substr($8, RSTART+4, RLENGTH-4);
                    }
                    {split($10,a,":"); gt=a[1]; cn=a[2]}
                    {if(cn>2) tp="DUP"; else if(cn<2) tp="DEL"; else tp="."}
                {if(cn==2) next; print $1, $2, end, tp";"gt";"sm, cn}
                }
                ' \
            > converted.bed
        touch converted.bed
        cat converted.bed
    >>>

    output {
        File bed = "converted.bed"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
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

task MergeAndSortBed {
    input {
        Array[File] bed_files
        File fai_idx
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + ceil(size(bed_files, "GiB"))

    command <<<
        set -eux

        cat ~{sep=" " bed_files} > merged.bed
        bedtools sort -faidx ~{fai_idx} -i merged.bed > merged.and.sorted.bed
    >>>

    output {
        File bed = "merged.and.sorted.bed"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
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
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
