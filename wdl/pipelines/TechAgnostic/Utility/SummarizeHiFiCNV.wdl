version 1.0

import "../../../structs/Structs.wdl"

import "../../../tasks/Utility/Finalize.wdl" as FF

workflow SummarizeHiFiCNV {
    meta {
        description:
        "Summarizes the outputs of HiFiCNV across multiple samples into a single bed file of CNV calls."
    }
    input {
        Array[Map[String, File]] samples_hificnv_results
        File ref_map_file
        String gcs_out_root_dir
        String cohort_name
    }

    output {
        File merged_bed = FinalizeToFile.gcs_path
    }

    scatter (mm in samples_hificnv_results) {
        call ConvertToBed { input: input_vcf = mm['vcf']}
    }

    Map[String, String] ref_map = read_map(ref_map_file)
    call MergeAndSortBed { input: bed_files = ConvertToBed.bed,  fai_idx = ref_map['fai'] }

    String workflow_name = 'HiFiCNV'
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_name}"
    call FF.FinalizeToFile { input: outdir = outdir, file = MergeAndSortBed.bed, name = "~{cohort_name}.HiFiCNV.summary.bed" }
}

task ConvertToBed {
    input {
        File input_vcf
        RuntimeAttr? runtime_attr_override
    }

    command <<<
    set -eux

        smid=$(zgrep "^#CHROM" "~{input_vcf}" | awk -F '\t' '{print $10}')
        bcftools view -f "PASS,." "~{input_vcf}" \
            | awk -F '\t' 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]{1,2}|X|Y|M)$/ {print}' \
            | awk -v sm="$smid" -F '\t' 'BEGIN{OFS="\t"}
                {
                    match($8, /SVTYPE=[A-Z]+/); {
                        tp=substr($8, RSTART+7, RLENGTH-7);
                    }
                    match($8, /END=[0-9]+/); {
                        end=substr($8, RSTART+4, RLENGTH-4);
                    }
                    {split($10,a,":"); gt=a[1]; cn=a[2]}
                print $1, $2, end, tp";"gt";"sm, cn
                }
                ' \
            > converted.bed

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
        preemptible_tries:  1,
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
    output {
        File bed = "merged.and.sorted.bed"
    }

    command <<<
    set -eux

        cat ~{sep=" " bed_files} \
        > merged.bed

        bedtools sort -faidx ~{fai_idx} -i merged.bed \
        > merged.and.sorted.bed
    >>>


    #########################
    Int disk_size = 10 + ceil(size(bed_files, "GiB"))

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
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
