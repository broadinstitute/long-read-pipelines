version 1.0

workflow MethBat{
    meta{
        description: "a workflow that using MethBat to extract methylation patterns"
    }

    input{
        File combined_bed
        File? combined_bed_index
        File hap1_bed
        File? hap1_bed_index
        File hap2_bed
        File? hap2_bed_index
        String prefix
        String sex
        File CpG_region_file
        String outputprefix
        Boolean call_segmentation
    }

    if (call_segmentation) {
        call MethBat_segmentation{input: 
            combined_bed=combined_bed, 
            combined_bed_index=combined_bed_index, 
            hap1_bed = hap1_bed,
            hap1_bed_index = hap1_bed_index,
            hap2_bed = hap2_bed,
            hap2_bed_index = hap2_bed_index,
            prefix = prefix, 
            outputprefix = outputprefix
        }
    }

    call MethBat_profile {input: 
        combined_bed = combined_bed, 
        combined_bed_index=combined_bed_index, 
        hap1_bed = hap1_bed,
        hap1_bed_index = hap1_bed_index,
        hap2_bed = hap2_bed,
        hap2_bed_index = hap2_bed_index,
        prefix = prefix ,
        sample_id = prefix,
        sex = sex,
        CpG_region = CpG_region_file
    }

    output{
        File? region_bed = MethBat_segmentation.region_bed
        File? asm_bedgraph = MethBat_segmentation.asm_bedgraph
        File? combined_methyl_bedgraph = MethBat_segmentation.combined_methyl_bedgraph
        File? hap1_methyl_bedgraph = MethBat_segmentation.hap1_methyl_bedgraph
        File? hap2_methyl_bedgraph = MethBat_segmentation.hap2_methyl_bedgraph
        File profile = MethBat_profile.profile
        String stat_methbat_methylated_count = MethBat_profile.stat_methbat_methylated_count
        String stat_methbat_unmethylated_count = MethBat_profile.stat_methbat_unmethylated_count
        String stat_methbat_asm_count = MethBat_profile.stat_methbat_asm_count
    }
}

task MethBat_segmentation{
    input{
        File combined_bed
        File? combined_bed_index
        File hap1_bed
        File? hap1_bed_index
        File hap2_bed
        File? hap2_bed_index
        String prefix
        String outputprefix
        RuntimeAttr? runtime_attr_override
    }

    String root_directory = "/mnt/disks/cromwell_root"

    command <<<
        set -euxo pipefail

        mkdir -p input
        mv ~{combined_bed} input
        if [[ "~{defined(combined_bed_index)}" == "true" ]]; then
            mv "~{combined_bed_index}" input
        fi

        mv ~{hap1_bed} input
        if [[ "~{defined(hap1_bed_index)}" == "true" ]]; then
            mv "~{hap1_bed_index}" input
        fi
        mv ~{hap2_bed} input
        if [[ "~{defined(hap2_bed_index)}" == "true" ]]; then
            mv "~{hap2_bed_index}" input
        fi

        ls input

        methbat segment --input-prefix input/~{prefix} \
                        --output-prefix ~{outputprefix} \
                        --enable-haplotype-segmentation
    >>>

    output{
        File region_bed = "~{outputprefix}.meth_regions.bed"
        File asm_bedgraph = "~{outputprefix}.asm.bedgraph"
        File combined_methyl_bedgraph = "~{outputprefix}.combined_methyl.bedgraph"
        File hap1_methyl_bedgraph = "~{outputprefix}.hap1.bedgraph"
        File hap2_methyl_bedgraph = "~{outputprefix}.hap2.bedgraph"
        
    }


    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/methylationanalysis:v1"
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

task MethBat_profile{
    input{
        File combined_bed
        File? combined_bed_index
        File hap1_bed
        File? hap1_bed_index
        File hap2_bed
        File? hap2_bed_index
        File CpG_region
        String prefix
        String sex
        String sample_id
        RuntimeAttr? runtime_attr_override
    }

    String root_directory = "/mnt/disks/cromwell_root"

    command <<<
        mkdir -p input
        mkdir -p output
        mv ~{combined_bed} input
        if [[ "~{defined(combined_bed_index)}" == "true" ]]; then
            mv "~{combined_bed_index}" input
        fi
        mv ~{hap1_bed} input
        if [[ "~{defined(hap1_bed_index)}" == "true" ]]; then
            mv "~{hap1_bed_index}" input
        fi
        mv ~{hap2_bed} input
        if [[ "~{defined(hap2_bed_index)}" == "true" ]]; then
            mv "~{hap2_bed_index}" input
        fi

        ls input

        methbat profile \
            --input-prefix input/~{prefix} \
            --input-regions ~{CpG_region} \
            --output-region-profile ./output/~{sample_id}.meth_region_stats.tsv \
            --profile-label ALL \
            --profile-label ~{sex} 

        # count for three most interesting methylation summary_labels
        awk '$5=="Methylated" {print}' "./output/~{sample_id}.meth_region_stats.tsv" \
        | wc -l > methylated_count.txt
        awk '$5=="Unmethylated" {print}' "./output/~{sample_id}.meth_region_stats.tsv" \
        | wc -l > unmethylated_count.txt
        awk '$5=="AlleleSpecificMethylation" {print}' "./output/~{sample_id}.meth_region_stats.tsv" \
        | wc -l > asm_count.txt
    >>> 

    output{
        File profile = "output/~{sample_id}.meth_region_stats.tsv"
        String stat_methbat_methylated_count = read_string("methylated_count.txt")
        String stat_methbat_unmethylated_count = read_string("unmethylated_count.txt")
        String stat_methbat_asm_count = read_string("asm_count.txt")
    }


    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/methylationanalysis:v1"
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


struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

struct DataTypeParameters {
    Int num_shards
    String map_preset
}

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
        boot_disk_gb:       10,
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