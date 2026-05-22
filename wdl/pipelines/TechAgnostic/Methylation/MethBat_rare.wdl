version 1.0

workflow MethBat_rare{
    meta{
        description: "a workflow that using MethBat to extract methylation patterns"
    }

    input{

        File background_manifest
        File case_combined_bed
        File? case_combined_bed_index
        File case_hap1_bed
        File? case_hap1_bed_index
        File case_hap2_bed
        File? case_hap2_bed_index
        String case_sample_id
        String case_sex
        File CpG_region_file
        String outputprefix
    }

    call Generate_cohort_file{input:
        manifest_file = background_manifest,
        prefix = outputprefix
    }

    call MethBat_profile {input: 
        combined_bed = case_combined_bed, 
        combined_bed_index = case_combined_bed_index, 
        hap1_bed = case_hap1_bed,
        hap1_bed_index = case_hap1_bed_index,
        hap2_bed = case_hap2_bed,
        hap2_bed_index = case_hap2_bed_index,
        prefix = outputprefix ,
        sample_id = case_sample_id,
        sex = case_sex,
        CpG_region = Generate_cohort_file.cohort_file
    }


    output{
        File output_cohort = Generate_cohort_file.cohort_file
        File rare_profile = MethBat_profile.profile
    }
}


task ProcessManifest{
    input{
        File manifest_file
        RuntimeAttr? runtime_attr_override
    }

    String root_directory = "/mnt/disks/cromwell_root"

    command <<<
        cat ~{manifest_file} | cut -f 1 > sampleid.txt
        cat ~{manifest_file} | cut -f 2 > combined_bed.txt
        cat ~{manifest_file} | cut -f 3 > hap1_bed.txt
        cat ~{manifest_file} | cut -f 4 > hap2_bed.txt
        cat ~{manifest_file} | cut -f 5 > sex.txt

    >>> 

    output{
        Array[File] background_combined_bed = read_lines("combined_bed.txt")
        Array[File] background_hap1_bed = read_lines("hap1_bed.txt")

        Array[File] background_hap2_bed = read_lines("hap2_bed.txt")

        Array[String] background_sample_ids = read_lines("sampleid.txt")
        Array[String] background_sex_list = read_lines("sex.txt")

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
        mv ~{combined_bed_index} input
        mv ~{hap1_bed} input
        mv ~{hap1_bed_index} input
        mv ~{hap2_bed} input
        mv ~{hap2_bed_index} input

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
        cpu_cores:          1,
        mem_gb:             4,
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

task Generate_cohort_file{
    input{
        File manifest_file
        String prefix
        RuntimeAttr? runtime_attr_override
    }

    String root_directory = "/mnt/disks/cromwell_root"

    command <<<
        set -euxo pipefail

        mkdir -p output
        cohort_tsv="output/~{prefix}.cohort.tsv"

        cat ~{manifest_file} | cut -f 1 > sampleid.txt
        cat ~{manifest_file} | cut -f 2 > profilelist.txt
        cat ~{manifest_file} | cut -f 3 > sex.txt

        profile_paths_file="profilelist.txt"
        sample_ids_file="sampleid.txt"
        sex_labels_file="sex.txt"

        # copy profiles listed in manifest to local folder
        mkdir -p input_profiles
        while IFS= read -r profile_path; do
            [[ -z "${profile_path}" ]] && continue
            if [[ "${profile_path}" == gs://* ]]; then
                gsutil cp "${profile_path}" input_profiles/
            else
                cp "${profile_path}" input_profiles/
            fi
        done < "${profile_paths_file}"


        n_profiles=$(wc -l < "${profile_paths_file}")
        n_samples=$(wc -l < "${sample_ids_file}")
        n_sex=$(wc -l < "${sex_labels_file}")
        if [[ "${n_profiles}" -ne "${n_samples}" || "${n_profiles}" -ne "${n_sex}" ]]; then
            echo "profile_list, sample_list, and sex must have the same length" >&2
            exit 1
        fi

        printf "identifier\tfilename\tlabels\n" > "${cohort_tsv}"

        paste "${sample_ids_file}" "${profile_paths_file}" "${sex_labels_file}" | \
        while IFS=$'\t' read -r sample_entry profile_path sex_label; do
            identifier=$(basename "${sample_entry}")
            local_profile_path="input_profiles/$(basename "${profile_path}")"
            printf "%s\t%s\t%s\n" "${identifier}" "${local_profile_path}" "${sex_label};UNAFFECTED" >> "${cohort_tsv}"
        done


        methbat build \
            --input-collection output/~{prefix}.cohort.tsv \
            --output-profile output/~{prefix}.output.cohort.tsv
    >>> 

    output{
        File cohort_file = "output/~{prefix}.output.cohort.tsv"
    }


    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
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

task MethBat_compare{
    input{
        File cohort_profile
        String prefix
        String baseline
        String compare
        RuntimeAttr? runtime_attr_override
    }

    String root_directory = "/mnt/disks/cromwell_root"

    command <<<
        methbat compare \
            --input-profile ~{cohort_profile} \
            --output-comparison ~{prefix}.output.tsv \
            --baseline-category ~{baseline} \
            --compare-category ~{compare}
    >>> 

    output{
        File profile = "~{prefix}.output.tsv"
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