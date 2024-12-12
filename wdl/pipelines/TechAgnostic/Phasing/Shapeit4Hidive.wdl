version 1.0

workflow Shapeit4Hidive {   

    input {
        File joint_vcf
        File joint_vcf_tbi
        File genetic_mapping_tsv_for_shapeit4
        File regionlist
        File chromosomelist
        Int shapeit4_num_threads
        Int shapeit4_memory
        String shapeit4_extra_args
        String output_prefix

    }

    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit4)

    Array[String] region_list = read_lines(regionlist)

    Array[String] chromosome_list = read_lines(chromosomelist)

    scatter (j in range(length(region_list))) {
        String chromosome = chromosome_list[j]
        call Shapeit4 { input:
            vcf_input = joint_vcf,
            vcf_index = joint_vcf_tbi,
            mappingfile = genetic_mapping_dict[chromosome],
            region = region_list[j],
            prefix = output_prefix + "." + chromosome + ".shard-" + j + ".phased",
            num_threads = shapeit4_num_threads,
            memory = shapeit4_memory,
            extra_args = shapeit4_extra_args
        }
    }

    call LigateVcfs { input:
        vcfs = Shapeit4.phased_bcf,
        prefix = output_prefix + "." + ".phased.ligated"
    }

    output {

        File phased_vcf_gz = LigateVcfs.ligated_vcf_gz
        File phased_vcf_gz_tbi = LigateVcfs.ligated_vcf_gz_tbi

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


task LigateVcfs {

    input {
        Array[File] vcfs
        Array[File]? vcf_idxs
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(vcfs, "GB")) + 1

    command <<<
        set -euxo pipefail
        if ! ~{defined(vcf_idxs)}; then
            for ff in ~{sep=' ' vcfs}; do bcftools index $ff; done
        fi

        wget https://github.com/odelaneau/shapeit5/releases/download/v5.1.1/ligate_static
        chmod +x ligate_static

        ./ligate_static --input ~{write_lines(vcfs)} --output ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File ligated_vcf_gz = "~{prefix}.vcf.gz"
        File ligated_vcf_gz_tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Shapeit4 {
    input{
        File vcf_input
        File vcf_index
        File mappingfile
        String region
        String prefix
        Int num_threads
        Int memory
        String extra_args

        RuntimeAttr? runtime_attr_override
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    command <<<
        # add AN AC tag

        # export MONITOR_MOUNT_POINT="/cromwell_root/"
        # bash /opt/vm_local_monitoring_script.sh &> resources.log &
        # job_id=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')

        shapeit4.2 --input ~{vcf_input} \
                --map ~{mappingfile} \
                --region ~{region} \
                --sequencing \
                --output ~{prefix}.bcf \
                --thread ~{num_threads} \
                ~{extra_args}

        # if ps -p "${job_id}" > /dev/null; then kill "${job_id}"; fi
    >>>

    output{
        # File resouce_monitor_log = "resources.log"
        File phased_bcf = "~{prefix}.bcf"
    }

    #Int disk_size = 100 + ceil(2 * size(vcf_input, "GiB"))

 #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_threads,
        mem_gb:             memory,
        disk_gb:            100,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "hangsuunc/shapeit4:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}