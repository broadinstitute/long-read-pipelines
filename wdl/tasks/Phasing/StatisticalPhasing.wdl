version 1.0

import "../../structs/Structs.wdl"

task Shapeit4 {
    input{
        File vcf_input
        File vcf_index
        File mappingfile
        String region
        Int num_threads
        Int memory

        RuntimeAttr? runtime_attr_override
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    command <<<
        # add AN AC tag

        # export MONITOR_MOUNT_POINT="/cromwell_root/"
        # bash /opt/vm_local_monitoring_script.sh &> resources.log &
        # job_id=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')

        shapeit4 --input ~{vcf_input} \
                --map ~{mappingfile} \
                --region ~{region} \
                --use-PS 0.0001 \
                --sequencing \
                --output ~{region}_scaffold.bcf \
                --thread ~{num_threads} \
                --log phased.log
        
        # if ps -p "${job_id}" > /dev/null; then kill "${job_id}"; fi
    >>>

    output{
        # File resouce_monitor_log = "resources.log"
        File scaffold_vcf = "~{region}_scaffold.bcf"
    }

    #Int disk_size = 100 + ceil(2 * size(vcf_input, "GiB"))

 #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_threads,
        mem_gb:             memory,
        disk_gb:            500,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "hangsuunc/hiphase:1.3.0"
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


task Shapeit4_phaseSVs {
    input{
        File vcf_input
        File vcf_index
        File scaffold_vcf
        File mappingfile
        String region
        Int num_threads
        Int memory

        RuntimeAttr? runtime_attr_override
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    command <<<
        # add AN AC tag
        bcftools index ~{scaffold_vcf}

        shapeit4 \
        --input ~{vcf_input} \
        --scaffold ~{scaffold_vcf} \
        --map ~{mappingfile} \
        --region ~{region} \
        --use-PS 0.0001 \
        --sequencing \
        --output ~{region}_finalsv_scaffold.bcf \
        --thread ~{num_threads} \
        --log phased.log
    
    >>>

    output{
        File final_phased_vcf = "~{region}_finalsv_scaffold.bcf"
    }

#########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_threads,
        mem_gb:             memory,
        disk_gb:            500,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "hangsuunc/hiphase:1.3.0"
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


task MeasurePhasingSwitchErrorRate {
    input{
        File truth_bcf
        File truth_bcf_index
        File test_bcf
        File test_bcf_index
        String region
        String outputprefix
        Int num_threads
    }
    command <<<
        switch_static --validation ~{truth_bcf} --estimation ~{test_bcf} --region ~{region} --output ~{outputprefix} --thread ~{num_threads}
    >>>

    output{
        Array[File] output_files = glob("*")
    }

    Int disk_size = 100 + ceil(2 * (size(truth_bcf, "GiB") + size(test_bcf, "GiB")))

    runtime {
        cpu: 8
        memory: "32 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "lindonkambule/shapeit5_2023-05-05_d6ce1e2:v5.1.1"
    }
}
