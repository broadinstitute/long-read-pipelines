version 1.0

task Shapeit4 {
    input{
        File vcf_input
        File vcf_index
        File mappingfile
        String region
        Int num_threads
    }
    command <<<
        # add AN AC tag

        export MONITOR_MOUNT_POINT="/cromwell_root/"
        bash /opt/vm_local_monitoring_script.sh &> resources.log &
        job_id=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')

        shapeit4 --input ~{vcf_input} --map ~{mappingfile} --region ~{region} --use-PS 0.0001 --sequencing --output ~{region}_scaffold.bcf --thread ~{num_threads} --log phased.log
        
        if ps -p "${job_id}" > /dev/null; then kill "${job_id}"; fi
    >>>

    output{
        File resouce_monitor_log = "resources.log"
        File scaffold_vcf = "~{region}_scaffold.bcf"
    }

    #Int disk_size = 100 + ceil(2 * size(vcf_input, "GiB"))

    runtime {
        cpu: 40
        memory: "500 GiB"
        disks: "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "hangsuunc/hiphase:0.7.2"
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

    #Int disk_size = 100 + ceil(2 * size(vcf_input, "GiB"))

    runtime {
        cpu: 40
        memory: "240 GiB"
        disks: "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "hangsuunc/hiphase:0.7.2"
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
