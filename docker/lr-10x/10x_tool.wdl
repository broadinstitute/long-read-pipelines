workflow ten_x_tool {
    String docker_image
    File bam_file
    File bam_index
    File adapter_sequence
    File reverse_adapter_sequence
    File whitelist_10x
    File whitelist_illumina
    File analysis_script
    Int? read_end_length = 80

    Int? mem_gb
    Int? preemptible_attempts
    Int? disk_space_gb
    Int? cpu
    Int? boot_disk_size_gb

    call ten_x_tool_task {
        input:
            docker_image = docker_image,
            bam_file = bam_file,
            bam_index = bam_index,
            adapter_sequence = adapter_sequence,
            reverse_adapter_sequence = reverse_adapter_sequence,
            whitelist_10x = whitelist_10x,
            whitelist_illumina = whitelist_illumina,
            analysis_script = analysis_script,
            read_end_length = read_end_length,

            mem_gb = mem_gb,
            preemptible_attempts = preemptible_attempts,
            disk_space_gb = disk_space_gb,
            cpu = cpu,
            boot_disk_size_gb = boot_disk_size_gb
    }
    output {
      File rg_stats          = ten_x_tool_task.rg_stats
      File barcode_stats          = ten_x_tool_task.barcode_stats
      File starcode          = ten_x_tool_task.starcode
      File output_bam          = ten_x_tool_task.output_bam
    }
}

task ten_x_tool_task {
    String docker_image
    File bam_file
    File bam_index
    File adapter_sequence
    File reverse_adapter_sequence
    File whitelist_10x
    File whitelist_illumina
    File analysis_script
    Int read_end_length

    Int? boot_disk_size_gb
    Int? cpu
    Int? disk_space_gb
    Int? mem_gb
    Int? preemptible_attempts


    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 3 * 1024

    Float reads_size_gb = size(bam_file, "GiB") + size(bam_index, "GiB")
    Int default_disk_space_gb = ceil((reads_size_gb * 2) + 20)

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb
    Int command_mem = machine_mem - 1024

    String timing_output_file = "timingInformation.txt"

    String output_name = basename(bam_file)

    command {
        set -e

        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ${timing_output_file}

        source activate 10x_tool

        bwa index ${adapter_sequence}
        bwa index ${reverse_adapter_sequence}

        python ${analysis_script} \
            --bam=${bam_file} \
            --adapter=${adapter_sequence} \
            --reverse-adapter=${reverse_adapter_sequence} \
            --whitelist-10x=${whitelist_10x} \
            --whitelist-illumina=${whitelist_illumina} \
            --name=${output_name} \
            --read-end-length=${read_end_length}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ${timing_output_file}
        elapsedTime=`python -c "print( $endTime - $startTime )"`
        echo "Elapsed Time: $elapsedTime" >> ${timing_output_file}

    }
    runtime {
        docker: docker_image
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: 0
        cpu: select_first([cpu, 1])
    }
    output {
      File rg_stats          = "${output_name}_rg_stats.tsv"
      File barcode_stats          = "${output_name}_barcode_stats.tsv"
      File starcode          = "${output_name}_starcode.tsv"
      File output_bam          = "${output_name}.bam"
    }
}