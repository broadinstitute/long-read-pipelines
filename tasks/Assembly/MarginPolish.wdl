version 1.0

import "../../structs/Structs.wdl"
# Run the MarginPolish reads polisher on a given set of data to produce a better assembly.
#
# This WDL assumes that the marginPolish utility is on the path of the given docker image.
# For an image that meets this criterion, see (on dockerhub):
#   jonnsmith/lrma_alignment_toolbox:latest
#
# For more information on the MarginPolish tool, see the official GitHub repository:
#   https://github.com/UCSC-nanopore-cgl/MarginPolish
#
# Description of inputs:
#
#   Required:
#     String docker                - Docker image in which to run this task.
#     File bam_file                - Bam file containing reads aligned to the given assembly that should be polished.
#     File bam_index               - Index for input bam_file.
#     File assembly_fasta          - FASTA file containing the assembly/reference to which the given reads are aligned.
#     File config_json             - MarginPolish runtime configuration JSON file.
#
#   Optional:
#     String? interval             - Interval string over which to run this task (contig:start-end).
#     String? log_level            - Set the log level [default = info]
#     String? feature_type         - Feature type for output features of chunks for HELEN.  Valid types:
#                                      splitRleWeight:  [default] run lengths split into chunks
#                                      simpleWeight:    weighted likelihood from POA nodes (non-RLE)
#     Int? split_rle_weight_max_rl - Max run length (for 'splitRleWeight' type only) [default = 10]
#     Boolean? true_reference_bam  - True reference aligned to ASSEMBLY_FASTA, for HELEN features.  Setting this parameter will include labels in output.
#
#   Runtime:
#     Int  mem                     - Amount of memory to give to the machine running each task in this workflow.
#     Int  preemptible_attempts    - Number of times to allow each task in this workflow to be preempted.
#     Int  disk_space_gb           - Amount of storage disk space (in Gb) to give to each machine running each task in this workflow.
#     Int  cpu                     - Number of CPU cores to give to each machine running each task in this workflow.
#     Int  boot_disk_size_gb       - Amount of boot disk space (in Gb) to give to each machine running each task in this workflow.
#
workflow MarginPolish {

    input {
        String docker = "us.gcr.io/broad-dsp-lrma/lr-shasta-marginpolish-helen:0.0.1"

        File bam_file
        File bam_index

        File assembly_fasta

        File config_json

        String? interval
        String? log_level

        String? feature_type
        Int? split_rle_weight_max_rl
        Boolean? true_reference_bam

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    call MarginPolishTask {
        input:
            docker_image              = docker,

            bam_file                  = bam_file,
            bam_index                 = bam_index,

            assembly_fasta            = assembly_fasta,

            config_json               = config_json,

            interval                  = interval,
            log_level                 = log_level,
            feature_type              = feature_type,
            split_rle_weight_max_rl   = split_rle_weight_max_rl,
            true_reference_bam        = true_reference_bam,

            mem_gb                    = mem_gb,
            preemptible_attempts      = preemptible_attempts,
            disk_space_gb             = disk_space_gb,
            cpu                       = cpu,
            boot_disk_size_gb         = boot_disk_size_gb
    }

    output {
        File polished_assembly        = MarginPolishTask.polished_assembly
        File images_for_helen         = MarginPolishTask.images_for_helen
        File poa_tar                  = MarginPolishTask.poa_tar
        File repeat_count_tar         = MarginPolishTask.repeat_count_tar
        File timing_info              = MarginPolishTask.timing_info
    }
}

task MarginPolishTask {

    # ------------------------------------------------
    # Input args:
    input {
        # Docker:
        String docker_image

        # Required:
        File bam_file
        File bam_index

        File assembly_fasta

        File config_json

        # Optional:
        String? interval
        String? log_level

        String? feature_type
        Int? split_rle_weight_max_rl
        Boolean? true_reference_bam

        # Runtime Options:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    # ------------------------------------------------
    # Process input args:

    String bam_extension = ".bam"
    String bam_file_base_name = basename( bam_file, bam_extension )

    String true_reference_bam_arg = if defined(true_reference_bam) && true_reference_bam then "-u" else ""

    String output_base_name = bam_file_base_name + '.margin_polish'

    String poa_dir_name = 'poa'
    String poa_output_base_name = poa_dir_name + '/' + output_base_name

    String repeat_count_dir_name = 'repeat_counts'
    String repeat_counts_output_base_name = repeat_count_dir_name + '/' + output_base_name

    String helen_images_tar_name = output_base_name + '.helen_images.tar'

    String timing_output_file = bam_file_base_name + ".timingInformation.txt"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements:
    Int default_ram_mb = 3 * 1024

    Float reads_size_gb = size(bam_file, "GiB") + size(bam_index, "GiB")
    Int default_disk_space_gb = ceil((reads_size_gb * 2) + 20)

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb
    Int command_mem = machine_mem - 1024

    # ------------------------------------------------
    # Misc. helpers:
    String dollar = "$"

    # ------------------------------------------------
    # Run our command:
    command {

        # Find out how many threads to use:
        num_threads=~{cpu}
        found_count=false
        which lscpu &> /dev/null
        r=$?
        if [[ $r -eq 0 ]] ; then
            num_threads=$( lscpu | grep ^CPU\(s\) | sed 's#^CPU(s):[\t ]*##' )
            found_count=true
        fi

        if ! $found_count ; then
            if [[ -f /proc/cpuinfo ]] ; then
                num_threads=$( cat /proc/cpuinfo | grep ^processor | tail -n1 | sed 's#processor[ \t]*:[\t ]*##' )
            fi
        fi

        ############################

        # Now we can run marginPolish:
        set -e
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        # Make directories for the POA and RepeatCounts:
        mkdir ~{poa_dir_name} ~{repeat_count_dir_name}

        marginPolish \
            ~{bam_file} \
            ~{assembly_fasta} \
            ~{config_json} \
            -t $num_threads \
            ~{"-r " + interval} \
            -f \
            -i ~{repeat_counts_output_base_name} \
            -j ~{poa_output_base_name} \
            -o ~{output_base_name} \
            ~{"-a" + log_level} \
            ~{"-F" + feature_type} \
            ~{"-L" + split_rle_weight_max_rl} \
            ~{true_reference_bam_arg}

        String? feature_type
        Int? split_rle_weight_max_rl
        Boolean? true_reference_bam

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ~{timing_output_file}
        elapsedTime=`python -c "print( $endTime - $startTime )"`
        echo "Elapsed Time: $elapsedTime" >> ~{timing_output_file}

        # Tar the images for HELEN, the POA dir, and the Repeat Count dir
        # (in parallel and without compressing them to save time):
        tar -cf ~{helen_images_tar_name} *.T[0-9]*.h5 &
        tar -cf ~{poa_dir_name}.tar ~{poa_dir_name} &
        tar -cf ~{repeat_count_dir_name}.tar ~{repeat_count_dir_name} &
        wait
    }

    # ------------------------------------------------
    # Runtime settings:
     runtime {
         docker: docker_image
         memory: machine_mem + " MB"
         disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
         bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
         preemptible: select_first([preemptible_attempts, 0])
         cpu: select_first([cpu, 1])
     }

    # ------------------------------------------------
    # Outputs:
    output {
        File polished_assembly      = "${output_base_name}.fa"
        File images_for_helen       = "${helen_images_tar_name}"
        File poa_tar                = "${poa_dir_name}.tar"
        File repeat_count_tar       = "${repeat_count_dir_name}.tar"
        File timing_info            = "${timing_output_file}"
    }
 }

