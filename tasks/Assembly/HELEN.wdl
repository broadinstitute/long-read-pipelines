version 1.0

import "../../structs/Structs.wdl"
# Run the HELEN reads polisher on a given set of data to produce a better assembly.
#
# This WDL assumes that the HELEN utilities are in the following location on
# For an image that meets this criterion, see (on dockerhub):
#   jonnsmith/lrma_alignment_toolbox:latest
#
# For more information on the HELEN tool, see the official GitHub repository:
#   https://github.com/kishwarshafin/helen
#
# Valid values for GPUs are (as of 20191016):
#    nvidia-tesla-p100
#    nvidia-tesla-k80
#
# Note: This assumes that the reads have already been processed by the `MarginPolish` tool.
#
# Description of inputs:
#
#   Required:
#     String docker                        - Docker image in which to run this task.
#     File helen_images_from_margin_polish - Tar file containing images created by marginPolish for use with HELEN.
#     File model_pickle                    - Python Pickle file contianing the model to use for polishing with HELEN.
#
#   Optional:
#     Int? batch_size                      - Batch size for testing.  Set to 512 or 1024 for a balanced execution time. (=512)
#
#   Runtime:
#     Int? num_gpus                        - How many GPUs to use for accelerated processing.
#     String? gpu_type                     - Which kind of GPU to use for accelerated processing.
#     Int  mem                             - Amount of memory to give to the machine running each task in this workflow.
#     Int  preemptible_attempts            - Number of times to allow each task in this workflow to be preempted.
#     Int  disk_space_gb                   - Amount of storage disk space (in Gb) to give to each machine running each task in this workflow.
#     Int  cpu                             - Number of CPU cores to give to each machine running each task in this workflow.
#     Int  boot_disk_size_gb               - Amount of boot disk space (in Gb) to give to each machine running each task in this workflow.
#
workflow HELEN {

    input {
        String docker = "us.gcr.io/broad-dsp-lrma/lr-shasta-marginpolish-helen:0.0.1"

        File helen_images_from_margin_polish
        File model_pickle

        Int? batch_size

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
        Int? num_gpus
        String? gpu_type
    }

    call CallConsensusTask {
        input:
            docker_image                    = docker,

            helen_images_from_margin_polish = helen_images_from_margin_polish,
            model_pickle                    = model_pickle,

            batch_size                      = batch_size,

            mem_gb                          = mem_gb,
            preemptible_attempts            = preemptible_attempts,
            disk_space_gb                   = disk_space_gb,
            cpu                             = cpu,
            boot_disk_size_gb               = boot_disk_size_gb,
            num_gpus                        = num_gpus,
            gpu_type                        = gpu_type
    }

    call StitchTask {
        input:
            docker_image              = docker,

            consensus_hdf5_file       = CallConsensusTask.consensus_hdf5_file,

            mem_gb                    = mem_gb,
            preemptible_attempts      = preemptible_attempts,
            disk_space_gb             = disk_space_gb,
            cpu                       = cpu,
            boot_disk_size_gb         = boot_disk_size_gb
    }

    output {
        File timing_info              = StitchTask.timing_info
    }
}

task CallConsensusTask {

    # ------------------------------------------------
    # Input args:
    input {
        # Docker:
        String docker_image

        # Required:
        File helen_images_from_margin_polish
        File model_pickle

        # Optional:
        Int? batch_size

        # Runtime Options:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
        Int? num_gpus
        String? gpu_type
    }

    # ------------------------------------------------
    # Process input args:

    String input_file_base_name = basename( basename( helen_images_from_margin_polish, '.gz' ), '.tar' )

    String output_base_name = input_file_base_name + '.HELEN_prediction'
    String output_file_dir = "helen_output"

    String timing_output_file = input_file_base_name + ".timingInformation.txt"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements:
    Int default_ram_mb = 3 * 1024

    Float reads_size_gb = size(helen_images_from_margin_polish, "GiB") + size(model_pickle, "GiB")
    Int default_disk_space_gb = ceil((reads_size_gb * 2) + 20)

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb
    Int command_mem = machine_mem - 1024

    # ------------------------------------------------
    # Misc. helpers:
    String dollar = "$"

    # ------------------------------------------------
    # Task-Specific runtime settings:
    String helen_install_path = '/helen'
    String batch_size_arg = select_first([batch_size, 512])
    String gpu_arg = if (defined(num_gpus) && num_gpus > 0) then "-g" else ""

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

        # Untar the inputs:

        # Check if we need to make a directory:
        image_dir="input_image_dir"
        mkdir $image_dir
        cd $image_dir
        tar -xf ~{helen_images_from_margin_polish}
        if [[ $( tar -tf ~{helen_images_from_margin_polish} | head -n 1 | grep -c '/' ) -ne 0 ]] ; then
            mv $( find . -type f -name \*T[0-9]*.h5 |xargs echo ) .
        fi
        cd -

        # Now we can run the HELEN consensus caller:
        set -e
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        python3 ~{helen_install_path}/call_consensus.py \
            -i $image_dir \
            -o ~{output_file_dir} \
            -b ~{batch_size_arg} \
            -m ~{model_pickle} \
            -p ~{output_base_name} \
            -w 0 \
            -t $num_threads \
            ~{gpu_arg}


        # NOTE: The 2 BLANK LINES ABOVE THIS ARE REQUIRED!!!!!!!!!!
        #       This is because the gpu_arg could actually be empty.

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ~{timing_output_file}
        elapsedTime=`python -c "print( $endTime - $startTime )"`
        echo "Elapsed Time: $elapsedTime" >> ~{timing_output_file}
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
        gpuType: select_first([gpu_type, ''])
        gpuCount: select_first([num_gpus, 0])
     }

    # ------------------------------------------------
    # Outputs:
    output {
      File consensus_hdf5_file    = "${output_file_dir}/${output_base_name}.hdf"
      File timing_info            = "${timing_output_file}"
    }
 }

task StitchTask {

    # ------------------------------------------------
    # Input args:
    input {
        # Docker:
        String docker_image

        # Requied input files:
        File consensus_hdf5_file

        # Runtime Options:
        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    # ------------------------------------------------
    # Process input args:

    String input_file_base_name = basename( consensus_hdf5_file, '.hdf' )

    String output_base_name = input_file_base_name + 'stitched'
    String output_dir_name = output_base_name

    String timing_output_file = input_file_base_name + ".timingInformation.txt"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements:
    Int default_ram_mb = 3 * 1024

    Float reads_size_gb = size(consensus_hdf5_file, "GiB")
    Int default_disk_space_gb = ceil((reads_size_gb * 2) + 20)

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb
    Int command_mem = machine_mem - 1024

    # ------------------------------------------------
    # Misc. helpers:
    String dollar = "$"

    String helen_install_path = '/helen'

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

        # Make a place for output:
        mkdir ~{output_dir_name}

        # Now we can run marginPolish:
        set -e
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        python3 ~{helen_install_path}/stitch.py \
            -i ~{consensus_hdf5_file} \
            -t $num_threads \
            -o ~{output_dir_name} \
            -p ~{output_base_name}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ~{timing_output_file}
        elapsedTime=`python -c "print( $endTime - $startTime )"`
        echo "Elapsed Time: $elapsedTime" >> ~{timing_output_file}
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
      File polished_consensus     = "${output_dir_name}/${output_base_name}.fa"
      File timing_info            = "${timing_output_file}"
    }
 }
