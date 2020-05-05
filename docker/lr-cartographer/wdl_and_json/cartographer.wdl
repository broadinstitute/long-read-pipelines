version 1.0
# Run Cartographer on a dataset to map out known and unknown segments of reads.
#
# Description of inputs:
#
#   Required:
#     File reads_fasta           - SAM/BAM/FASTA/FASTQ file containing reads for which to determine the layout.
#     File segments_fasta        - FASTA file containing the sequences and names of known possible segments in the reads.
#     File layout_order          - Text file containing the expected order of the known segments in each read.
#
#   Optional:
#
#     Float min_qual             - Minimum quality for good alignment.
#     Int min_bases              - Minimum number of bases for an alignment to be retained.
#     Float prec_known           - Probability of recombination for known segment alignment.
#     Float prec_unknown         - Probability of recombination for UNKNOWN segment alignment.
#     Boolean use_mosaic_aligner - If True, will use MosaicAligner instead of Tesserae for alignments.
#     Int barcode_length         - Length of the barcode sequence.
#     File barcode_file          - File containing a list of all possible barcodes.  If specified, must also specify --barcode_length.
#     Boolean plot               - If true, will create plots for each alignment in the plot folder.
#
#   Runtime:
#     Int  mem                   - Amount of memory to give to the machine running each task in this workflow.
#     Int  preemptible_attempts  - Number of times to allow each task in this workflow to be preempted.
#     Int  disk_space_gb         - Amount of storage disk space (in Gb) to give to each machine running each task in this workflow.
#     Int  cpu                   - Number of CPU cores to give to each machine running each task in this workflow.
#     Int  boot_disk_size_gb     - Amount of boot disk space (in Gb) to give to each machine running each task in this workflow.
#
workflow Cartographer {
    input {
        File reads_fasta
        File segments_fasta
        File layout_order

        Float? min_qual
        Int? min_bases
        Float? prec_known
        Float? prec_unknown
        Boolean? use_mosaic_aligner
        Int? barcode_length
        File? barcode_file
        Boolean? plot

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    String docker_image = "us.gcr.io/broad-dsp-lrma/lr-cartographer:0.0.1"

    call SplitFastaFileTask {
        input:
            docker_image              = docker_image,
            reads_fasta               = reads_fasta,
    }

    scatter ( f in SplitFastaFileTask.reads_files ) {
        call CartographerTask {
            input:
                docker_image              = docker_image,

                reads_fasta               = f,
                segments_fasta            = segments_fasta,
                layout_order              = layout_order,

                min_qual                  = min_qual,
                min_bases                 = min_bases,
                prec_known                = prec_known,
                prec_unknown              = prec_unknown,
                use_mosaic_aligner        = use_mosaic_aligner,
                barcode_length            = barcode_length,
                barcode_file              = barcode_file,
                plot                      = plot,

                mem_gb                    = mem_gb,
                preemptible_attempts      = preemptible_attempts,
                disk_space_gb             = disk_space_gb,
                cpu                       = cpu,
                boot_disk_size_gb         = boot_disk_size_gb
        }
    }

    output {
      Array[File] cartographs = CartographerTask.cartograph
      Array[File] log_files   = CartographerTask.log_file
      Array[File] timing_info = CartographerTask.timing_info
    }
}

task CartographerTask {

    input {
        # ------------------------------------------------
        # Input args:
        # Required:

        # Runtime Options:
        String docker_image

        File reads_fasta
        File segments_fasta
        File layout_order

        Float? min_qual
        Int? min_bases
        Float? prec_known
        Float? prec_unknown
        Boolean? use_mosaic_aligner
        Int? barcode_length
        File? barcode_file
        Boolean? plot

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    # ------------------------------------------------
    # Process input args:

    String log_file_name = "cartographer.log"
    String timing_output_file = "cartographer.timingInformation.txt"
    String memory_log_file = "cartographer.memory_log.txt"

    String min_qual_arg = if defined(min_qual) then " --minqual " else ""
    String min_bases_arg = if defined(min_bases) then " --minbases " else ""
    String prec_known_arg = if defined(prec_known) then " --prec_known " else ""
    String prec_unknown_arg = if defined(prec_unknown) then " --prec_unknown " else ""
    String mosaic_aligner_arg = if defined(use_mosaic_aligner) then " --MOSAIC " else ""
    String barcode_length_arg = if defined(barcode_length) then " --barcode_length " else ""
    String barcode_file_arg = if defined(barcode_file) then " --barcode_file " else ""
    String plot_arg = if defined(plot) then " --plot " else ""

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    Float input_files_size_gb = size(reads_fasta, "GiB") + size(segments_fasta, "GiB") + size(layout_order, "GiB")

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 2048
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb

    # ------------------------------------------------
    # Run our command:
    command {
        # Set up memory logging daemon:
        while true ; do
            date
            date +%s
            cat /proc/meminfo
            sleep 10
        done >> ~{memory_log_file} &
        mem_pid=$!

#        set -e

        # Do the real work here:
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        /usr/bin/time -v /cartographer.py -v \
            -r ~{reads_fasta} \
            -s ~{segments_fasta} \
            -l ~{layout_order} \
            ~{min_qual_arg}~{default="" sep=" --minqual " min_qual} \
            ~{min_bases_arg}~{default="" sep=" --minbases " min_qual} \
            ~{prec_known_arg}~{default="" sep=" --prec_known " prec_known} \
            ~{prec_unknown_arg}~{default="" sep=" --prec_unknown " prec_unknown} \
            ~{mosaic_aligner_arg} \
            ~{barcode_length_arg}~{default="" sep=" --barcode_length " barcode_length} \
            ~{barcode_file_arg}~{default="" sep=" --barcode_file " barcode_file} \
            ~{plot_arg}~{default="" sep=" --plot " plot} \
            2>&1 | tee ~{log_file_name}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ~{timing_output_file}

        # Stop the memory daemon (and stop it hard):
        set +e
        kill -9 $mem_pid

        # Get and compute timing information:
        set +e
        elapsedTime=""
        which bc &> /dev/null ; bcr=$?
        which python3 &> /dev/null ; python3r=$?
        which python &> /dev/null ; pythonr=$?
        if [[ $bcr -eq 0 ]] ; then elapsedTime=`echo "scale=6;$endTime - $startTime" | bc`;
        elif [[ $python3r -eq 0 ]] ; then elapsedTime=`python3 -c "print( $endTime - $startTime )"`;
        elif [[ $pythonr -eq 0 ]] ; then elapsedTime=`python -c "print( $endTime - $startTime )"`;
        fi
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
      # Default output file name:
      File cartograph           = "cartograph_output.tsv"
      File log_file             = "${log_file_name}"
      File timing_info          = "${timing_output_file}"
      File memory_log           = "${memory_log_file}"
    }
 }

 task SplitFastaFileTask {
    input {
        # ------------------------------------------------
        # Input args:
        # Required:

        # Runtime Options:
        String docker_image

        File reads_fasta

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    # ------------------------------------------------
    # Process input args:

    String timing_output_file = "split_fasta_file.timingInformation.txt"

    Int split_size = 5400

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    Float input_files_size_gb = size(reads_fasta, "GiB")

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 2048
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb

    # ------------------------------------------------
    # Run our command:
    command {
        set -e
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        # The tool takes about 2 seconds per read, so we split the reads into ~6h chunks here:
        split -a4 --additional-suffix .fasta -d -l ~{split_size} ~{reads_fasta} reads_chunk_

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ~{timing_output_file}

        # Get and compute timing information:
        set +e
        elapsedTime=""
        which bc &> /dev/null ; bcr=$?
        which python3 &> /dev/null ; python3r=$?
        which python &> /dev/null ; pythonr=$?
        if [[ $bcr -eq 0 ]] ; then elapsedTime=`echo "scale=6;$endTime - $startTime" | bc`;
        elif [[ $python3r -eq 0 ]] ; then elapsedTime=`python3 -c "print( $endTime - $startTime )"`;
        elif [[ $pythonr -eq 0 ]] ; then elapsedTime=`python -c "print( $endTime - $startTime )"`;
        fi
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
      Array[File] reads_files   = glob("reads_chunk_*")
      File timing_info          = "${timing_output_file}"
    }
 }
