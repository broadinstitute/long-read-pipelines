version 1.0

task ExtractBoundedReadSectionsTask {

    meta {
        description : "Run extract_bounded_read_sections on a dataset to map out known and unknown segments of reads."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        # ------------------------------------------------
        # Input args:
        # Required:

        File reads_file
        File segments_fasta
        File boundaries_file

        String aligner = "BWA_MEM"

        Float min_qual = 7.0

        Int? max_read_length
        Int? min_bases
        Float? prec_known
        Float? prec_unknown

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    parameter_meta {
        reads_file : "SAM/BAM/FASTA/FASTQ file containing reads for which to determine the layout."
        segments_fasta : "FASTA file containing unique segments for which to search in the given BAM files.   These segments are used as delimiters in the reads.  Read splitting uses these delimiters and the boundaries file."
        boundaries_file : "Text file containing two comma-separated segment names from the segments_fasta on each line.  These entries define delimited sections to be extracted from the reads and treated as individual array elements."

        aligner : "[optional] The alignment method to use to find the bounded regions (TESSERAE, MOSAIC_ALIGNER, SMITH_WATERMAN, NEEDLEMAN_WUNSCH, BWA_MEM, BWA_ALN) (Default: BWA_MEM)"

        min_qual : "[optional] Minimum quality for good alignment.  (Default: 7.0)"

        max_read_length : "[optional] The read length beyond which a read will not be processed."
        min_bases : "[optional] Minimum number of bases for an alignment to be retained."
        prec_known : "[optional] Probability of recombination for known segment alignment."
        prec_unknown : "[optional] Probability of recombination for UNKNOWN segment alignment."

        mem_gb : "[optional] Amount of memory to give to the machine running each task in this workflow."
        preemptible_attempts : "[optional] Number of times to allow each task in this workflow to be preempted."
        disk_space_gb : "[optional] Amount of storage disk space (in Gb) to give to each machine running each task in this workflow."
        cpu : "[optional] Number of CPU cores to give to each machine running each task in this workflow."
        boot_disk_size_gb : "[optional] Amount of boot disk space (in Gb) to give to each machine running each task in this workflow."
    }

    # Docker image:
    String docker_image = "us.gcr.io/broad-dsp-lrma/lr-cartographer:0.0.3"

    # ------------------------------------------------
    # Process input args:

    String log_file_name = "extract_bounded_read_sections.log"
    String timing_output_file = "timingInformation.txt"
    String memory_log_file = "memory_log.txt"

    String max_read_length_arg = if defined(max_read_length) then " --max_read_length " else ""
    String min_bases_arg = if defined(min_bases) then " --minbases " else ""
    String prec_known_arg = if defined(prec_known) then " --prec_known " else ""
    String prec_unknown_arg = if defined(prec_unknown) then " --prec_unknown " else ""

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    Float input_files_size_gb = size(reads_file, "GiB") + size(segments_fasta, "GiB") + size(boundaries_file, "GiB")

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 8192
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb

    # ------------------------------------------------
    # Run our command:
    command {
        # Set up memory logging daemon:
        MEM_LOG_INTERVAL_s=5
        DO_MEMORY_LOG=true
        while $DO_MEMORY_LOG ; do
            date
            date +%s
            cat /proc/meminfo
            sleep $MEM_LOG_INTERVAL_s
        done >> ~{memory_log_file} &
        mem_pid=$!

        set -e

        # Do the real work here:
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        /usr/bin/time -v /cartographer/extract_bounded_read_sections.py -v \
            -r ~{reads_file} \
            -s ~{segments_fasta} \
            -b ~{boundaries_file} \
            ~{max_read_length_arg}~{default="" sep=" --max_read_length " max_read_length} \
            --minqual ~{min_qual} \
            ~{min_bases_arg}~{default="" sep=" --minbases " min_bases} \
            ~{prec_known_arg}~{default="" sep=" --prec_known " prec_known} \
            ~{prec_unknown_arg}~{default="" sep=" --prec_unknown " prec_unknown} \
            --aligner ~{aligner} \
            2>&1 | tee ~{log_file_name}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ~{timing_output_file}

        # Stop the memory daemon softly.  Then stop it hard if it's not cooperating:
        set +e
        DO_MEMORY_LOG=false
        sleep $(($MEM_LOG_INTERVAL_s  * 2))
        kill -0 $mem_pid &> /dev/null
        if [ $? -ne 0 ] ; then
            kill -9 $mem_pid
        fi

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
      File extracted_reads            = "extracted_bounded_sub_reads.fasta"
      File rejected_reads             = "extracted_bounded_sub_reads.rejected.fasta"
      File raw_marker_alignments      = "extracted_bounded_sub_reads.raw_marker_alignments.txt"
      File initial_section_alignments = "extracted_bounded_sub_reads.initial_section_alignments.txt"
      File final_section_alignments   = "extracted_bounded_sub_reads.final_section_alignments.txt"
      File log_file                   = "${log_file_name}"
      File timing_info                = "${timing_output_file}"
      File memory_log                 = "${memory_log_file}"
    }
 }

task SplitSequenceOnDelimiters {

    meta {
        description : "Run split_sequence_on_delimiters on a dataset to split each read by the given delimiters."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        # ------------------------------------------------
        # Input args:
        # Required:

        File reads_file
        File delimiters_fasta

        Float min_qual = 7.0

        Int? max_read_length
        Int? min_bases

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    parameter_meta {
        reads_file : "SAM/BAM/FASTA/FASTQ file containing reads for which to determine the layout."
        delimiters_fasta : "FASTA file containing unique delimiters by which to split the reads in the given reads_file."

        min_qual : "[optional] Minimum quality for good alignment.  (Default: 7.0)"

        max_read_length : "[optional] The read length beyond which a read will not be processed."
        min_bases : "[optional] Minimum number of bases for an alignment to be retained."

        mem_gb : "[optional] Amount of memory to give to the machine running each task in this workflow."
        preemptible_attempts : "[optional] Number of times to allow each task in this workflow to be preempted."
        disk_space_gb : "[optional] Amount of storage disk space (in Gb) to give to each machine running each task in this workflow."
        cpu : "[optional] Number of CPU cores to give to each machine running each task in this workflow."
        boot_disk_size_gb : "[optional] Amount of boot disk space (in Gb) to give to each machine running each task in this workflow."
    }

    # Docker image:
    String docker_image = "us.gcr.io/broad-dsp-lrma/lr-cartographer:0.0.3"

    # ------------------------------------------------
    # Process input args:

    String log_file_name = "split_sequence_on_delimiters.log"
    String timing_output_file = "timingInformation.txt"
    String memory_log_file = "memory_log.txt"

    String max_read_length_arg = if defined(max_read_length) then " --max_read_length " else ""
    String min_bases_arg = if defined(min_bases) then " --minbases " else ""

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    Float input_files_size_gb = size(reads_file, "GiB") + size(delimiters_fasta, "GiB")

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 8192
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb

    # ------------------------------------------------
    # Run our command:
    command {
        # Set up memory logging daemon:
        MEM_LOG_INTERVAL_s=5
        DO_MEMORY_LOG=true
        while $DO_MEMORY_LOG ; do
            date
            date +%s
            cat /proc/meminfo
            sleep $MEM_LOG_INTERVAL_s
        done >> ~{memory_log_file} &
        mem_pid=$!

        set -e

        # Do the real work here:
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        /usr/bin/time -v /cartographer/split_sequence_on_delimiters.py -v \
            -r ~{reads_file} \
            -d ~{delimiters_fasta} \
            ~{max_read_length_arg}~{default="" sep=" --max_read_length " max_read_length} \
            --minqual ~{min_qual} \
            ~{min_bases_arg}~{default="" sep=" --minbases " min_bases} \
            2>&1 | tee ~{log_file_name}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ~{timing_output_file}

        # Stop the memory daemon softly.  Then stop it hard if it's not cooperating:
        set +e
        DO_MEMORY_LOG=false
        sleep $(($MEM_LOG_INTERVAL_s  * 2))
        kill -0 $mem_pid &> /dev/null
        if [ $? -ne 0 ] ; then
            kill -9 $mem_pid
        fi

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
      File split_sequences            = "split_sequences.bam"
      File rejected_reads             = "split_sequences.rejected.bam"
      File seq_annotations            = "split_sequences.annotations.tsv"
      File log_file                   = "${log_file_name}"
      File timing_info                = "${timing_output_file}"
      File memory_log                 = "${memory_log_file}"
    }
 }

 task GetApproxRawSubreadArrayLengths {

    meta {
        description : "Get the approximate array lengths from raw PacBio subreads."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        # ------------------------------------------------
        # Input args:
        # Required:

        File reads_file
        File delimiters_fasta

        Float min_qual = 7.0

        Array[String]? ignore_seqs

        String prefix = "approx_raw_subread_array_lengths"

        String alignment_algorithm = "BWA_ALN"

        Int? max_read_length
        Int? min_bases

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    parameter_meta {
        reads_file : "SAM/BAM/FASTA/FASTQ file containing reads for which to determine the layout."
        delimiters_fasta : "FASTA file containing unique delimiters by which to split the reads in the given reads_file."

        min_qual : "[optional] Minimum quality for good alignment.  (Default: 7.0)"

        ignore_seqs : "[optional] Sequences to ignore and not use as delimiters."

        prefix : "[optional] Base name for the output data file.  (Default: approx_raw_subread_array_lengths)"

        alignment_algorithm : "[optional] Alignment to use to discover delimiters.  (Default: BWA_ALN)"

        max_read_length : "[optional] The read length beyond which a read will not be processed."
        min_bases : "[optional] Minimum number of bases for an alignment to be retained."

        mem_gb : "[optional] Amount of memory to give to the machine running each task in this workflow."
        preemptible_attempts : "[optional] Number of times to allow each task in this workflow to be preempted."
        disk_space_gb : "[optional] Amount of storage disk space (in Gb) to give to each machine running each task in this workflow."
        cpu : "[optional] Number of CPU cores to give to each machine running each task in this workflow."
        boot_disk_size_gb : "[optional] Amount of boot disk space (in Gb) to give to each machine running each task in this workflow."
    }

    # Docker image:
    String docker_image = "us.gcr.io/broad-dsp-lrma/lr-cartographer:0.0.4"

    # ------------------------------------------------
    # Process input args:

    String ignore_seqs_arg = if defined(ignore_seqs) then " --ignore " else ""

    String log_file_name = "get_approx_raw_subread_array_lengths.log"
    String timing_output_file = "timingInformation.txt"
    String memory_log_file = "memory_log.txt"

    String max_read_length_arg = if defined(max_read_length) then " --max_read_length " else ""
    String min_bases_arg = if defined(min_bases) then " --minbases " else ""

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    Float input_files_size_gb = size(reads_file, "GiB") + size(delimiters_fasta, "GiB")

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 8192
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb

    # ------------------------------------------------
    # Run our command:
    command {
        # Set up memory logging daemon:
        MEM_LOG_INTERVAL_s=5
        DO_MEMORY_LOG=true
        while $DO_MEMORY_LOG ; do
            date
            date +%s
            cat /proc/meminfo
            sleep $MEM_LOG_INTERVAL_s
        done >> ~{memory_log_file} &
        mem_pid=$!

        set -e

        # Do the real work here:
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        /usr/bin/time -v /cartographer/get_approx_raw_subread_array_lengths.py -v \
            -b ~{reads_file} \
            -d ~{delimiters_fasta} \
            ~{ignore_seqs_arg}~{default="" sep=" --ignore " ignore_seqs} \
            -o ~{prefix}.tsv \
            --aligner ~{alignment_algorithm} \
            --minqual ~{min_qual} \
            2>&1 | tee ~{log_file_name}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ~{timing_output_file}

        # Stop the memory daemon softly.  Then stop it hard if it's not cooperating:
        set +e
        DO_MEMORY_LOG=false
        sleep $(($MEM_LOG_INTERVAL_s  * 2))
        kill -0 $mem_pid &> /dev/null
        if [ $? -ne 0 ] ; then
            kill -9 $mem_pid
        fi

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
      File approx_subread_array_lengths = "~{prefix}.tsv"
      File log_file                     = "${log_file_name}"
      File timing_info                  = "${timing_output_file}"
      File memory_log                   = "${memory_log_file}"
    }
 }