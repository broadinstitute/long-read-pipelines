version 1.0

task FlairQuant {

    meta {
        description : "Quantify transcript isoforms using FLAIR quantify."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        # ------------------------------------------------
        # Input args:
        # Required:

        File fasta_tar_gz
        File transcript_isoforms_fasta
        String out_base_name = "reads"

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    parameter_meta {
        fasta_tar_gz : "tar.gz file contianing FASTA files, each with reads from the same sample / cell barcode to quantify for transcripts.  Corresponds to `fasta_tar_gz_out` from Preprocessing_Tasks.SplitBamBySampleAndCellBarcodeTask."
        transcript_isoforms_fasta : "FASTA file containing isoforms sequences to quantify."
        out_base_name : "[optional] Base name for the resulting output file.  Ideally should be a unique identifier for this dataset."

        mem_gb : "[optional] Amount of memory to give to the machine running each task in this workflow."
        preemptible_attempts : "[optional] Number of times to allow each task in this workflow to be preempted."
        disk_space_gb : "[optional] Amount of storage disk space (in Gb) to give to each machine running each task in this workflow."
        cpu : "[optional] Number of CPU cores to give to each machine running each task in this workflow."
        boot_disk_size_gb : "[optional] Amount of boot disk space (in Gb) to give to each machine running each task in this workflow."
    }

    # ------------------------------------------------
    # Process input args:

    String timing_output_file = "timingInformation.txt"
    String memory_log_file = "memory_log.txt"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # Triple our input file size for our disk size - should be plenty.
    Float input_files_size_gb = 3*(size(fasta_tar_gz, "GiB") + size(transcript_isoforms_fasta, "GiB"))

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 4096
    Int default_disk_space_gb = ceil((input_files_size_gb * 2) + 1024)

    Int default_boot_disk_size_gb = 15

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb

    # ------------------------------------------------
    # Run our command:
    command <<<
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

        echo "Extracting fasta files from: ~{fasta_tar_gz} ..."
        tar -zxf ~{fasta_tar_gz}
        echo "Done."

        echo "Generating reads manifest file..."
        reads_manifest_file="flair_reads_manifest.tsv"
        ls *.fasta | grep -v '~{transcript_isoforms_fasta}' | awk -F"_" 'BEGIN{OFS="\t"}{z=$0;sub(".*/", "", $0);sub(".fasta", "", $0); print $2,"A",$3,z}' > $reads_manifest_file
        echo "Done"

        echo "Quantifying flair results..."

        python /flair/flair.py quantify -r $reads_manifest_file -i ~{transcript_isoforms_fasta} -o ~{out_base_name}.flair.counts_matrix.tsv

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
    >>>

    # ------------------------------------------------
    # Runtime settings:
     runtime {
         docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.1"
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
      File count_matrix         = "${out_base_name}.flair.counts_matrix.tsv"

      File timing_info          = "${timing_output_file}"
      File memory_log           = "${memory_log_file}"
    }
}

task AlignCellDropletReadsTask {
    input {
        File read_file
        File transcripts_fasta
    }

    String out_file = basename(read_file, ".fasta") + ".sam"

    # 10x for the decompressed file size
    # 2x for potential for FASTQ and SAM files (from file conversion).
    Int disk_size = 10*ceil(size(read_file, "GB"))*2 + ceil(size(transcripts_fasta, "GB"))

    command <<<
        minimap2 -a -N 4 -t 4 ~{transcripts_fasta} ~{read_file} > ~{out_file}
    >>>

    output {
        File aligned_sam = "~{out_file}"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
        memory: 8 + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        boot_disk_gb: 10
        preemptible: 0
        cpu: 4
    }
}

task DiscoverIsoformsTask {

    meta {
        description : "Analyze transcript isoforms with Flair.  (https://github.com/BrooksLabUCSC/flair)"
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        # ------------------------------------------------
        # Input args:
        # Required:

        File unaligned_reads
        File aligned_reads
        File gtf_annotations

        File reference_fasta
        File reference_fasta_index
        File reference_seq_dict

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    parameter_meta {
        unaligned_reads : "FASTA/FASTQ file containing unaligned reads.  This should correspond to the given aligned_reads file."
        aligned_reads : "SAM/BAM file containing aligned reads.  This should correspond to the given unaligned_reads file."
        gtf_annotations : "GTF file containing gene annotations for the reference to which the given reads are aligned."

        reference_fasta : "FASTA file containing the reference sequence."
        reference_fasta_index : "Index for the given reference_fasta file."
        reference_seq_dict : "Sequence Dictionary for the given reference_fasta file."

        mem_gb : "[optional] Amount of memory to give to the machine running each task in this workflow."
        preemptible_attempts : "[optional] Number of times to allow each task in this workflow to be preempted."
        disk_space_gb : "[optional] Amount of storage disk space (in Gb) to give to each machine running each task in this workflow."
        cpu : "[optional] Number of CPU cores to give to each machine running each task in this workflow."
        boot_disk_size_gb : "[optional] Amount of boot disk space (in Gb) to give to each machine running each task in this workflow."
    }

    # Docker image:
    String docker_image = "us.gcr.io/broad-dsp-lrma/lr-cartographer:0.0.1"

    # ------------------------------------------------
    # Process input args:

    String out_base_name = sub(sub(basename(aligned_reads), "\.bam$", ""), "\.sam$", "")

    String log_file_name = "flair.log"
    String timing_output_file = "timingInformation.txt"
    String memory_log_file = "memory_log.txt"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    Float input_files_size_gb = size(unaligned_reads, "GiB") + size(aligned_reads, "GiB") + size(gtf_annotations, "GiB") + size(reference_fasta, "GiB") + size(reference_fasta_index, "GiB") + size(reference_seq_dict, "GiB")

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

        echo "Converting reads to bed file: ~{aligned_reads} ..."
        python /flair/bin/bam2Bed12.py -i ~{aligned_reads} > ~{out_base_name}.bed12

        echo "Correcting bed file regions ..."
        python /flair/flair.py correct \
            -q ~{out_base_name}.bed12 \
            -g ~{reference_fasta} \
            -f ~{gtf_annotations} \
            -o ~{out_base_name}.flair

        echo "Collapsing the corrected regions into results ..."
        python /flair/flair.py collapse \
            -q ~{out_base_name}.flair_all_corrected.bed \
            -r ~{unaligned_reads} \
            -g ~{reference_fasta} \
            -f ~{gtf_annotations} \
            -o ~{out_base_name}.flair.collapse

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
      File flair_bed12                     = "${out_base_name}.bed12"
      File flair_isoforms_inconsistent_bed = "${out_base_name}.flair_all_inconsistent.bed"
      File flair_isoforms_corrected_bed    = "${out_base_name}.flair_all_corrected.bed"
      File flair_isoforms_collapsed_bed    = "${out_base_name}.flair.collapse.isoforms.bed"
      File flair_isoforms_collapsed_fasta  = "${out_base_name}.flair.collapse.isoforms.fa"
      File flair_isoforms_collapsed_gtf    = "${out_base_name}.flair.collapse.isoforms.gtf"

      File log_file             = "${log_file_name}"
      File timing_info          = "${timing_output_file}"
      File memory_log           = "${memory_log_file}"
    }
 }

task QuantifyAndVisualizeFlairResults {

    meta {
        description : "Quantify and Visualize results from FLAIR tool output (i.e. from AnalyzeWithFlairTask).  (https://github.com/BrooksLabUCSC/flair)"
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        # ------------------------------------------------
        # Input args:
        # Required:

        File unaligned_reads
        File flair_collapsed_isoforms_fasta
        File flair_collapsed_isoforms_bed

        File gtf_annotations
        File reference_fasta
        File reference_fasta_index
        File reference_seq_dict

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    parameter_meta {
        unaligned_reads : "FASTA/FASTQ file containing unaligned reads.  This should correspond to the given aligned_reads file."
        flair_collapsed_isoforms_fasta : "A collapsed isoforms fasta file created by `flair collapse`."
        flair_collapsed_isoforms_bed : "A collapsed isoforms bed file created by `flair collapse`."

        gtf_annotations : "GTF file containing gene annotations for the reference to which the given reads are aligned."
        reference_fasta : "FASTA file containing the reference sequence."
        reference_fasta_index : "Index for the given reference_fasta file."
        reference_seq_dict : "Sequence Dictionary for the given reference_fasta file."

        mem_gb : "[optional] Amount of memory to give to the machine running each task in this workflow."
        preemptible_attempts : "[optional] Number of times to allow each task in this workflow to be preempted."
        disk_space_gb : "[optional] Amount of storage disk space (in Gb) to give to each machine running each task in this workflow."
        cpu : "[optional] Number of CPU cores to give to each machine running each task in this workflow."
        boot_disk_size_gb : "[optional] Amount of boot disk space (in Gb) to give to each machine running each task in this workflow."
    }

    # Docker image:
    String docker_image = "us.gcr.io/broad-dsp-lrma/lr-cartographer:0.0.1"

    # ------------------------------------------------
    # Process input args:

    String out_base_name = sub(sub(basename(unaligned_reads), "\.bam$", ""), "\.sam$", "")

    String log_file_name = "flair.log"
    String timing_output_file = "timingInformation.txt"
    String memory_log_file = "memory_log.txt"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    Float input_files_size_gb = size(unaligned_reads, "GiB") + size(flair_collapsed_isoforms_fasta, "GiB") + size(flair_collapsed_isoforms_bed, "GiB") + size(gtf_annotations, "GiB") + size(reference_fasta, "GiB") + size(reference_fasta_index, "GiB") + size(reference_seq_dict, "GiB")

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 4096
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

        echo "Quantifying flair results..."
        reads_manifest_file=flair_reads_mainfest.tsv
        echo -e "~{out_base_name}\tcondition_1\tbatch_1\t$~{unaligned_reads}" > $reads_manifest_file
        python /flair/flair.py quantify -r $reads_manifest_file -i ~{flair_collapsed_isoforms_fasta} -o ~{out_base_name}.flair.counts_matrix.tsv

        echo "Generating differential expression graphs..."
        python /flair/flair.py diffExp -q ~{out_base_name}.flair.counts_matrix.tsv -o ~{out_base_name}_flair_diff_exp -of
        tar -zcf ~{out_base_name}_flair_diff_exp.tar.gz ~{out_base_name}_flair_diff_exp

        echo "Calling diff splice..."
        python /flair/flair.py diffSplice -i ~{flair_collapsed_isoforms_bed} -q ~{out_base_name}.flair.counts_matrix.tsv -o ~{out_base_name}.flair.diffsplice

        echo "Calling Predict Productivity..."
        python /flair/bin/predictProductivity.py -i ~{flair_collapsed_isoforms_bed} -g ~{gtf_annotations} -f ~{reference_fasta} --longestORF > ~{out_base_name}.productivity.bed

        echo "Calling Mark Intron Retention..."
        python /flair/bin/mark_intron_retention.py ~{flair_collapsed_isoforms_bed} ~{out_base_name}.intron_retention_out_isoforms.psl ~{out_base_name}.intron_retention_out_coords.txt

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
      File counts_matrix                   = "${out_base_name}.flair.counts_matrix.tsv"
      File differential_expression_tar_gz  = "${out_base_name}_flair_diff_exp.tar.gz"
      File diff_splice_isoforms            = "${out_base_name}.flair.diffsplice"
      File productivity_bed                = "${out_base_name}.productivity.bed"
      File intron_retention_out_isoforms   = "${out_base_name}.intron_retention_out_isoforms.psl"
      File intron_retention_out_coords     = "${out_base_name}.intron_retention_out_coords.txt"

      File log_file             = "${log_file_name}"
      File timing_info          = "${timing_output_file}"
      File memory_log           = "${memory_log_file}"
    }
 }


task FlairPlotIsoformUsage {

    meta {
        description : "Plot Isoforms Usage for a Gene from Flair Results.  (https://github.com/BrooksLabUCSC/flair)"
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }


    input {
        # ------------------------------------------------
        # Input args:
        # Required:

        String gene_name
        File flair_collapsed_isoforms_bed
        File flair_quantified_counts_matrix

        File gtf_annotations
        File reference_fasta
        File reference_fasta_index
        File reference_seq_dict

        Int? mem_gb
        Int? preemptible_attempts
        Int? disk_space_gb
        Int? cpu
        Int? boot_disk_size_gb
    }

    parameter_meta {
        gene_name : "The name of the gene for which to create plots."
        flair_collapsed_isoforms_bed : "A collapsed isoforms bed file created by `flair collapse`."
        flair_quantified_counts_matrix : "A count matrix file created by `flair quantify`."

        gtf_annotations : "GTF file containing gene annotations for the reference to which the given reads are aligned."
        reference_fasta : "FASTA file containing the reference sequence."
        reference_fasta_index : "Index for the given reference_fasta file."
        reference_seq_dict : "Sequence Dictionary for the given reference_fasta file."

        mem_gb : "[optional] Amount of memory to give to the machine running each task in this workflow."
        preemptible_attempts : "[optional] Number of times to allow each task in this workflow to be preempted."
        disk_space_gb : "[optional] Amount of storage disk space (in Gb) to give to each machine running each task in this workflow."
        cpu : "[optional] Number of CPU cores to give to each machine running each task in this workflow."
        boot_disk_size_gb : "[optional] Amount of boot disk space (in Gb) to give to each machine running each task in this workflow."
    }

    # Docker image:
    String docker_image = "us.gcr.io/broad-dsp-lrma/lr-cartographer:0.0.1"

    # ------------------------------------------------
    # Process input args:

    String log_file_name = "flair.log"
    String timing_output_file = "timingInformation.txt"
    String memory_log_file = "memory_log.txt"

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    Float input_files_size_gb = size(flair_quantified_counts_matrix, "GiB") + size(flair_collapsed_isoforms_bed, "GiB")

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 4096
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

        python plot_isoform_usage.py ~{flair_collapsed_isoforms_bed} ~{flair_quantified_counts_matrix} ~{gene_name}

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
      File isoforms_plot        = "${gene_name}_isoforms.png"
      File usage_plot           = "${gene_name}_usage.png"

      File log_file             = "${log_file_name}"
      File timing_info          = "${timing_output_file}"
      File memory_log           = "${memory_log_file}"
    }
 }
