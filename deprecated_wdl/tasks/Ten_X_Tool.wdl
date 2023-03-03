version 1.0

import "Utils.wdl" as Utils
import "PBUtils.wdl" as PB
import "Structs.wdl"

# TODO: Merge this file and `AnnotateAdapters.wdl`

workflow AnnotateBarcodesAndUMIsWorkflow {

    meta {
        description : "Run the 10x tool on a bam file containing 10x library prepared reads to annotate each read with the adapter locations."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File bam_file
        File? bam_index
        File head_adapter_fasta
        File tail_adapter_fasta
        Int read_end_length

        File? whitelist_10x
        File? whitelist_illumina

        Int? poly_t_length
        Int? barcode_length
        Int? umi_length

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {

        bam_file : "BAM file containing reads from a 10x library preparation as sequenced off of a PacBio Sequel II."
        bam_index : "[optional] Index file for the given bam_file."

        head_adapter_fasta : "FASTA file containing the sequence that each transcript should start with.  Typically this will be the 10x adapter sequence from the 10x library prep."
        tail_adapter_fasta : "FASTA file containing the sequence that each transcript should end with.  Typically this will be the Template Switch Oligo (TSO) sequence from the 10x library prep."
        read_end_length : "Number of bases at either end of the read to check for the adapters."

        whitelist_10x : "[optional] Barcode whitelist from the 10x library prep.  If provided, only reads matching one of these barcodes will be annotated and output."
        whitelist_illumina : "[optional] Additional barcode whitelist from a parallel Illumina sequencing run.  If this option is provided, you must also supply the `whitelist_10x` barcode file.  When both this and the `whitelist_10x` barcode file are supplied, only reads matching barcodes in this file will be annotated and output."

        poly_t_length : "[optional] Expected length of the poly-T region in this library preparation (default: 10)."
        barcode_length : "[optional] Length of the barcode region in this library preparation (default: 16)."
        umi_length : "[optional] Length of the UMI region in this library preparation (default: 12)."

        runtime_attr_override : "[optional] Runtime attributes struct with which to override the docker container runtime.."
    }

    call PB.PBIndex as PBIndexSubreadShard { input: bam = bam_file }
    call PB.ShardLongReads {
        input:
            unaligned_bam = bam_file,
            unaligned_pbi = PBIndexSubreadShard.pbi,
    }

    scatter (reads_file in ShardLongReads.unmapped_shards) {
        call AnnotateBarcodesAndUMIs {
            input:
                bam_file = reads_file,
                head_adapter_fasta = head_adapter_fasta,
                tail_adapter_fasta = tail_adapter_fasta,
                read_end_length = read_end_length,
                whitelist_10x = whitelist_10x,
                whitelist_illumina = whitelist_illumina,
                poly_t_length = poly_t_length,
                barcode_length = barcode_length,
                umi_length = umi_length,
                runtime_attr_override = runtime_attr_override
        }
    }
    
    output {
        Array[File] output_bam        = AnnotateBarcodesAndUMIs.output_bam
        Array[File] barcode_stats     = AnnotateBarcodesAndUMIs.barcode_stats
        Array[File] starcode          = AnnotateBarcodesAndUMIs.starcode
        Array[File] stats             = AnnotateBarcodesAndUMIs.stats
        Array[File] timing_info       = AnnotateBarcodesAndUMIs.timing_info
    }
}

task AnnotateBarcodesAndUMIs {
    input {
        File bam_file
        File? bam_index
        File head_adapter_fasta
        File tail_adapter_fasta
        Int read_end_length

        Boolean raw_extract_only = false

        File? whitelist_10x
        File? whitelist_illumina

        File? illumina_barcoded_bam

        Int? poly_t_length
        Int? barcode_length
        Int? umi_length

        RuntimeAttr? runtime_attr_override
    }

    # ------------------------------------------------
    # Set runtime options:
    String whitelist_10x_arg = if defined(whitelist_10x) then " --whitelist-10x " else ""
    String whitelist_ilmn_arg = if defined(whitelist_illumina) then " --whitelist-illumina " else ""

    String illumina_barcoded_bam_arg = if defined(illumina_barcoded_bam) then " --illumina-bam " else ""

    String poly_t_len_arg = if defined(poly_t_length) then " --poly-t-length " else ""
    String barcode_len_arg = if defined(barcode_length) then " --barcode-length " else ""
    String umi_len_arg = if defined(umi_length) then " --umi-length " else ""

    String do_raw_arg = if raw_extract_only then " --raw " else ""

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int disk_size_gb = 40 + ceil(size(bam_file, "GiB")* 8)

    String timing_output_file = "timingInformation.txt"

    String output_name = basename(bam_file) + ".10x_annotated"

    String memory_log_file = "memory_use.txt"

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
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        source activate 10x_tool

        python /lrma/tool.py \
            --bam=~{bam_file} \
            --adapter=~{head_adapter_fasta} \
            --reverse-adapter=~{tail_adapter_fasta} \
            --name=~{output_name} \
            --read-end-length=~{read_end_length} \
            --record-umis \
            ~{do_raw_arg} \
            ~{whitelist_10x_arg}~{default="" sep=" --whitelist-10x " whitelist_10x} \
            ~{whitelist_ilmn_arg}~{default="" sep=" --whitelist-illumina " whitelist_illumina} \
            ~{illumina_barcoded_bam_arg}~{default="" sep=" --illumina-bam " illumina_barcoded_bam} \
            ~{poly_t_len_arg}~{default="" sep=" --poly-t-length " poly_t_length} \
            ~{barcode_len_arg}~{default="" sep=" --barcode-length " barcode_length} \
            ~{umi_len_arg}~{default="" sep=" --umi-length " umi_length}

        endTime=`date +%s.%N`
        echo "EndTime: $endTime" >> ~{timing_output_file}

        # Make sure that the reqired output files all exist:
        if [ ! -e "~{output_name}_starcode_confidence_factor_barcode_counts.tsv" ] ; then
            touch ~{output_name}_starcode_confidence_factor_barcode_counts.tsv
        fi
        if [ ! -e "~{output_name}_starcode.tsv" ] ;then
            touch ~{output_name}_starcode.tsv
        fi

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

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             16,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-10x:0.1.14"
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


    output {
      File output_bam          = "${output_name}.bam"
      File barcode_stats       = "${output_name}_barcode_stats.tsv"
      File starcode            = "${output_name}_starcode.tsv"
      File stats               = "${output_name}_stats.tsv"
      File raw_starcode_counts = "${output_name}_starcode_confidence_factor_barcode_counts.tsv"
      File memory_info         = "${memory_log_file}"
      File timing_info         = "${timing_output_file}"
    }
}

task CorrectBarcodesWithStarcodeSeedCounts {
    input {
        File bam_file
        File starcode_seeds_tsv

        String prefix = "sample"

        File? whitelist_10x
        String? extra_parameters

        RuntimeAttr? runtime_attr_override
    }

    String whitelist_10x_arg = if defined(whitelist_10x) then " --whitelist-10x " else ""

    # ------------------------------------------------
    # Get machine settings:

    # You may have to change the following two parameter values depending on the task requirements
    Int disk_size_gb = ceil(( (size(bam_file, "GiB") + size(starcode_seeds_tsv, "GiB")) * 8) + 40)

    String timing_output_file = "timingInformation.txt"
    String memory_log_file = "memory_use.txt"

    String extra_params_arg = if defined(extra_parameters) then " ~{extra_parameters} " else ""

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
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        source activate 10x_tool

        python /lrma/tool_starcode_seeded.py \
            --bam=~{bam_file} \
            --counts=~{starcode_seeds_tsv} \
            --name=~{prefix} \
            ~{extra_params_arg} \
            ~{whitelist_10x_arg}~{default="" sep=" --whitelist-10x " whitelist_10x}


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

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             16,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-10x:0.1.16"
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

    output {
      File output_bam          = "${prefix}_corrected_barcodes.bam"
      File memory_info         = "${memory_log_file}"
      File timing_info         = "${timing_output_file}"
    }
}

task ExtractIlmnBarcodeConfScores {
    input {
        File bam_file

        String prefix = "sample"

        RuntimeAttr? runtime_attr_override
    }

    # You may have to change the following two parameter values depending on the task requirements
    Int disk_size_gb = ceil(size(bam_file, "GiB") * 8) + 40

    String timing_output_file = "timingInformation.txt"
    String memory_log_file = "memory_use.txt"

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
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        source activate 10x_tool

        python /lrma/extract_ilmn_bc_conf_scores.py \
            --bam=~{bam_file} \
            --prefix=~{prefix}

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

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             16,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-10x:0.1.14"
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

    output {
      File conf_score_tsv = "${prefix}_ilmn_barcode_conf_scores.tsv"
      File memory_info    = "${memory_log_file}"
      File timing_info    = "${timing_output_file}"
    }
}



task ExtractCbcAndUmiFromAnnotatedReadForUmiTools {
    input {
        File annotated_bam_file

        Int? boot_disk_size_gb
        Int? cpu
        Int? disk_space_gb
        Int? mem_gb
        Int? preemptible_attempts
    }

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 16 * 1024

    Float reads_size_gb = size(annotated_bam_file, "GiB")
    Int default_disk_space_gb = ceil((reads_size_gb * 2) + 20)

    Int default_boot_disk_size_gb = 10

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb

    String timing_output_file = "timingInformation.txt"

    String output_name = basename(annotated_bam_file, ".bam") + ".cbc_umi_extracted.bam"

    command {
        set -e
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        source activate 10x_tool

        python3 /lrma/extract_cbc_and_umi_from_annotated_read.py \
            --bam=~{annotated_bam_file} \
            --out-name=~{output_name}

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
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-10x:0.1.13"
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: 0
        cpu: select_first([cpu, 1])
    }
    output {
      File output_bam        = "~{output_name}"
      File timing_info       = "~{timing_output_file}"
    }
}

task RestoreAnnotationstoAlignedBam {
    input {
        File annotated_bam_file
        File aligned_bam_file

        Array[String] tags_to_ignore = ["RG"]

        Int? boot_disk_size_gb
        Int? cpu
        Int? disk_space_gb
        Int? mem_gb
        Int? preemptible_attempts
    }

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 32 * 1024

    Float reads_size_gb = size(annotated_bam_file, "GiB") + size(aligned_bam_file, "GiB")
    Int default_disk_space_gb = 10 * ceil((reads_size_gb * 10) + 20)

    Int default_boot_disk_size_gb = 100

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb

    String timing_output_file = "timingInformation.txt"

    String memory_log_file = "memory_use.txt"
    String output_name = basename(aligned_bam_file, ".bam") + ".AnnotationsRestored.bam"

    String ignore_tags_arg = if (length(tags_to_ignore) != 0 ) then "--ignore-tags " else ""

    command <<<

        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

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
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        source activate 10x_tool

        python3 /lrma/restore_annotations_to_aligned_bam.py \
            --bam ~{annotated_bam_file} \
            --aligned-bam ~{aligned_bam_file} \
            ~{ignore_tags_arg} ~{default="" sep=" " tags_to_ignore} \
            --out-name ~{output_name}

        samtools index -@${np} ~{output_name}

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

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-10x:0.1.15"
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: 0
        cpu: select_first([cpu, 1])
    }
    output {
        File timing_info       = "~{timing_output_file}"
        File memory_info       = "~{memory_log_file}"
        File output_bam        = "~{output_name}"
        File output_bam_index  = "~{output_name}.bai"
    }
}

task CopyContigNameToReadTag {
    meta {
        description : "Copies the contig name to a tag inside the read (tag: XG).  (This is useful for umi-tools.)"
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }
    input {
        File aligned_bam_file
        String prefix = "out"

        Int? boot_disk_size_gb
        Int? cpu
        Int? disk_space_gb
        Int? mem_gb
        Int? preemptible_attempts
    }

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 16 * 1024

    Float reads_size_gb = size(aligned_bam_file, "GiB")
    Int default_disk_space_gb = ceil((reads_size_gb * 10) + 20)

    Int default_boot_disk_size_gb = 10

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb

    String timing_output_file = "timingInformation.txt"

    String output_name = prefix + ".bam"

    command <<<
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        set -e
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        source activate 10x_tool

        python3 /lrma/copy_contig_name_to_tag.py \
            -b ~{aligned_bam_file} \
            -o ~{output_name}

        samtools index -@${np} ~{output_name}

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
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-10x:0.1.13"
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: 0
        cpu: select_first([cpu, 1])
    }
    output {
      File output_bam        = "~{output_name}"
      File output_bam_index  = "~{output_name}.bai"
      File timing_info       = "~{timing_output_file}"
    }
}


task TagSirvUmiPositionsFromLongbowAnnotatedArrayElement {
    meta {
        description : "Extracts the UMI from each read in the given bam file into the ZU tag."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }
    input {
        File bam_file
        String prefix = "out"

        Int? boot_disk_size_gb
        Int? cpu
        Int? disk_space_gb
        Int? mem_gb
        Int? preemptible_attempts
    }

    # ------------------------------------------------
    # Get machine settings:
    Boolean use_ssd = false

    # You may have to change the following two parameter values depending on the task requirements
    Int default_ram_mb = 16 * 1024

    Float reads_size_gb = size(bam_file, "GiB")
    Int default_disk_space_gb = ceil((reads_size_gb * 2) + 20)

    Int default_boot_disk_size_gb = 10

    # Mem is in units of GB but our command and memory runtime values are in MB
    Int machine_mem = if defined(mem_gb) then mem_gb * 1024 else default_ram_mb

    String timing_output_file = "timingInformation.txt"

    String output_name = prefix + ".bam"

    command {
        set -e
        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        source activate 10x_tool

        python3 /lrma/tag_mas_sirv_umi_positions.py \
            -b ~{bam_file} \
            -o ~{output_name}

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
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-10x:0.1.13"
        memory: machine_mem + " MB"
        disks: "local-disk " + select_first([disk_space_gb, default_disk_space_gb]) + if use_ssd then " SSD" else " HDD"
        bootDiskSizeGb: select_first([boot_disk_size_gb, default_boot_disk_size_gb])
        preemptible: 0
        cpu: select_first([cpu, 1])
    }
    output {
      File output_bam        = "~{output_name}"
      File timing_info       = "~{timing_output_file}"
    }
}
