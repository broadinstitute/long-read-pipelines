version 1.0

import "../Structs.wdl"

task RunSalmonQuantTask {

    meta {
        description : "Quantify transcripts from RNA transcript reads using SALMON."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File reads_fasta
        File salmon_index_tar_gz

        String? extra_args
    }

    parameter_meta {
        reads_fasta : "FASTA/FASTQ file containing unaligned reads."
        salmon_index_tar_gz : "SALMON index file corresponding to the transcripts FASTA file used to quantify the transcripts in the reads (TAR.GZ format)."

        extra_args : "[optional] Extra arguments to pass to SALMON after the rest of the command line has been given."
    }

    String out_prefix = basename(reads_fasta, ".fasta")

    String extra_args_arg = select_first([extra_args, ""])

    # 10x for the decompressed file size
    # 3x for salmon index files
    # 10 for baseline size
    Int disk_size = 10*ceil(size(reads_fasta, "GB"))*2 + (3 * ceil(size(salmon_index_tar_gz, "GB"))) + 10

    command <<<
        index_folder=$( tar -tvf ~{salmon_index_tar_gz}  | head -n1 | awk '{print $NF}' )
        tar -xf ~{salmon_index_tar_gz}

        salmon quant -lSR --dumpEq  -i $index_folder -r ~{reads_fasta} --validateMappings -o out_folder ~{extra_args_arg}

        # In the event that the salmon process does not align anything,
        # we still want to succeed, but we want to make our files empty:
        for f in "out_folder/quant.sf" "out_folder/cmd_info.json" "out_folder/lib_format_counts.json" "out_folder/aux_info/ambig_info.tsv" "out_folder/aux_info/eq_classes.txt.gz" "out_folder/aux_info/expected_bias.gz" "out_folder/aux_info/fld.gz" "out_folder/aux_info/meta_info.json" "out_folder/aux_info/observed_bias.gz" "" "out_folder/aux_info/observed_bias_3p.gz" "" "out_folder/logs/salmon_quant.log" ; do
          if [ ! -e $f ] ; then
            echo "Run did not produce output file: $f"
            echo "Creating empty file."
            touch $f
          fi
        done

        # Rename the output files by the fasta name:
        mv out_folder/quant.sf out_folder/~{out_prefix}_quant.sf
        mv out_folder/cmd_info.json out_folder/~{out_prefix}_cmd_info.json
        mv out_folder/lib_format_counts.json out_folder/~{out_prefix}_lib_format_counts.json
        mv out_folder/aux_info/ambig_info.tsv out_folder/aux_info/~{out_prefix}_ambig_info.tsv
        mv out_folder/aux_info/eq_classes.txt.gz out_folder/aux_info/~{out_prefix}_eq_classes.txt.gz
        mv out_folder/aux_info/expected_bias.gz out_folder/aux_info/~{out_prefix}_expected_bias.gz
        mv out_folder/aux_info/fld.gz out_folder/aux_info/~{out_prefix}_fld.gz
        mv out_folder/aux_info/meta_info.json out_folder/aux_info/~{out_prefix}_meta_info.json
        mv out_folder/aux_info/observed_bias.gz out_folder/aux_info/~{out_prefix}_observed_bias.gz
        mv out_folder/aux_info/observed_bias_3p.gz out_folder/aux_info/~{out_prefix}_observed_bias_3p.gz
        mv out_folder/logs/salmon_quant.log out_folder/logs/~{out_prefix}_salmon_quant.log
    >>>

    output {
        File quant_file = "out_folder/~{out_prefix}_quant.sf"
        File cmd_info = "out_folder/~{out_prefix}_cmd_info.json"
        File lib_format_counts = "out_folder/~{out_prefix}_lib_format_counts.json"

        File ambig_info  = "out_folder/aux_info/~{out_prefix}_ambig_info.tsv"
        File eq_classes  = "out_folder/aux_info/~{out_prefix}_eq_classes.txt.gz"
        File expected_bias  = "out_folder/aux_info/~{out_prefix}_expected_bias.gz"
        File fld  = "out_folder/aux_info/~{out_prefix}_fld.gz"
        File meta_info  = "out_folder/aux_info/~{out_prefix}_meta_info.json"
        File observed_bias  = "out_folder/aux_info/~{out_prefix}_observed_bias.gz"
        File observed_bias_3p  = "out_folder/aux_info/~{out_prefix}_observed_bias_3p.gz"

        File log = "out_folder/logs/~{out_prefix}_salmon_quant.log"
    }

    runtime {
        # This should not require much in the way of memory, space, or cpus since we're doing 1 read at a time.
        # Let's try 2 cores. 4gb mem (e2-small machine - $0.01055/hr non-preemptible)
        docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.2"
        memory: 4 + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        boot_disk_gb: 10
        preemptible: 0
        cpu: 2
    }
}

task Alevin {

    meta {
        description : "Run Alevin on a dataset for single cell transcriptome analysis."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File mates1_fastq
        File mates2_fastq
        File cb_whitelist

        File transcript_fasta
        File tgmap

        Int index_kmer_size = 31
        Boolean is_gencode_transcriptome = false

        String library_type = "ISR"
        Int end = 5
        Int umi_length = 10
        Int barcode_length = 16

        String prefix = "alevin_output"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 4*ceil(size(transcript_fasta, "GB")) + 2*ceil(size(mates1_fastq, "GB")) + 2*ceil(size(mates2_fastq, "GB")) + 3*ceil(size(cb_whitelist, "GB"))

    String tgmap_file = if defined(tgmap) then tgmap else ""

    command <<<
        set -euxo pipefail

        # Make a salmon index for the given transcript fasta:
        gencode_arg=""
        if ~{is_gencode_transcriptome} ; then
            gencode_arg="--gencode"
        fi
        salmon index -i salmon_index -k ~{index_kmer_size} $gencode_arg -p 24 -t ~{transcript_fasta}

        # Run alevin:
        salmon alevin \
            --libType ~{library_type} \
            --end ~{end} \
            --umiLength ~{umi_length} \
            --barcodeLength ~{barcode_length} \
            --output ~{prefix}_k~{index_kmer_size} \
            --whitelist ~{cb_whitelist} \
            --mates1 ~{mates1_fastq} \
            --mates2 ~{mates2_fastq} \
            --index salmon_index \
            --tgMap ~{tgmap} \
            --threads 24 \
            --dumpMtx \
            --dumpBarcodeEq \
            --dumpArborescences \
            --dumpfq \
            --dumpUmiGraph \
            --dumpEq \
            > ~{prefix}_k~{index_kmer_size}.fastq

        # Zip up the output folder:
        echo "Zipping up the output..."
        tar -zcf ~{prefix}_k~{index_kmer_size}.tar.gz ~{prefix}_k~{index_kmer_size}
        echo "Done"
    >>>

    output {
        File output_tar_gz = "~{prefix}_k~{index_kmer_size}.tar.gz"
        File output_fastq = "~{prefix}_k~{index_kmer_size}.fastq"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          24,
        mem_gb:             64,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.2"
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
}

task ConvertQuantFileToCountMatrix {

    meta {
        description : "Convert a set of files produced by `salmon quant` (via RunSalmonQuantTask) into a single count matrix containing all the data."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        Array[File] quant_files
    }

    parameter_meta {
        quant_files : "Array of salmon quant.sf files to convert into a single count matrix."
    }

    # Get the directory in which the first file lives:
    File first_file = quant_files[0]

    # 5x the total size of all input files - probably overkill
    Int disk_size = 5*ceil(size(quant_files, "GB"))

    command {
        # Run the script on the current directory, which should have the quant.sf files in it.
        # This will combine all the counts into a single TSV file, then it will convert that TSV to
        # an h5ad file for convenience of importing into scanpy later.

        d=$(dirname ~{first_file})
        echo "Reading quant files from: $d"

        echo "Contents of $d :"
        ls -la $d

        /python_scripts/process_quant_files_into_tsv.py $d
    }

    output {
        File count_matrix_tsv = "differential_expression-known_isoforms.tsv"
        File count_matrix_h5ad = "differential_expression-known_isoforms.h5ad"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.2"
        memory: 16 + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        boot_disk_gb: 10
        preemptible: 0
        cpu: 8
    }
}

