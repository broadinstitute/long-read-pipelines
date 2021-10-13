version 1.0

task SplitBamBySampleAndCellBarcodeTask {

    meta {
        description : "Convert a single annotated (via the 10x tool), aligned bam file into individual FASTA files named by sample name and cell barcode.  Also produces a manifest file for FLAIR to easily quantify output."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File aligned_annotated_bam
        String output_base_name = "reads"
    }

    parameter_meta {
        aligned_annotated_bam : "Bam file containing aligned reads that have been annotated with the 10x tool."
        output_base_name : "[optional] base name to give to every output file.  Should correspond to some unique identifier from this dataset."
    }

    # 10x the total size of the input bam (uncompressed reads)
    # 1x for the file itself
    # 1x for wiggle-room
    # 2x for tar/gz-ing the output:
    Int disk_size = ((10+1+1)*2) * ceil(size(aligned_annotated_bam, "GB"))

    String fasta_tar_gz_name = "fasta_files_by_sample_and_barcode.tar.gz"

    command {
        /python_scripts/split_annotated_reads_by_sample_and_cell_barcode.py -b ~{aligned_annotated_bam} -o ~{output_base_name}
        tar -zcf ~{fasta_tar_gz_name} *.fasta
    }

    output {
        File flair_manifest = "${output_base_name}_flair_reads_manifest.tsv"
        Array[File] sample_cell_barcode_fasta_files = glob("*.fasta")
        File fasta_tar_gz_out = "~{fasta_tar_gz_name}"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.1"
        memory: 16 + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        boot_disk_gb: 10
        preemptible: 0
        cpu: 8
    }
}

task DownsampleToIsoSeqEquivalent {

    meta {
        description : "Downsample a given MAS-seq array element bam file into one containing only 1 read per ZMW (equivalent to IsoSeq)."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File array_element_bam
        String prefix = "downsampled_masseq"
    }

    parameter_meta {
        array_element_bam : "Bam file containing aligned reads that have been annotated with the 10x tool."
        prefix : "[optional] base name to give to every output file.  Should correspond to some unique identifier from this dataset."
    }

    Int disk_size = 10 + 20 * ceil(size(array_element_bam, "GB"))

    String out_name = basename(array_element_bam, ".bam") + ".ZMW_downsampled.bam"

    command {
        /python_scripts/downsample_masseq_by_zmw.py ~{array_element_bam}

        # TODO: THIS IS A HACK - FIX IT LATER
        mv ~{out_name} ~{prefix}.bam
    }

    output {
        File downsampled_bam = "${prefix}.bam"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.6"
        memory: 16 + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        boot_disk_gb: 10
        preemptible: 0
        cpu: 2
    }
}

task DemuxMasSeqDataByIndex {

    meta {
        description : "This workflow will split MAS-seq data that is indexed with a 10bp sequence at the 3' end."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File array_bam
    }

    parameter_meta {
        array_bam : "Bam file containing annotated MAS-seq array reads that contain a 10bp index near the 3' end.."
    }

    Int disk_size = 10 + 20 * ceil(size(array_bam, "GB"))

    String base_out_name = basename(array_bam, ".bam")

    command {
        /python_scripts/mas_seq_demux_by_index.py ~{array_bam} 1> demux_by_index.log
    }

    output {
        File demux_i1 = base_out_name + ".demux_i1.bam"
        File demux_i2 = base_out_name + ".demux_i2.bam"
        File demux_i3 = base_out_name + ".demux_i3.bam"
        File demux_i4 = base_out_name + ".demux_i4.bam"
        File demux_ambiguous = base_out_name + ".demux_ambiguous_indices.bam"
        File log_file = "demux_by_index.log"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.8"
        memory: 4 + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        boot_disk_gb: 10
        preemptible: 0
        cpu: 2
    }
}

task MergeDemuxMasSeqByIndexLogs {

    meta {
        description : "This workflow will merge logs from the DemuxMasSeqDataByIndex task."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        Array[File] demux_logs
    }

    parameter_meta {
        demux_logs : "Log files from DemuxMasSeqDataByIndex task."
    }

    Int disk_size = 10 + 20 * ceil(size(demux_logs, "GB"))
    
    String out_log_name = "merged_demux_log.log"
    
    command <<<

        OUT_LOG="~{out_log_name}"
        
        t="t.txt"
        a="a.txt"
        ap="ap.txt"
        i1="i1.txt"
        i2="i2.txt"
        i3="i3.txt"
        i4="i4.txt"
        tt="tt.txt"
        tpr="tpr.txt"
        trps="trps.txt"

        rm -f $OUT_LOG

        for f in ~{sep=' ' demux_logs} ; do

            grep "^Ambiguous read:" $f >> $OUT_LOG

            # Get the data from the log file:
            tail -n 11 $f | head -n1 | awk '{print $NF}' >> $t
            tail -n 10 $f | head -n1 | awk '{print $NF}' >> $a
            tail -n 9 $f | head -n1 | awk '{print $NF}' | tr -d '%' >> $ap

            tail -n 5 $f | head -n1 | awk '{print $1}' >> $i1
            tail -n 5 $f | head -n1 | awk '{print $2}' >> $i2
            tail -n 5 $f | head -n1 | awk '{print $3}' >> $i3
            tail -n 5 $f | head -n1 | awk '{print $4}' >> $i4

            tail -n 3 $f | head -n1 | awk '{print $NF}' | tr -d 's' >> $tt
            tail -n 2 $f | head -n1 | awk '{print $NF}' | tr -d 's' >> $tpr
            tail -n 1 $f | head -n1 | awk '{print $NF}' | tr -d 's' >> $trps
        done

        awk 'BEGIN{s=0};{s+=$1};END{printf("Total reads seen:     %d\n", s)}' $t >> $OUT_LOG
        awk 'BEGIN{s=0};{s+=$1};END{printf("Num ambiguous reads:  %d\n", s)}' $a >> $OUT_LOG
        awk 'BEGIN{s=0};{s+=$1};END{printf("Ambiguity percentage: %2.3f%%\n", s/NR)}' $ap >> $OUT_LOG
        echo "" >> $OUT_LOG
        echo -e "1\t2\t3\t4" >> $OUT_LOG
        echo "Reads demuxed by index:" >> $OUT_LOG
        awk 'BEGIN{s=0};{s+=$1};END{printf("%d\t", s)}' $i1 >> $OUT_LOG
        awk 'BEGIN{s=0};{s+=$1};END{printf("%d\t", s)}' $i2 >> $OUT_LOG
        awk 'BEGIN{s=0};{s+=$1};END{printf("%d\t", s)}' $i3 >> $OUT_LOG
        awk 'BEGIN{s=0};{s+=$1};END{printf("%d", s)}' $i4 >> $OUT_LOG
        echo "" >> $OUT_LOG
        echo "" >> $OUT_LOG
        awk 'BEGIN{s=0};{s+=$1};END{printf("Time (total elapsed):    %2.3fs\n", s/NR)}' $tt >> $OUT_LOG
        awk 'BEGIN{s=0};{s+=$1};END{printf("Time (per read):         %2.3fs\n", s/NR)}' $tpr >> $OUT_LOG
        awk 'BEGIN{s=0};{s+=$1};END{printf("Time (reads per second): %2.3f\n", s/NR)}' $trps >> $OUT_LOG
        echo "" >> $OUT_LOG

    >>>

    output {
        File merged_log = out_log_name
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.8"
        memory: 4 + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        boot_disk_gb: 10
        preemptible: 0
        cpu: 2
    }
}