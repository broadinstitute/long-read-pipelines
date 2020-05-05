version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/PBUtils.wdl" as PB

workflow QuantifySIRVs {

    meta {
        description : "Quantify SIRV isoforms simple counting of aligned reads to transcripts.  Assumes each read is a single, complete transcript."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        # ------------------------------------------------
        # Input args:
        # Required:

        File reads_bam
        File sirv_tsv
        String prefix = ""

        Int num_redads_per_shard = 200000
    }

    parameter_meta {
        reads_bam : "Bam file containing transcript reads to be aligned to the given transcript isoforms file."
        sirv_tsv : "TSV file containing the SIRVs to quantify and their concentrations  (format: Name|Concentration|PolyA_Length|Sequence)."
        prefix : "[optional] Prefix to prepend to the output file name."

        num_shards : "[optional] Number of shards to parallelize across.  [Default: 20]"
    }
    call Utils.ShardLongReads { input: unmapped_files = [ reads_bam ], num_reads_per_split = num_redads_per_shard}

    scatter (subreads in ShardLongReads.unmapped_shards) {
        call QuantifySIRVsTask {
            input:
                reads_bam = subreads,
                sirv_tsv  = sirv_tsv,
                prefix    = prefix + "_"
        }
    }

    call Utils.MergeBams as MergeAlignedReads { input: bams = QuantifySIRVsTask.aligned_reads, prefix = "~{prefix}_aligned" }
    call Utils.MergeBams as MergeUnmappedReads { input: bams = QuantifySIRVsTask.unmapped_reads, prefix = "~{prefix}_unaligned" }
    call Utils.MergeTsvFiles as MergeCountTables { input: tsv_files = QuantifySIRVsTask.count_table, prefix = "~{prefix}_counts" }

    # ------------------------------------------------
    # Outputs:
    output {
      # Default output file name:
      File aligned_reads        = MergeAlignedReads.merged_bam
      File aligned_reads_index  = MergeAlignedReads.merged_bai
      File unmapped_reads       = MergeUnmappedReads.merged_bam
      File unmapped_reads_index = MergeUnmappedReads.merged_bai
      File count_table          = MergeCountTables.merged_tsv
    }
}


task QuantifySIRVsTask {

    meta {
        description : "Quantify SIRV isoforms simple counting of aligned reads to transcripts.  Assumes each read is a single, complete transcript."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        # ------------------------------------------------
        # Input args:
        # Required:

        File reads_bam
        File sirv_tsv
        String prefix = ""
    }

    parameter_meta {
        reads_bam : "Bam file containing transcript reads to be aligned to the given transcript isoforms file."
        sirv_tsv : "TSV file containing the SIRVs to quantify and their concentrations  (format: Name|Concentration|PolyA_Length|Sequence)."
        prefix : "[optional] Prefix to prepend to the output file name."
    }

    # ------------------------------------------------
    # Process input args:

    String timing_output_file = "timingInformation.txt"
    String memory_log_file = "memory_log.txt"

    String out_base_name="~{prefix}quantified_transcript_reads"

    # ------------------------------------------------
    # Get machine settings:

    Int input_files_size_gb = ceil(10*(size(reads_bam, "GiB") + size(sirv_tsv, "GiB")) + 2)

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

        startTime=`date +%s.%N`
        echo "StartTime: $startTime" > ~{timing_output_file}

        ##########################################################################################
        # Do the real work here:

        # Create a fasta file from our TSV:
        tail -n+2 ~{sirv_tsv} |  awk 'BEGIN{FS="\t"}{printf ">%s\n%s\n", $1, $4}' > SIRVs.fasta

        # Create a python file to run containing our info:
        cat >count_sam_transcripts.py <<'END_SCRIPT'
#!/usr/bin/env python
import sys
import csv

sirv_concentration_dict = dict()
sirv_freq_tsv = None
if len(sys.argv) > 1:
    sirv_freq_tsv = sys.argv[1]
if sirv_freq_tsv:
    with open(sirv_freq_tsv) as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)
        total_concentration = 0
        for row in reader:
            # Only some rows have frequency information:
            if row[1]:
                sirv_concentration_dict[row[0]] = float(row[1])
                total_concentration += float(row[1])
            else:
                sirv_concentration_dict[row[0]] = 0

actual_tx_count = 0
tx_qual_list_dict = dict()
for line in sys.stdin.readlines():
    sp = line.split("\t")

    tx_name = sp[2]
    qual = int(sp[4])

    try:
        tx_qual_list_dict[tx_name].append(qual)
    except KeyError:
        tx_qual_list_dict[tx_name] = [qual]

    actual_tx_count += 1

print("Transcript_Name\tCount\tExpected_Count\tConcentration_(fraction_of_library)\tExpected_Concentration\tMean_Quality")
for tx_name, quals in dict(sorted(tx_qual_list_dict.items(), key=lambda item: len(item[1]), reverse=True)).items():
    expected_concentration_percentage = sirv_concentration_dict[tx_name]/total_concentration
    expected_count = actual_tx_count*expected_concentration_percentage

    actual_concentration_percentage = len(quals)/actual_tx_count

    print(f"{tx_name}\t{len(quals)}\t{expected_count:2.3f}\t{actual_concentration_percentage:2.10f}\t{expected_concentration_percentage:2.10f}\t{(sum(quals)/len(quals)):2.3f}")

END_SCRIPT
        chmod +x count_sam_transcripts.py

        echo "Indexing input transcript isoforms fasta..."
        samtools faidx SIRVs.fasta

        echo "Aligning reads and filtering for only primary alignments..."
        samtools fastq ~{reads_bam} \
            | minimap2 -ayYL --MD --eqx -x asm20 -t4 SIRVs.fasta - \
            | samtools view -hb -F 2304 - \
        > ~{out_base_name}.bam

        echo "Getting unmapped reads..."
        samtools view -hb -f4 ~{out_base_name}.bam > ~{out_base_name}.unmapped.bam

        echo "Quantifying alignments..."
        samtools view -F 2308 ~{out_base_name}.bam | ./count_sam_transcripts.py ~{sirv_tsv} > ~{out_base_name}transcript_count_table.tsv


        ##########################################################################################

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
        cpu:                    2
        memory:                 4 + " GiB"
        disks:                  "local-disk " +  input_files_size_gb + " HDD"
        bootDiskSizeGb:         10
        preemptible:            0
        maxRetries:             1
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
    }

    # ------------------------------------------------
    # Outputs:
    output {
      # Default output file name:
      File aligned_reads  = "~{out_base_name}.bam"
      File unmapped_reads = "~{out_base_name}.unmapped.bam"
      File count_table    = "~{out_base_name}transcript_count_table.tsv"

      File timing_info    = "${timing_output_file}"
      File memory_log     = "${memory_log_file}"
    }
}