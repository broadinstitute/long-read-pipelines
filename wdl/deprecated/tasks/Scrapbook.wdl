version 1.0

import "../../structs/Structs.wdl"

task CoverageTrack {
    input {
        File bam
        File bai
        String chr
        String start
        String end

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB"))

    command <<<
        set -euxo pipefail

        samtools depth -a ~{bam} -r ~{chr}:~{start}-~{end} | bgzip > ~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz
        tabix -p bed ~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz
    >>>

    output {
        File coverage = "~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz"
        File coverage_tbi = "~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
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

task FilterMQ0Reads {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB"))
    String prefix = basename(bam, ".bam")

    command <<<
        set -euxo pipefail

        samtools view -q 1 -b ~{bam} > ~{prefix}.no_mq0.bam
        samtools index ~{prefix}.no_mq0.bam
    >>>

    output {
        File no_mq0_bam = "~{prefix}.no_mq0.bam"
        File no_mq0_bai = "~{prefix}.no_mq0.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
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

task ComputeBedCoverage {
    input {
        File bam
        File bai
        File bed
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB") + size(bed, "GB"))

    command <<<
        set -euxo pipefail

        bedtools coverage -b ~{bed} -a ~{bam} -nobuf | gzip > ~{prefix}.txt.gz
        zcat ~{prefix}.txt.gz | awk '{ sum += sprintf("%f", $15*$16) } END { printf("%f\n", sum) }' > ~{prefix}.count.txt
    >>>

    output {
        File coverage = "~{prefix}.txt.gz"
        Float counts = read_float("~{prefix}.count.txt")
        File counts_file = "~{prefix}.count.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
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

task BamToBed {
    input {
        File bam
        File bai

        RuntimeAttr? runtime_attr_override
    }

    String bed = basename(bam, ".bam") + ".bed"
    Int disk_size = 4*ceil(size(bam, "GB") + size(bai, "GB"))

    command <<<
        set -euxo pipefail

        bedtools bamtobed -i ~{bam} > ~{bed}
    >>>

    output {
        File bedfile = bed
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
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

task NanoPlotFromUBam {

    meta {
        description: "NanoPlot from an unaligned BAM"
    }

    parameter_meta {
        bam: "BAM file"
        runtime_attr_override: "Runtime attributes to override"
    }

    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        NanoPlot -t ${num_core} \
            -c orangered \
            --N50 \
            --tsv_stats \
            --ubam "~{bam}"

        cat NanoStats.txt | \
            grep -v -e '^Metrics' -e '^highest' -e '^longest' | \
            sed 's/ >/_/' | \
            sed 's/://' | \
            awk '{ print $1 "\t" $2 }' | \
            tee map.txt
    >>>

    #number_of_reads 991
    #number_of_bases 12949457.0
    #median_read_length      13705.0
    #mean_read_length        13067.1
    #read_length_stdev       9581.3
    #n50     18618.0
    #mean_qual       0.0
    #median_qual     0.0
    #Reads_Q5        0
    #Reads_Q7        0
    #Reads_Q10       0
    #Reads_Q12       0
    #Reads_Q15       0

    output {
        File stats = "NanoStats.txt"
        Map[String, Float] stats_map = read_map("map.txt")

        Array[File] plots = glob("*.png")
        File Non_weightedHistogramReadlength = "Non_weightedHistogramReadlength.png"
        File Non_weightedLogTransformed_HistogramReadlength = "Non_weightedLogTransformed_HistogramReadlength.png"
        File WeightedHistogramReadlength = "WeightedHistogramReadlength.png"
        File WeightedLogTransformed_HistogramReadlength = "WeightedLogTransformed_HistogramReadlength.png"
        File Yield_By_Length = "Yield_By_Length.png"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/biocontainers/nanoplot:1.35.5--pyhdfd78af_0"
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

task SAMtoPAF {
    input {
        File sam_formatted_file
        File? index

        RuntimeAttr? runtime_attr_override
    }
    meta {
        description: "Convert SAM-formatted alignment to PAF format using minimap2's paftools.js"
    }
    parameter_meta {
        sam_formatted_file: "SAM-formated input file to be converted to PAF (note currently we only support SAM or BAM, not CRAM)"
        index:              "[optional] index for sam_formatted_file"
    }

    String prefix = basename(basename(sam_formatted_file, ".bam"), ".sam") # we have hack like this because WDL stdlib doesn't provide endsWith stuff

    Int disk_size = 2*ceil(size(sam_formatted_file, "GB"))

    command <<<
        set -eu

        MM2_VERSION="2.24"

        filename=$(basename -- ~{sam_formatted_file})
        extension="${filename##*.}"
        if [[ "$extension" == "sam" ]]; then
            /minimap2-${MM2_VERSION}_x64-linux/k8 \
                /minimap2-${MM2_VERSION}_x64-linux/paftools.js \
                sam2paf \
                -L \
                ~{sam_formatted_file} \
                > ~{prefix}".paf"
        elif [[ "$extension" == "bam" ]]; then
            samtools view -h ~{sam_formatted_file} | \
                /minimap2-${MM2_VERSION}_x64-linux/k8 \
                /minimap2-${MM2_VERSION}_x64-linux/paftools.js \
                sam2paf \
                -L \
                - \
                > ~{prefix}".paf"
        else
            echo "Currently we only support SAM or BAM (not CRAM)." && exit 1;
        fi
    >>>

    output {
        File pat_formatted_file = "~{prefix}.paf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
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

task RaconPolish {

    meta {
        description: "Polish a draft assembly with long reads using Racon. Recommended to run a few times."
    }
    parameter_meta {
        reads:          "long reads to polish the draft assembly with"
        draft_assembly: "draft to be polished"
        n_rounds: "Number of times to run Racon"
    }

    input {
        File reads
        File draft_assembly

        Int n_rounds
    }

    Int mem_size = 4 * ceil(size(reads, "GB") + size(draft_assembly, "GB"))
    Int disk_size = mem_size

    command <<<
        set -euxo pipefail

        cp ~{draft_assembly} input_draft.fasta
        for i in {1..~{n_rounds}}
        do
          minimap2 -ax map-ont input_draft.fasta ~{reads} > aln.sam
          racon -c 1 -m 8 -x -6 -g -8 -w 500 -t 8 ~{reads} aln.sam input_draft.fasta > polished_${i}_draft.fasta
          cp polished_${i}_draft.fasta input_draft.fasta
        done
    >>>

    output {
        File final_polished_assembly = "input_draft.fasta"
        Array[File] incremental_polished_assemblies = glob("polished_*_draft.fasta")
    }

    runtime {
        cpu:                    8
        # Racon has a high memory requirement. Not sure what it is exactly but you need at least
        # the size of the generated alignment file and more
        memory:                 mem_size + " GiB"
        disks:                  "local-disk " + disk_size + " HDD"
        bootDiskSizeGb:         30
        preemptible:            0
        maxRetries:             0
        gpuType:                "nvidia-tesla-t4"
        gpuCount:               1
        nvidiaDriverVersion:    "418.152.00"
        zones:                  ["us-east1-c"]
        cpuPlatform:            "Intel Haswell"
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-racon:0.1.0"
    }

}

task CompressAndIndex {
    meta {
        description: "Convert a BCF file to a vcf.bgz file and index it."
    }

    input {
        File joint_bcf

        Int num_cpus = 8
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 3*ceil(size(joint_bcf, "GB"))

    command <<<
        set -x

        bcftools view ~{joint_bcf} | bgzip -@ ~{num_cpus} -c > ~{prefix}.g.vcf.bgz
        tabix -p vcf ~{prefix}.g.vcf.bgz
    >>>

    output {
        File joint_gvcf = "~{prefix}.g.vcf.bgz"
        File joint_gvcf_tbi = "~{prefix}.g.vcf.bgz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             4*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Run_Group {

    meta {
        description : "Run umi-tools group on a bam file."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        aligned_transcriptome_reads: "Aligned reads in bam format."
        aligned_transcriptome_reads_index: "Index for aligned reads in bam format."
        gene_tag: "Tag for gene name in bam file."
        cell_barcode_tag: "Tag for cell barcode in bam file."
        umi_tag: "Tag for UMI in bam file."
        do_per_cell: "Whether to do per-cell grouping."
        prefix: "Prefix for output files."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        File aligned_transcriptome_reads
        File aligned_transcriptome_reads_index

        String gene_tag = "XG"
        String cell_barcode_tag = "CB"
        String umi_tag = "ZU"

        Boolean do_per_cell = true

        String prefix = "umi_tools_group"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 20*ceil(size(aligned_transcriptome_reads, "GB") + size(aligned_transcriptome_reads_index, "GB"))

    String per_cell_args = if do_per_cell then " --per-cell --cell-tag " + cell_barcode_tag + " " else ""

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

        set -euxo pipefail

        # Run umi-tools group:
        umi_tools group \
          --buffer-whole-contig \
          --no-sort-output \
          --per-gene \
          ~{per_cell_args} \
          --gene-tag ~{gene_tag} \
          --extract-umi-method tag \
          --umi-tag ~{umi_tag} \
          -I ~{aligned_transcriptome_reads} \
          --group-out=~{prefix}.tsv \
          --output-bam \
          --log=~{prefix}.log > ~{prefix}.bam


        # Stop the memory daemon softly.  Then stop it hard if it's not cooperating:
        set +e
        DO_MEMORY_LOG=false
        sleep $(($MEM_LOG_INTERVAL_s  * 2))
        kill -0 $mem_pid &> /dev/null
        if [ $? -ne 0 ] ; then
            kill -9 $mem_pid
        fi
    >>>

    output {
        File output_bam = "~{prefix}.bam"
        File output_tsv = "~{prefix}.tsv"
        File memory_log = "~{memory_log_file}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             64,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
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

task SplitBamBySampleAndCellBarcodeTask {

    meta {
        description : "Convert a single annotated (via the 10x tool), aligned bam file into individual FASTA files named by sample name and cell barcode.  Also produces a manifest file for FLAIR to easily quantify output."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        aligned_annotated_bam : "Bam file containing aligned reads that have been annotated with the 10x tool."
        output_base_name : "[optional] base name to give to every output file.  Should correspond to some unique identifier from this dataset."
    }

    input {
        File aligned_annotated_bam
        String output_base_name = "reads"
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
        memory: 16 + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        boot_disk_gb: 25
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

    parameter_meta {
        array_element_bam : "Bam file containing aligned reads that have been annotated with the 10x tool."
        prefix : "[optional] base name to give to every output file.  Should correspond to some unique identifier from this dataset."
    }

    input {
        File array_element_bam
        String prefix = "downsampled_masseq"
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
        memory: 32 + " GiB"  # Need a lot of ram here because we keep a set of ZMWs in memory
        disks: "local-disk " + disk_size + " HDD"
        boot_disk_gb: 25
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
        echo "Reads demuxed by index:" >> $OUT_LOG
        echo -e "1\t2\t3\t4" >> $OUT_LOG
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
        memory: 4 + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        boot_disk_gb: 25
        preemptible: 0
        cpu: 2
    }
}

task SplitBamByContig {
    meta {
        description : "Split a given bam file into separate files by contig."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        bam : "Bamfile to be split by contig."
        prefix : "Prefix for ouput files."
    }

    input {
        File bam

        String prefix = "reads"

        RuntimeAttr? runtime_attr_override
    }


    Int disk_size_gb = 10 + 5*ceil(size(bam, "GB"))

    String contig_list = "contig_list.txt"

    command <<<
        /python_scripts/split_reads_by_contig.py ~{bam} ~{prefix}

        ls ~{prefix}.*.bam | sed 's#~{prefix}.\(.*\).bam#\1#' > ~{contig_list}
    >>>

    output {
        Array[File] contig_bams = glob("~{prefix}*.bam")
        Array[String] contig_names = read_lines(contig_list)
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
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

task AggregateUmiAdjustmentStats {

    meta {
        description : "Aggregates the UMI adjustment stats from the given longbow UMI adjustment log files."
    }

    parameter_meta {
        longbow_umi_adjustment_log_files: "The longbow UMI adjustment log files to aggregate."
        out_name: "The name of the output file."
    }

    # TODO: FINISHME
    input {
        Array[File] longbow_umi_adjustment_log_files

        String out_name = "longbow_umi_adjustment_stats.txt"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(longbow_umi_adjustment_log_files, "GB"))

    # YES, this SHOULD be a proper tool, but right now it isn't.
    command <<<

        f=~{write_lines(longbow_umi_adjustment_log_files)}

        mv $f THE_UMI_FILE_LIST.txt

python << CODE
import os

stats_dict = dict()
line_key = "STATS: "

with open("THE_UMI_FILE_LIST.txt", 'r') as umi_file_list_file:
    for line in umi_file_list_file:
        stats_file = line.strip()
        with open(stats_file, 'r') as f:
            for line in f:
                if line_key in line:
                    line = line.strip()
                    s = line[line.find(line_key) + len(line_key):]
                    key, remainder = [t.strip() for t in s.split(":")]
                    if "/" in remainder:
                        count = int(remainder[:remainder.find("/")])
                        tot = int(remainder[remainder.find("/")+1:remainder.find(" ")])
                    else:
                        count = int(remainder)
                        tot = None

                    try:
                        c, t = stats_dict[key]
                        if tot is not None:
                            tot += t
                        stats_dict[key] = (count + c, tot)
                    except KeyError:
                        stats_dict[key] = (count, tot)

k_len = 0
for k in stats_dict.keys():
    if len(k) > k_len:
        k_len = len(k)

k_prefix = list(stats_dict.keys())[0]
k_prefix = k_prefix[:k_prefix.find(" ")]
with open("~{out_name}", 'w') as f:
    for k, v in stats_dict.items():

        if not k.startswith(k_prefix):
            f.write("\n")
            k_prefix = k[:k.find(" ")]

        k_spacing = k_len - len(k)

        count, tot = v
        if tot is None:
            f.write(f"{k}:{' '*k_spacing} {count}\n")
        else:
            f.write(f"{k}:{' '*k_spacing} {count}/{tot} ({100.0*count/tot:2.4f}%)\n")

CODE
    >>>

    output {
        File stats = out_name
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,             # Decent amount of CPU and Memory because network transfer speed is proportional to VM "power"
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,             # This shouldn't take very long, but it's nice to have things done quickly, so no preemption here.
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.27"
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

task MergeBarcodeCounts {
    meta {
        description : "Merge all counts for each unique barcode in the given TSV file.  Assumes file is unheadered and have two columns: BARCODE COUNT.  Merging performed by adding all COUNTs for each BARCODE."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        barcode_count_tsv: "The TSV file containing the barcode counts to merge."
        prefix: "The prefix to use for the output file."
    }

    input {
        File barcode_count_tsv
        String prefix = "merged_counts"

        RuntimeAttr? runtime_attr_override
    }

    # 20 gb - baseline storage for safety
    # 1x for the file itself
    # 2x for the results and some wiggle room.
    Int disk_size_gb = 20 + (3 * ceil(size(barcode_count_tsv, "GB")))

    command {
        /python_scripts/merge_barcode_counts.py ~{barcode_count_tsv}
        if [[ "~{prefix}.tsv" != "merged_counts.tsv" ]] ; then
            mv merged_counts.tsv "~{prefix}.tsv"
        fi
    }

    output {
        File merged_counts = "~{prefix}.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
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

task CreateCountMatrixAnndataFromTsv {

    meta {
        description : "Creates a python anndata object from the given countmatrix tsv.  Expects the input to have been generated by CreateCountMatrixFromAnnotatedBam.  The resulting anndata object can be directly read into scanpy for single-cell analysis."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        count_matrix_tsv: "The TSV file containing the count matrix to convert to anndata."
        genome_annotation_gtf_file: "The GTF file containing the genome annotation."
        force_anndata_gencode_overwrite: "Forces the values in `gene_names` and `gene_ids` to overwrite `gencode_overlap_gene_names` and `gencode_overlap_gene_ids` respectively."
        prefix: "The prefix to use for the output file."
        equivalence_class_definitions: "The equivalence class definitions file to use for the anndata object."
        overlap_intervals: "The overlap intervals file to use for the anndata object."
        overlap_interval_label: "The overlap interval label to use for the anndata object."
        gencode_reference_gtf_file: "The gencode reference GTF file to use for the anndata object."
    }

    input {
        File count_matrix_tsv
        File genome_annotation_gtf_file

        Boolean force_anndata_gencode_overwrite = false

        String prefix = "umi_tools_group"

        File? equivalence_class_definitions

        File? overlap_intervals
        String? overlap_interval_label
        File? gencode_reference_gtf_file

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 20 + 4*ceil(size(count_matrix_tsv, "GB")) + 4*ceil(size(genome_annotation_gtf_file, "GB")) + 2*ceil(size(equivalence_class_definitions, "GB"))

    String overlap_intervals_arg = if defined(overlap_intervals)  then " --overlap-intervals " else ""
    String overlap_interval_label_arg = if defined(overlap_interval_label) then " --overlap-interval-label " else ""
    String gencode_reference_gtf_file_arg = if defined(gencode_reference_gtf_file) then " --gencode-reference-gtf " else ""

    String force_gencode_overwrite_flag = if force_anndata_gencode_overwrite then " --force-overwrite-gencode-overlaps " else ""

    String eq_class_arg = if defined(equivalence_class_definitions)  then " --eq-class-defs-tsv " else ""

    command <<<
        set -euxo pipefail
        /python_scripts/create_count_matrix_anndata_from_tsv.py \
            -t ~{count_matrix_tsv} \
            -g ~{genome_annotation_gtf_file} \
            ~{overlap_intervals_arg}~{default="" overlap_intervals} \
            ~{overlap_interval_label_arg}~{default="" overlap_interval_label} \
            ~{gencode_reference_gtf_file_arg}~{default="" gencode_reference_gtf_file} \
            ~{eq_class_arg} ~{default="" equivalence_class_definitions} \
            ~{force_gencode_overwrite_flag} \
            -o ~{prefix}
    >>>

    output {
        File transcript_gene_count_anndata_h5ad = "~{prefix}_tx_gene_counts_adata.h5ad"
        Array[File] pickles = glob("*.pickle")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             32,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
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

task SubsetCountsMatrixByGenes {

    meta {
        description : "Subsets a count matrix TSV file to contain only the transcripts from the given list of genes.  Assumes the count matrix was produced by comparison with Gencode (due to data formatting) and that the table is a TSV with samples as rows and transcripts as columns."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        count_matrix_tsv : "TSV file containing the counts of each transcript expressed in a sample.  (Transcripts in columns.  Samples in rows.  One header row.)"
        gene_names : "Array of gene names for which to keep data from the given count matrix TSV."
    }

    input {
        File count_matrix_tsv
        Array[String] gene_names
    }


    # We're subsetting the file, so we should be able to get away with very little space here.
    # 1x for the file itself
    # 2x for the results and some wiggle room.
    Int disk_size = 3 * ceil(size(count_matrix_tsv, "GB"))

    command {
        /python_scripts/subset_count_matrix_by_gene.py ~{count_matrix_tsv} ~{sep=' ' gene_names}
    }

    output {
        File subset_count_matrix_tsv = "count_matrix_subset_by_gene.tsv"
        File subset_count_matrix_h5ad = "count_matrix_subset_by_gene.h5ad"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
        memory: 16 + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        boot_disk_gb: 25
        preemptible: 0
        cpu: 8
    }
}

task GetDefaultDir {

    meta {
        description: "Get the default directory for a given workflow"
    }

    parameter_meta {
        workflow_name: "The name of the workflow"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        String workflow_name

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        NAME=$(cat gcs_localization.sh | grep 'source bucket' | sed 's/# Localize files from source bucket //' | sed 's/ to container.*//' | sed "s/'//g")

        echo "gs://$NAME/results/~{workflow_name}"
    >>>

    output {
        String path = read_string(stdout())
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task PrepareManifest {

    meta {
        description: "Prepare a manifest file for a given workflow"
    }

    parameter_meta {
        files: "The files to include in the manifest"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        Array[String] files

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        echo ~{sep=' ' files} | sed 's/ /\n/g' > manifest.txt
    >>>

    output {
        File manifest = "manifest.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task EchoManifest {

    meta {
        description: "Echo the contents of a manifest file"
    }

    parameter_meta {
        manifest: "The manifest file to echo"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File manifest

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        cat ~{manifest}
    >>>

    output {
        File out = stdout()
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task SortBam {

    meta {
        description: "Sort a BAM file"
    }

    parameter_meta {
        input_bam: "The BAM file to sort"
        prefix: "The prefix for the output BAM file"
        runtime_attr_override: "Override the default runtime attributes"
    }

    input {
        File input_bam
        String prefix = "sorted"

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        samtools sort -@$num_core -o ~{prefix}.bam ~{input_bam}
        samtools index ~{prefix}.bam
    >>>

    output {
        File sorted_bam = "~{prefix}.bam"
        File sorted_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            10,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task FilterListOfStrings {

    meta {
        description : "Filter a given list of files by a query term (essentially pipes the query into grep)."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        list_to_filter: "Array of strings to filter by the query."
        query: "Term to use to filter the given strings."
    }

    input {
        Array[String] list_to_filter
        String query

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        \grep "~{query}" ~{write_lines(list_to_filter)} > filtered_list.txt
    >>>

    output {
        Array[String] filtered_list = read_lines("filtered_list.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "ubuntu:hirsute-20210825"
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

task FilterReadsBySamFlags {

    meta {
        description : "Filter reads based on sam flags.  Reads with ANY of the given flags will be removed from the given dataset."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        bam:   "BAM file to be filtered."
        sam_flags: "Flags for which to remove reads.  Reads with ANY of the given flags will be removed from the given dataset."
        prefix : "[Optional] Prefix string to name the output file (Default: filtered_reads)."
    }

    input {
        File bam
        String sam_flags

        String extra_args = ""

        String prefix = "filtered_reads"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + ceil(11 * size(bam, "GiB"))

    command <<<

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        samtools view -h -b -F ~{sam_flags} -@$np ~{extra_args} ~{bam} > ~{prefix}.bam
        samtools index -@$np ~{prefix}.bam
    >>>

    output {
        File output_bam = "~{prefix}.bam"
        File output_bam_index = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
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

task GrepCountUniqueBamRecords {

    meta {
        description : "Count the number of unique records in a bam file that match a given regex."
    }

    parameter_meta {
        bam: "BAM file to be filtered."
        samfilter: "[Optional] Extra arguments to pass into samtools view."
        regex: "Regex to match against the bam file."
        invert: "[Optional] Invert the regex match."
        prefix: "[Optional] Prefix string to name the output file (Default: sum)."
    }

    input {
        String bam
        String samfilter = ""
        String regex
        Boolean invert = false
        String prefix = "sum"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(2 * size(bam, "GiB"))
    String arg = if invert then "-v" else ""

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view ~{samfilter} ~{bam} | grep ~{arg} ~{regex} | > ~{prefix}.txt
    >>>

    output {
        Int num_records = read_int("~{prefix}.txt")
        File num_records_file = "~{prefix}.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
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

task FilterReadsWithTagValues {

    meta {
        description : "Filter reads from a BAM file based on the values of a given tag."
    }

    parameter_meta {
        bam:   "Input BAM file from which to remove a tag with certain values."
        tag:   "Name of the tag to target for potential removal."
        value_to_remove:   "Tag value to use to remove reads.  Reads will be removed if they have the given tag with this value."
        prefix: "[default-valued] prefix for output BAM"
    }

    input {
        File bam
        String tag
        String value_to_remove
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20 + 11*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        java -jar /usr/picard/picard.jar \
            FilterSamReads \
                --VALIDATION_STRINGENCY SILENT \
                --FILTER excludeTagValues \
                --TAG ~{tag} \
                --TAG_VALUE ~{value_to_remove} \
                -I ~{bam} \
                -O ~{prefix}.bam
    >>>

    output {
        File output_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "broadinstitute/picard:2.23.7"
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

# A utility to subset a BAM to specifed loci
task ExcludeRegionsFromBam {

    meta {
        description : "Exclude regions from a BAM file."
    }

    parameter_meta {
        bam:   "Input BAM file to subset."
        bai:   "Index for input BAM file."
        loci: "Genomic loci to exclude."
        prefix: "[default-valued] prefix for output BAM and BAI file names."
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        File bai
        Array[String] loci
        String prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        echo ~{sep=',' loci} | sed 's/,/\n/g' | sed 's/[:-]/\t/g' > regions.bed
        samtools view -L regions.bed -hbU ~{prefix}.bam -o /dev/null ~{bam}
        samtools index ~{prefix}.bam
    >>>

    output {
        File subset_bam = "~{prefix}.bam"
        File subset_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

# A utility to select the first N reads from a BAM file
task SelectFirstNReads {

    meta {
        description : "Select the first N reads from a BAM file."
    }

    parameter_meta {
        bam: {
                 description: "bam to subset",
                 localization_optional: true
             }
        n: "number of reads to select"
        prefix: "prefix for output bam file"
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        Int n
        String prefix = "selected"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(bam, "GB"))

    command <<<
        set -x

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        ((samtools view -H ~{bam}) && (samtools view ~{bam} | head -n ~{n})) | samtools view -b > ~{prefix}.bam
    >>>

    output {
        File selected_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

task SplitBam {

    meta {
        description: "Takes a BAM file as input, extracts the reference sequence names from the header, filters out any unwanted sequences, and creates new BAM files for each remaining sequence."
    }

    parameter_meta {
        bam:    "bam to split"
        bai:    "index for bam file"
        filter: "contigs to ignore"
    }

    input {
        File bam
        File bai
        Array[String] filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA']

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size([bam, bai], "GB"))

    command <<<
        set -euxo pipefail

        samtools view -H ~{bam} | \
            grep '^@SQ' | \
            grep -v -e '^@HD' ~{true='-e' false='' length(filter) > 0} ~{sep=" -e " filter} | \
            awk '{ print $2 }' | \
            sed 's/SN://' |
            parallel -j+0 "samtools view -bh -o {}.bam ~{bam} {} && samtools index {}.bam"
    >>>

    output {
        Array[File] subset_bams = glob("*.bam")
        Array[File] subset_bais = glob("*.bai")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task FilterBamOnTag {

    meta {
        description: "Filters a BAM file on a given tag and expression"
    }

    parameter_meta {
        bam:        "input BAM file"
        prefix:     "prefix for output bam"
        tag:        "tag to filter on"
        expression: "expression to filter on"
    }

    input {
        File bam
        String prefix = "out"
        String tag
        String expression

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        bamtools filter -in ~{bam} -out ~{prefix}.bam -tag "~{tag}":"~{expression}"
        samtools index ~{prefix}.bam
    >>>

    output {
        File filtered_bam = "~{prefix}.bam"
        File filtered_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
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

task ListBamContigs {

    meta {
        description: "Utility to list contigs in a BAM file"
    }

    parameter_meta {
        bam:    "input BAM from which available contigs should be listed"
    }

    input {
        String bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} | grep '^@SQ' | awk '{ print $2 }' | sed 's/SN://' > chrs.txt
    >>>

    output {
        Array[String] contigs = read_lines("chrs.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task ShardReads {

    meta {
        description: "Shard a bam file by number of reads"
    }

    parameter_meta {
        bam: "bam file to shard"
        bam_index: "bam index file"
        prefix: "prefix for output files"
        num_shards: "number of shards to create"
        runtime_attr_override: "override the runtime attributes"
    }

    input {
        File bam
        File bam_index

        String prefix = "shard"

        Int num_shards = 10

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(4 * size(bam, "GiB"))

    String sharded_bam_folder = "sharded_bams"

    command <<<
        num_reads=$(samtools idxstats ~{bam} | awk 'BEGIN{s=0}{s+=$3;s+=$4}END{print s}')

        mkdir sharded_bams

        java -jar /usr/picard/picard.jar \
            SplitSamByNumberOfReads \
                --VALIDATION_STRINGENCY SILENT \
                -I ~{bam} \
                -O ~{sharded_bam_folder} \
                -OUT_PREFIX ~{prefix} \
                -N_FILES ~{num_shards} \
                -TOTAL_READS ${num_reads}
    >>>

    output {
        Array[File] shards = glob("~{sharded_bam_folder}/*.bam")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9.gamma"
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

task FailWithWarning {

    meta {
        description: "Fail with a warning"
    }

    parameter_meta {
        warning: "The warning message to print"
    }

    input {
        String warning
    }
    command <<<
        set -e

        echo "~{warning}"
        echo "~{warning}" 1>&2
        false
    >>>
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            10,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
    runtime {
        cpu:                    default_attr.cpu_cores
        memory:                 select_first([default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         default_attr.boot_disk_gb
        preemptible:            default_attr.preemptible_tries
        maxRetries:             default_attr.max_retries
        docker:                 default_attr.docker
    }
}

task SplitDelimitedString {

    meta {
        description: "Split a delimited string into an array"
    }

    parameter_meta {
        s: "The string to split"
        separate: "The delimiter to split on"
    }

    input {
        String s
        String separate
    }

    command <<<
        set -eux

        echo ~{s} | tr ~{separate} '\n' > result.txt
    >>>

    output {
        Array[String] arr = read_lines("result.txt")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task IndexVCF {

    meta {
        description: "Indexing vcf.gz. Note: do NOT use remote index as that's buggy."
    }

    parameter_meta {
        vcf: "VCF file to be indexed"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File vcf
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(vcf, ".vcf.gz")
    Int proposed_disk = 3*ceil(size(vcf, "GB")) + 1
    Int disk_size = if (proposed_disk > 100) then proposed_disk else 100

    command <<<
        cp ~{vcf} ~{prefix}.vcf.gz && \
            tabix -p vcf ~{prefix}.vcf.gz && \
            find ./ -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g'
    >>>

    output {
        File tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             3,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
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

task FindSequencingSummaryFiles {

    meta {
        description: "Find sequencing summary files in an ONT basecall directory."
    }

    parameter_meta {
        gcs_input_dir: "GCS directory containing sequencing summary files."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        String gcs_input_dir

        RuntimeAttr? runtime_attr_override
    }

    String indir = sub(gcs_input_dir, "/$", "")

    command <<<
        for summary_file in $(gsutil ls "~{indir}/**sequencing_summary*.txt*")
        do
            DIR=$(dirname $summary_file)
            echo ${DIR}

            gsutil ls "${DIR}" | grep fastq_pass && gsutil ls "${DIR}" | grep fast5_pass

            if [ $? -eq 0 ]; then
                FASTQ_COUNT=$(gsutil ls "${DIR}/fastq_pass/*.fastq*" | wc -l)
                FAST5_COUNT=$(gsutil ls "${DIR}/fast5_pass/*.fast5*" | wc -l)

                echo "${FASTQ_COUNT} ${FAST5_COUNT}"

                if [ ${FASTQ_COUNT} -eq ${FAST5_COUNT} ]; then
                    echo $summary_file >> summaries.txt
                else
                    echo "# fastq != # fast5.  Skipped."
                fi
            else
                echo "No passing fastq and fast5 files.  Skipped."
            fi

            echo ""
        done
    >>>

    output {
        Array[String] summary_files = read_lines("summaries.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task MASIsoSeqReport {

    meta {
        description : "Create a report for a given MAS-ISO-Seq run which summarizes the results using a given Jupyter Notebook template."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        notebook_template : "Jupyter notebook MASSeq template to run with the given data to produce a MASSeq report."

        sample_name : "Name of the MAS-seq sample being analyzed in this report."

        subreads_stats : "Samtools stats file created from the raw subreads from the PacBio instrument."
        ccs_reads_stats : "Samtools raw stats file created from the aligned CCS corrected reads from the PacBio instrument."
        array_elements_stats : "Samtools raw stats file created from the aligned MASSeq array elements."

        ccs_report_file : "CCS report file from the CCS run for the data from the PacBio instrument."
        raw_ccs_bam_file : "Unaligned reads file in BAM format from the CCS process (pre-array splitting)."

        array_element_bam_file : "Transcriptome aligned reads file in BAM format containing aligned MASSeq array elements as individual reads."
        array_elements_genome_aligned : "Genome aligned reads file in BAM format containing aligned MASSeq array elements as individual reads."
        ccs_rejected_bam_file : "Bam file containing all subreads from zmws that were rejected by CCS."

        annotated_bam_file : "Bam file containing ccs corrected reads with annotated sections in the SG tag."

        longbow_passed_reads_file : "Bam file containing all reads that passed the longbow filter for the model used in this run (both ccs passed and reclaimed)."
        longbow_failed_reads_file : "Bam file containing alll reads that failed the longbow filter for the model used in this run (both ccs passed and reclaimed)."

        longbow_passed_ccs_reads : "Bam file containing ccs corrected reads that passed the longbow filter for the model used in this run (CCS Corrected reads ONLY)."
        longbow_failed_ccs_reads : "Bam file containing ccs corrected reads that failed the longbow filter for the model used in this run (CCS Corrected reads ONLY)."
        ccs_reclaimable_reads : "Bam file containing ccs rejected reads that are deemed to be reclaimable."
        ccs_reclaimed_reads : "Bam file containing ccs rejected reads that have been reclaimed."
        ccs_rejected_longbow_failed_reads : "Bam file containing ccs reclaimable reads that did not pass longbow filtering and were not reclaimed."
        raw_array_elements : "Bam file containing the raw unaligned array elements created from the longbow_passed_reads_file."
        ccs_reclaimed_array_elements : "Bam file containing the unaligned array elements created from reclaimed CCS reads."

        zmw_stats_json_gz : "ZMW stats json.gz file from the PacBio instrument."

        zmw_subread_stats_file : "[optional] File containing statistics about the subreads from each ZMW (created by collect_zmw_subread_stats.py in the PBUtils docker container)."
        polymerase_read_lengths_file : "[optional] File containing the lengths of each polymerase read from the sequencer (as created by collect_polymerase_read_lengths.py)"
        approx_raw_subread_array_lengths : "[optional] File containing the approximate array length information from the raw (pre-ccs) subreads file  (created by get_approx_raw_subread_array_lengths.py in the Cartographer docker container)."

        ten_x_metrics_file : "[optional] Stats file from the 10x tool run for the data in this MASSeq run.  If not supplied stats will not be displayed in the resulting report."
        mas_seq_model : "Built-in mas-seq model to use."

        workflow_dot_file : "DOT file containing the representation of this workflow used to create and analyze the data.  This is included in the QC reports (the DOT file can be generated with womtool)."

        prefix : "[optional] Prefix to prepend to the name of the generated report."
        runtime_attr_override : "[optional] Runtime attributes struct with which to override the docker container runtime.."
    }

    input {
        File notebook_template

        String sample_name

        File subreads_stats
        File ccs_reads_stats
        File array_elements_stats
        File ccs_report_file

        File raw_ccs_bam_file
        File array_element_bam_file
        File array_elements_genome_aligned
        File ccs_rejected_bam_file

        File annotated_bam_file

        File longbow_passed_reads_file
        File longbow_failed_reads_file

        File longbow_passed_ccs_reads
        File longbow_failed_ccs_reads
        File ccs_reclaimable_reads
        File ccs_reclaimed_reads
        File ccs_rejected_longbow_failed_reads
        File raw_array_elements
        File ccs_reclaimed_array_elements

        File zmw_stats_json_gz

        File? zmw_subread_stats_file
        File? polymerase_read_lengths_file
        File? approx_raw_subread_array_lengths

        File? ten_x_metrics_file
        String mas_seq_model

        File workflow_dot_file

        String prefix = ""
        RuntimeAttr? runtime_attr_override
    }



    String nb_name = prefix + "report.ipynb"
    String html_out = prefix + "report.html"
    String pdf_out = prefix + "report.pdf"

    Int disk_size = 20 + 8*ceil((
            size(notebook_template, "GB") +
            size(subreads_stats, "GB") +
            size(ccs_reads_stats, "GB") +
            size(ccs_report_file, "GB") +
            size(raw_ccs_bam_file, "GB") +
            size(array_element_bam_file, "GB") +
            size(ccs_rejected_bam_file, "GB") +
            size(annotated_bam_file, "GB") +
            size(raw_ccs_bam_file, "GB") +
            size(zmw_subread_stats_file, "GB") +
            size(polymerase_read_lengths_file, "GB") +
            size(ten_x_metrics_file, "GB") +
            size(workflow_dot_file, "GB")
        ))

    # Handle the optional files:
    String ten_x_metrics_file_flag = if defined(ten_x_metrics_file) then "true" else "false"
    String zmw_subread_stats_file_flag = if defined(zmw_subread_stats_file) then "true" else "false"
    String polymerase_read_lengths_file_flag = if defined(polymerase_read_lengths_file) then "true" else "false"
    String approx_raw_subread_array_lengths_flag = if defined(approx_raw_subread_array_lengths) then "true" else "false"

    command <<<
        set -euxo pipefail

        # Set up memory logging daemon:
        MEM_LOG_INTERVAL_s=5
        DO_MEMORY_LOG=true
        while $DO_MEMORY_LOG ; do
            date
            date +%s
            cat /proc/meminfo
            sleep $MEM_LOG_INTERVAL_s
        done >> memory_log.txt &
        mem_pid=$!

        # Copy the notebook template to our current folder:
        cp ~{notebook_template} ~{nb_name}

        # Create a template to create the html report with collapsed code:
        echo "{%- extends 'full.tpl' -%}" > hidecode.tpl
        echo "" >> hidecode.tpl
        echo "{% block input_group %}" >> hidecode.tpl
        echo "    {%- if cell.metadata.get('nbconvert', {}).get('show_code', False) -%}" >> hidecode.tpl
        echo "        ((( super() )))" >> hidecode.tpl
        echo "    {%- endif -%}" >> hidecode.tpl
        echo "{% endblock input_group %}" >> hidecode.tpl

        # Set some environment variables for the notebook to read in:
        export DATE_RUN="$(date)"
        export WDL_NAME="PB10xMasSeqArraySingleFlowcell.wdl"
        export REPO_INFO="git@github.com:broadinstitute/long-read-pipelines.git"

        # Prepare the config file:
        rm -f mas-seq_qc_inputs.config

        echo "~{sample_name}" >> mas-seq_qc_inputs.config

        echo "~{subreads_stats}" >> mas-seq_qc_inputs.config
        echo "~{ccs_reads_stats}" >> mas-seq_qc_inputs.config
        echo "~{array_elements_stats}" >> mas-seq_qc_inputs.config
        echo "~{ccs_report_file}" >> mas-seq_qc_inputs.config

        echo "~{raw_ccs_bam_file}" >> mas-seq_qc_inputs.config
        echo "~{array_element_bam_file}" >> mas-seq_qc_inputs.config
        echo "~{array_elements_genome_aligned}" >> mas-seq_qc_inputs.config
        echo "~{ccs_rejected_bam_file}" >> mas-seq_qc_inputs.config

        echo "~{annotated_bam_file}" >> mas-seq_qc_inputs.config

        echo "~{longbow_passed_reads_file}" >> mas-seq_qc_inputs.config
        echo "~{longbow_failed_reads_file}" >> mas-seq_qc_inputs.config

        echo "~{longbow_passed_ccs_reads}" >> mas-seq_qc_inputs.config
        echo "~{longbow_failed_ccs_reads}" >> mas-seq_qc_inputs.config
        echo "~{ccs_reclaimable_reads}" >> mas-seq_qc_inputs.config
        echo "~{ccs_reclaimed_reads}" >> mas-seq_qc_inputs.config
        echo "~{ccs_rejected_longbow_failed_reads}" >> mas-seq_qc_inputs.config
        echo "~{raw_array_elements}" >> mas-seq_qc_inputs.config
        echo "~{ccs_reclaimed_array_elements}" >> mas-seq_qc_inputs.config

        echo "~{zmw_stats_json_gz}" >> mas-seq_qc_inputs.config

        if ~{zmw_subread_stats_file_flag} ; then
            echo "~{zmw_subread_stats_file}" >> mas-seq_qc_inputs.config
        else
            echo "NON-EXISTENT-PLACEHOLDER" >> mas-seq_qc_inputs.config
        fi
        if ~{polymerase_read_lengths_file_flag} ; then
            echo "~{polymerase_read_lengths_file}" >> mas-seq_qc_inputs.config
        else
            echo "NON-EXISTENT-PLACEHOLDER" >> mas-seq_qc_inputs.config
        fi
        if ~{approx_raw_subread_array_lengths_flag} ; then
            echo "~{approx_raw_subread_array_lengths}" >> mas-seq_qc_inputs.config
        else
            echo "NON-EXISTENT-PLACEHOLDER" >> mas-seq_qc_inputs.config
        fi

        if ~{ten_x_metrics_file_flag} ; then
            echo "~{ten_x_metrics_file}" >> mas-seq_qc_inputs.config
        else
            echo "NON-EXISTENT-PLACEHOLDER" >> mas-seq_qc_inputs.config
        fi
        echo "~{mas_seq_model}" >> mas-seq_qc_inputs.config

        echo "~{workflow_dot_file}" >> mas-seq_qc_inputs.config

        # Do the conversion:

        # Run the notebook and populate the notebook itself:
        jupyter nbconvert --execute ~{nb_name} --to notebook --inplace --no-prompt --no-input --clear-output --debug --ExecutePreprocessor.timeout=None

        # Convert the notebook output we created just above here to the HTML report:
        jupyter nbconvert ~{nb_name} --to html --no-prompt --no-input --debug --ExecutePreprocessor.timeout=None

        # Create a tar.gz of the figures directory:
        tar -zcf figures.tar.gz figures

        # Create a dummy pickle for process safety:
        touch dummy.pickle

        # Stop the memory daemon softly.  Then stop it hard if it's not cooperating:
        set +e
        DO_MEMORY_LOG=false
        sleep $(($MEM_LOG_INTERVAL_s  * 2))
        kill -0 $mem_pid &> /dev/null
        if [ $? -ne 0 ] ; then
            kill -9 $mem_pid
        fi
    >>>

    output {
        File populated_notebook = nb_name
        File html_report = html_out
        File figures_tar_gz = "figures.tar.gz"
        File generated_config = "mas-seq_qc_inputs.config"

        Array[File] pickles = glob("*.pickle")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-jupyter_interactive:0.0.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"  # LOCAL here is a local SSD - much faster and even money with normal disk if preemptible
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FindBams {

    meta {
        description: "Find all subreads.bam files in a GCS directory."
    }

    parameter_meta {
        gcs_input_dir: "GCS directory containing subreads.bam files."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        String gcs_input_dir

        RuntimeAttr? runtime_attr_override
    }

    String indir = sub(gcs_input_dir, "/$", "")

    command <<<
        set -euxo pipefail

        gsutil ls "~{indir}/**subreads.bam" > subread_bams.txt
    >>>

    output {
        Array[String] subread_bams = read_lines("subread_bams.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            1,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task ExtractUncorrectedReads {

    meta {
        description: "Extract uncorrected reads from subreads and consensus."
    }

    parameter_meta {
        subreads: "Input subreads BAM."
        consensus: "Input consensus BAM."
        prefix: "Prefix for output files."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        File subreads
        File consensus

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(subreads, "GB") + size(consensus, "GB"))

    command <<<
        set -euxo pipefail

        python3 /usr/local/bin/extract_uncorrected_reads.py -o ~{prefix}.bam ~{subreads} ~{consensus}
    >>>

    output {
        File uncorrected = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.29"
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

task PolishTranscripts {

    meta {
        description: "Polish transcripts."
    }

    parameter_meta {
        bam: "Input BAM file."
        subreads_bam: "Input subreads BAM file."
        subreads_pbi: "Input subreads PBI file."
        prefix: "Prefix for output BAM file."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        File bam
        File subreads_bam
        File subreads_pbi
        String prefix = "polished"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size([bam, subreads_bam], "GB"))

    command <<<
        set -euxo pipefail

        isoseq3 polish ~{bam} ~{subreads_bam} ~{prefix}.bam
    >>>

    output {
        File polished_bam = "~{prefix}.bam"
        File polished_fastq = "~{prefix}.hq.fastq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          24,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-pb:0.1.29"
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

task SummarizeXMLMetadata {

    meta {
        description: "Summarize XML metadata."
    }

    parameter_meta {
        xml: "Input XML metadata."
        runtime_attr_override: "Override default runtime attributes."
    }

    input {
        File xml

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(xml, "GB"))

    command <<<
        set -euxo pipefail

        cat ~{xml} | grep '<pbds:TotalLength>' | sed 's/<pbds:TotalLength>//g' | sed 's/<\/pbds:TotalLength>//' | sed 's/\s*//g' > xml_total_length.txt
        cat ~{xml} | grep '<pbds:NumRecords>' | sed 's/<pbds:NumRecords>//g' | sed 's/<\/pbds:NumRecords>//' | sed 's/\s*//g' > xml_num_records.txt
    >>>

    output {
        Float xml_total_length = read_float("xml_total_length.txt")
        Float xml_num_records = read_float("xml_num_records.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
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

task SplitSoftClippedReads {

    meta {
        description : "Split softclipped reads into multiple reads, each with a different amount of softclipping."
    }

    parameter_meta {
        reads_fastq : "Input reads in fastq format."
        reference_fasta : "Reference fasta file."
        rounds : "Number of rounds of softclipping to perform."
        clipping_threshold : "Threshold for softclipping."
    }

    input {
        File reads_fastq
        File reference_fasta
        Int rounds
        Int clipping_threshold
    }

    Int disk_size = (4 + rounds) * ceil(size(reads_fastq, "GB"))
    String basename = basename(reads_fastq, ".fastq")

    command <<<
        set -euxo pipefail

        for i in {1..~{rounds}}
        do
          if [[ $i -eq 1 ]];
          then
            input_fn=~{reads_fastq}
          else
            input_fn=~{basename}_softclipped_x$((i - 1)).fastq
          fi

          minimap2 --eqx -ax map-ont ~{reference_fasta} ${input_fn} \
            | python /soft_clipper.py --split-read-prefix=x${i} --clipping-threshold=~{clipping_threshold} \
            > ~{basename}_softclipped_x${i}.fastq

        done
    >>>

    output {
       Array[File] split_reads = glob("*_softclipped_*.fastq")
       File most_split_read = "~{basename}_softclipped_x~{rounds}.fastq"
    }

    runtime {
        cpu:                    8
        memory:                 "32 GiB"
        disks:                  "local-disk " +  disk_size + " HDD"
        bootDiskSizeGb:         25
        preemptible:            0
        maxRetries:             0
        docker:                 "quay.io/broad-long-read-pipelines/lr-softclipper:0.5.0"
    }
}

task SplitSoftClippedReadsAssisted {

    meta {
        description: "Split softclipped reads into multiple reads, each with a different amount of softclipping. This version uses an additional reference to help determine the correct amount of softclipping."
    }

    parameter_meta {
        reads_fastq : "Input reads in fastq format."
        reference_fasta : "Reference fasta file."
        rounds : "Number of rounds of softclipping to perform."
        clipping_threshold : "Threshold for softclipping."
        ref_conflict_threshold : "Threshold for softclipping."
        aid_reference_fasta : "Reference fasta file."
    }

    input {
        File reads_fastq
        File reference_fasta
        Int rounds
        Int clipping_threshold
        Int ref_conflict_threshold
        File aid_reference_fasta
    }

    Int disk_size = (14 + 6 + rounds) * ceil(size(reads_fastq, "GB"))
    String basename = basename(reads_fastq, ".fastq")

    command <<<
        set -euxo pipefail

        for i in {1..~{rounds}}
        do
          if [[ $i -eq 1 ]];
          then
            input_fn=~{reads_fastq}
          else
            input_fn=~{basename}_softclipped_x$((i - 1)).fastq
          fi

          minimap2 --eqx -ax map-ont ~{aid_reference_fasta} ${input_fn} > aid_ref.sam
          minimap2 --eqx -ax map-ont ~{reference_fasta} ${input_fn} > ref.sam

          cat ref.sam | python /soft_clipper.py \
               --clipping-threshold=~{clipping_threshold} \
               --split-read-prefix=x${i} \
               --ref=aid_ref.sam \
               --ref-diff-threshold=~{ref_conflict_threshold} \
               --write-ref-conflicts-prefix=conflicts_x${i} \
            > ~{basename}_softclipped_x${i}.fastq

        done
    >>>

    output {
       Array[File] split_reads = glob("*_softclipped_*.fastq")
       Array[File] conflicting_alignments = glob("conflicts_*.bam")
       File most_split_read = "~{basename}_softclipped_x~{rounds}.fastq"
    }

    runtime {
        cpu:                    8
        memory:                 "60 GiB"
        disks:                  "local-disk " +  disk_size + " HDD"
        bootDiskSizeGb:         25
        preemptible:            0
        maxRetries:             0
        docker:                 "quay.io/broad-long-read-pipelines/lr-softclipper:0.5.0"
    }
}

task WriteNamedFile {

    meta {
        description : "Write a file to the given directory with the given name."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        name : "Name of the file to write."
        outdir : "Google cloud path to the destination folder."
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
    }

    input {
        String name
        String outdir
        File? keyfile
    }

    command <<<
        set -euxo pipefail

        touch "~{name}"

        gsutil cp "~{name}" ~{outdir}
    >>>

    #########################

    runtime {
        cpu:                    1
        memory:                 1 + " GiB"
        disks: "local-disk " +  10 + " HDD"
        bootDiskSizeGb:         25
        preemptible:            2
        maxRetries:             2
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-finalize:0.1.2"
    }
}
