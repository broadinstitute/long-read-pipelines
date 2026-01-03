version 1.0

import "../../structs/Structs.wdl"

task BamToFq {
    input {
        File bam
        File? bam_index
        
        File? reference_fasta
        File? reference_fasta_index
        File? reference_dict

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    String ref_arg = if defined(reference_fasta) then " --reference " else ""

    Int disk_size = 10 + 20*ceil(size(bam, "GB"))

    command <<<

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        if [[ ${np} -gt 2 ]] ; then
            np=$((np-1))
        fi

        set -euxo pipefail

        # Have samtools sort use all but one of our processors:
        # NOTE: the `@` options is for ADDITIONAL threads, not the total number of threads.
        samtools sort -@$((np-1)) -n ~{ref_arg} ~{reference_fasta} ~{bam} -O bam -o tmp.bam
        
        # Have samtools bam2fq use all but one of our processors:
        # NOTE: the `@` options is for ADDITIONAL threads, not the total number of threads.
        samtools bam2fq -@$((np-1)) \
            -n \
            -s /dev/null \
            -c 2 \
            ~{ref_arg} ~{reference_fasta} \
            -1 ~{prefix}.end1.fq.gz \
            -2 ~{prefix}.end2.fq.gz \
            -0 ~{prefix}.unpaired.fq.gz \
            tmp.bam
    >>>

    output {
        File fq_end1 = "~{prefix}.end1.fq.gz"
        File fq_end2 = "~{prefix}.end2.fq.gz"
        File fq_unpaired = "~{prefix}.unpaired.fq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
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

task FixMate {
    input {
        File input_bam
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(input_bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools fixmate ~{input_bam} ~{prefix}.bam
    >>>

    output {
        File bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
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

task Bam2FqPicard {
    input {
        File bam
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        java -Xms8192m -Xmx30768m -jar /usr/picard/picard.jar \
            SamToFastq \
            INPUT=~{bam} \
            FASTQ=~{prefix}.fastq \
            INTERLEAVE=true \
            NON_PF=true
    >>>

    output {
        File fastq = "~{prefix}.fastq"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
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

task BwaMem2 {
    input {
        File fq_end1
        File fq_end2

        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File ref_0123
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac

        String? read_group

        String prefix = "out"

        Boolean mark_short_splits_as_secondary = false

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(fq_end1, "GB"))
                      + 4*ceil(size(fq_end2, "GB"))
                      + 4*ceil(size(ref_fasta, "GB"))
                      + 4*ceil(size(ref_fasta_index, "GB"))
                      + 4*ceil(size(ref_dict, "GB"))
                      + 4*ceil(size(ref_amb, "GB"))
                      + 4*ceil(size(ref_ann, "GB"))
                      + 4*ceil(size(ref_bwt, "GB"))
                      + 4*ceil(size(ref_pac, "GB"))
                      + 4*ceil(size(ref_0123, "GB"))

    String rg_arg = if defined(read_group) then " -R " else "" 
    String rg_val_wrapper = if defined(read_group) then "'" else ""

    command <<<
        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        if [[ ${np} -gt 2 ]] ; then
            np=$((np-1))
        fi

        # Breakdown of the arguments:
        # -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility) []
        # -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]
        # -t INT        number of threads [1]
        # -Y            use soft clipping for supplementary alignments
        # -R STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
        # -M            mark shorter split hits as secondary

        bwa-mem2 mem \
            -K 100000000 \
            -v 3 \
            -t ${np} \
            -Y \
            ~{rg_arg}~{rg_val_wrapper}~{default="" read_group}~{rg_val_wrapper} \
            ~{true='-M' false="" mark_short_splits_as_secondary} \
            ~{ref_fasta} \
            ~{fq_end1} \
            ~{fq_end2} | \
        samtools view -1 - > ~{prefix}.bam
    >>>

    output {
        File bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"
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

task Bowtie2 {

    parameter_meta {
        prefix: "Prefix for output files."

        fq_end1: "FASTQ file containing end 1 of reads."
        fq_end2: "FASTQ file containing end 2 of reads."

        ref_fasta: "Reference FASTA file for the genome of the organism to which to align reads."
        ref_fasta_index: "Reference FASTA index file for the genome of the organism to which to align reads."
        ref_bowtie_indices: "Bowtie2 indices for the given genome."

        runtime_attr_override: "Optional override for runtime attributes."
    }

    input {
        String prefix

        File fq_end1
        File fq_end2

        File ref_fasta
        File ref_fasta_index
        Array[File] ref_bowtie_indices

        Boolean skip_sort = false

        String? rg_id
        String? rg_pl
        String? rg_lb
        String? rg_sm

        RuntimeAttr? runtime_attr_override
    }

    Int bowtie_index_size = ceil(size(ref_bowtie_indices, "GB"))
    Int disk_size = 10 + 10 * (bowtie_index_size + ceil(size([fq_end1, fq_end2, ref_fasta, ref_fasta_index], "GB")))

    String rgid_cmd = if defined(rg_id) then " --rg-id " else ""
    String rg_cmd = if (defined(rg_pl) || defined(rg_lb) || defined(rg_sm)) then " --rg " else ""

    String rg_id_val = if defined(rg_id) then select_first([rg_id]) else ""
    String rg_pl_val = if defined(rg_pl) then "PL:" + select_first([rg_pl]) else ""
    String rg_lb_val = if defined(rg_lb) then "LB:" + select_first([rg_lb]) else ""
    String rg_sm_val = if defined(rg_sm) then "SM:" + select_first([rg_sm]) else ""

    command <<<
        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        if [[ ${np} -gt 2 ]] ; then
            np=$((np-1))
        fi

        # Move the bowtie2 index files to the same directory as the reference fasta file:
        ref_dir=$(dirname ~{ref_fasta})
        while read f ; do 
            indx_dirname=$( dirname ${f} )
            if [[ "${indx_dirname}" != "${ref_dir}" ]] ; then
                mv ${f} ${ref_dir}
            fi
        done < ~{write_lines(ref_bowtie_indices)}

        # Get the basename of the bowtie2 index file so we can reference it in the bowtie2 command:
        bowtie2_index_basename=$(echo "~{ref_bowtie_indices[0]}" | grep bt2 | head -n1 | sed 's@\.[a-zA-Z0-9]*\.bt2@@')

        echo "Aligning to genome:"
        bowtie2 -x ${bowtie2_index_basename} \
            --threads ${np} \
            -1 ~{fq_end1} \
            -2 ~{fq_end2} \
            ~{rgid_cmd} "~{rg_id_val}" \
            ~{rg_cmd} "~{rg_pl_val}" \
            ~{rg_cmd} "~{rg_lb_val}" \
            ~{rg_cmd} "~{rg_sm_val}" | \
        samtools view -bh --no-PG - > tmp.bam

        # Now sort the output:
        if ~{skip_sort} ; then
            mv tmp.bam ~{prefix}.bam
        else
            samtools sort -@$((np-1)) tmp.bam > ~{prefix}.bam
        fi

        samtools index -@$((np-1)) ~{prefix}.bam
    >>>

    output {
        File bam = "~{prefix}.bam"
        File bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/sr-alternate-tools:0.0.1"
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

task MergeBamAlignment {
    input {
        File aligned_bam
        File unaligned_bam

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(aligned_bam, "GB"))
                      + 4*ceil(size(unaligned_bam, "GB"))
                      + 4*ceil(size(ref_fasta, "GB"))
                      + 4*ceil(size(ref_fasta_index, "GB"))
                      + 4*ceil(size(ref_dict, "GB"))

    command <<<
        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        np=$((np-1))

        java -Dsamjdk.compression_level=2 -Xms8192m -Xmx30768m -jar /usr/picard/picard.jar \
            MergeBamAlignment \
            VALIDATION_STRINGENCY=SILENT \
            EXPECTED_ORIENTATIONS=FR \
            ATTRIBUTES_TO_RETAIN=X0 \
            ATTRIBUTES_TO_REMOVE=NM \
            ATTRIBUTES_TO_REMOVE=MD \
            ALIGNED_BAM=~{aligned_bam} \
            UNMAPPED_BAM=~{unaligned_bam} \
            OUTPUT=~{prefix}.bam \
            REFERENCE_SEQUENCE=~{ref_fasta} \
            SORT_ORDER="unsorted" \
            IS_BISULFITE_SEQUENCE=false \
            ALIGNED_READS_ONLY=false \
            CLIP_ADAPTERS=false \
            MAX_RECORDS_IN_RAM=2000000 \
            ADD_MATE_CIGAR=true \
            MAX_INSERTIONS_OR_DELETIONS=-1 \
            PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            PROGRAM_RECORD_ID="bwa-mem2" \
            PROGRAM_GROUP_VERSION="2.2.1" \
            PROGRAM_GROUP_COMMAND_LINE="bwa-mem2 mem -K 100000000 -p -v 3 -t 15 -Y" \
            PROGRAM_GROUP_NAME="bwa-mem2" \
            UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
            ALIGNER_PROPER_PAIR_FLAGS=true \
            UNMAP_CONTAMINANT_READS=true \
            ADD_PG_TAG_TO_READS=false
    >>>

    output {
        File bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
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

task MarkDuplicates {
    input {
        File input_bam

        String prefix

        # The program default for READ_NAME_REGEX is appropriate in nearly every case.
        # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
        # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
        String? read_name_regex

        Float? sorting_collection_size_ratio

        RuntimeAttr? runtime_attr_override
    }

    Int compression_level = 2

    Int disk_size = 1 + 4*ceil(size(input_bam, "GB"))

    # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
    # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
    # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

    command <<<
        set -euxo pipefail

        tot_mem_mb=$(free -m | grep '^Mem' | awk '{print $2}')
        java_memory_size_mb=$((tot_mem_mb-5120))

        java -Dsamjdk.compression_level=~{compression_level} -Xms${java_memory_size_mb}m -jar /usr/picard/picard.jar \
            MarkDuplicates \
            INPUT=~{input_bam} \
            OUTPUT=~{prefix}.bam \
            METRICS_FILE=~{prefix}.metrics.txt \
            VALIDATION_STRINGENCY=SILENT \
            ~{"READ_NAME_REGEX=" + read_name_regex} \
            ~{"SORTING_COLLECTION_SIZE_RATIO=" + sorting_collection_size_ratio} \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
            ASSUME_SORT_ORDER="queryname" \
            CLEAR_DT="false" \
            ADD_PG_TAG_TO_READS=false
    >>>

    output {
        File bam = "~{prefix}.bam"
        File metrics = "~{prefix}.metrics.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
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

task MarkDuplicatesAndSort {
    input {
        File input_bam

        String prefix

        # The program default for READ_NAME_REGEX is appropriate in nearly every case.
        # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
        # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
        String? read_name_regex

        Float? sorting_collection_size_ratio

        RuntimeAttr? runtime_attr_override
    }

    Int compression_level = 2

    Int disk_size = 1 + 4*ceil(size(input_bam, "GB"))

    # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
    # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
    # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

    command <<<
        set -euxo pipefail

        tot_mem_mb=$(free -m | grep '^Mem' | awk '{print $2}')
        java_memory_size_mb=$((tot_mem_mb-5120))

        java -Dsamjdk.compression_level=~{compression_level} -Xms${java_memory_size_mb}m -jar /usr/picard/picard.jar \
            MarkDuplicates \
            INPUT=~{input_bam} \
            OUTPUT=/dev/stdout \
            METRICS_FILE=~{prefix}.metrics.txt \
            VALIDATION_STRINGENCY=SILENT \
            ~{"READ_NAME_REGEX=" + read_name_regex} \
            ~{"SORTING_COLLECTION_SIZE_RATIO=" + sorting_collection_size_ratio} \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
            ASSUME_SORT_ORDER="queryname" \
            CLEAR_DT="false" \
            ADD_PG_TAG_TO_READS=false | \
        java -jar /usr/picard/picard.jar SortSam \
            INPUT=/dev/stdin \
            OUTPUT=~{prefix}.bam \
            --SORT_ORDER coordinate \
            CREATE_INDEX=true
    >>>

    output {
        File bam = "~{prefix}.bam"
        File bai = "~{prefix}.bam.bai"
        File metrics = "~{prefix}.metrics.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
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


# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
    input {
        File input_bam
        File input_bam_index

        File ref_dict
        File ref_fasta
        File ref_fasta_index

        File known_sites_vcf
        File known_sites_index

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(input_bam, "GB"))
                      + 4*ceil(size(input_bam_index, "GB"))
                      + 2*ceil(size(ref_dict, "GB"))
                      + 2*ceil(size(ref_fasta, "GB"))
                      + 2*ceil(size(ref_fasta_index, "GB"))
                      + 2*ceil(size(known_sites_vcf, "GB"))
                      + 2*ceil(size(known_sites_index, "GB"))

    parameter_meta {
        input_bam: {
            localization_optional: true
        }
    }

    command {

        set -euxo pipefail

        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal -Xms5000m" \
            BaseRecalibrator \
            -R ~{ref_fasta} \
            -I ~{input_bam} \
            --use-original-qualities \
            -O ~{prefix}.txt \
            --known-sites ~{known_sites_vcf}

    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.5.0.0"
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
        File recalibration_report = "~{prefix}.txt"
    }
}

task ApplyBQSR {
    input {
        File input_bam
        File input_bam_index

        File ref_dict
        File ref_fasta
        File ref_fasta_index

        File recalibration_report

        Boolean bin_base_qualities = true
        Boolean emit_original_quals = true

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int compression_level = 2
    Int java_memory_size_mb = 30768

    parameter_meta {
        input_bam: {
            localization_optional: true
        }
    }

    Int disk_size = 1 + 4*ceil(size(input_bam, "GB"))
                      + 4*ceil(size(input_bam_index, "GB"))
                      + 2*ceil(size(ref_dict, "GB"))
                      + 2*ceil(size(ref_fasta, "GB"))
                      + 2*ceil(size(ref_fasta_index, "GB"))
                      + 2*ceil(size(recalibration_report, "GB"))

    command <<<
        set -euxo pipefail

        gatk --java-options "-XX:+PrintFlagsFinal \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=~{compression_level} -Xms8192m -Xmx~{java_memory_size_mb}m" \
            ApplyBQSR \
            --create-output-bam-md5 \
            --add-output-sam-program-record \
            -R ~{ref_fasta} \
            -I ~{input_bam} \
            --use-original-qualities \
            -O ~{prefix}.bam \
            -bqsr ~{recalibration_report} \
            --emit-original-quals ~{emit_original_quals} \
            ~{true='--static-quantized-quals 10' false='' bin_base_qualities} \
            ~{true='--static-quantized-quals 20' false='' bin_base_qualities} \
            ~{true='--static-quantized-quals 30' false='' bin_base_qualities} \
            --allow-missing-read-group true \

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        samtools index -@${np} ~{prefix}.bam
    >>>
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        # docker:             "us.gcr.io/broad-gatk/gatk:4.5.0.0"
        # Temporary snapshot build for testing the fix for BQSR issue https://github.com/broadinstitute/gatk/issues/6242
        docker:             "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk-remote-builds:jonn-4dd794b3f4e4e4e6a86f309a5ea1b580bf774b7c-4.5.0.0-48-g4dd794b3f" 
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
        File recalibrated_bam = "~{prefix}.bam"
        File recalibrated_bai = "~{prefix}.bam.bai"
    }
}

task RunBaseRecalibratorAndApplyBQSR {
    input {
        File input_bam
        File input_bam_index

        File ref_dict
        File ref_fasta
        File ref_fasta_index

        File known_sites_vcf
        File known_sites_index

        Boolean bin_base_qualities = true
        Boolean emit_original_quals = true

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int compression_level = 2
    Int java_memory_size_mb = 30768

    parameter_meta {
        input_bam: { localization_optional: true }
        known_sites_vcf: { localization_optional: true }
    }

    Int disk_size = 5 + 4*ceil(size(input_bam, "GB"))
                      + 4*ceil(size(input_bam_index, "GB"))
                      + 2*ceil(size(ref_dict, "GB"))
                      + 2*ceil(size(ref_fasta, "GB"))
                      + 2*ceil(size(ref_fasta_index, "GB"))
                      + 2*ceil(size(known_sites_vcf, "GB"))
                      + 2*ceil(size(known_sites_index, "GB"))

    command <<<
        set -euxo pipefail

        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal -Xms5000m" \
            BaseRecalibrator \
            -R ~{ref_fasta} \
            -I ~{input_bam} \
            --use-original-qualities \
            -O ~{prefix}.baseRecalibratorReport.txt \
            --known-sites ~{known_sites_vcf}

        gatk --java-options "-XX:+PrintFlagsFinal \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=~{compression_level} -Xms8192m -Xmx~{java_memory_size_mb}m" \
            ApplyBQSR \
            --create-output-bam-md5 \
            --add-output-sam-program-record \
            -R ~{ref_fasta} \
            -I ~{input_bam} \
            --use-original-qualities \
            -O ~{prefix}.bam \
            -bqsr ~{prefix}.baseRecalibratorReport.txt \
            --emit-original-quals ~{emit_original_quals} \
            ~{true='--static-quantized-quals 10' false='' bin_base_qualities} \
            ~{true='--static-quantized-quals 20' false='' bin_base_qualities} \
            ~{true='--static-quantized-quals 30' false='' bin_base_qualities} \
            --allow-missing-read-group true \

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        samtools index -@${np} ~{prefix}.bam
    >>>
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        # docker:             "us.gcr.io/broad-gatk/gatk:4.5.0.0"
        # Temporary snapshot build for testing the fix for BQSR issue https://github.com/broadinstitute/gatk/issues/6242
        docker:             "us.gcr.io/broad-dsde-methods/broad-gatk-snapshots/gatk-remote-builds:jonn-4dd794b3f4e4e4e6a86f309a5ea1b580bf774b7c-4.5.0.0-48-g4dd794b3f" 
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
        File recalibrated_bam = "~{prefix}.bam"
        File recalibrated_bai = "~{prefix}.bam.bai"
        File recalibration_report = "~{prefix}.baseRecalibratorReport.txt"
    }
}

task RevertSam {
    input {
        File input_bam
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int compression_level = 2
    Int java_memory_size_mb = 30768

    Int disk_size = 1 + 20*ceil(size(input_bam, "GB"))

    # As documented on the GATK website:
    # https://gatk.broadinstitute.org/hc/en-us/articles/4403687183515--How-to-Generate-an-unmapped-BAM-from-FASTQ-or-aligned-BAM
    command {

        set -euxo pipefail

        java -Dsamjdk.compression_level=~{compression_level} -Xms~{java_memory_size_mb}m -jar /usr/picard/picard.jar \
            RevertSam \
            INPUT=~{input_bam} \
            OUTPUT=~{prefix}.bam \
            SANITIZE=true \
            MAX_DISCARD_FRACTION=0.005 \
            ATTRIBUTE_TO_CLEAR=XT \
            ATTRIBUTE_TO_CLEAR=XN \
            ATTRIBUTE_TO_CLEAR=AS \
            ATTRIBUTE_TO_CLEAR=OC \
            ATTRIBUTE_TO_CLEAR=OP \
            SORT_ORDER=queryname \
            RESTORE_ORIGINAL_QUALITIES=true \
            REMOVE_DUPLICATE_INFORMATION=true \
            REMOVE_ALIGNMENT_INFORMATION=true
    }

    output {
        File bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
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

task ComputeBamStats {
    input {
        File bam_file
        Int? qual_threshold

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(bam_file, "GB"))
    String qual_thresh_arg = if defined(qual_threshold) then " -q " else ""

    String qual_stats_file_decoration = if defined(qual_threshold) then ".q" + select_first([qual_threshold]) else ""

    String stats_file_name = basename(bam_file, ".bam") + ".stats_map" + qual_stats_file_decoration + ".txt"

    command <<<
        set -euxo pipefail

        python3 /python/compute_sr_stats.py \
            ~{qual_thresh_arg}~{default="" qual_threshold} \
            ~{bam_file} \
        | tee ~{stats_file_name}

    >>>

    output {
        Map[String, Float] results = read_map(stats_file_name)
        File results_file = stats_file_name
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.2"
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

task MergeVCFs {
    meta {
        description: "Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs"
    }

    input {
        Array[File] input_vcfs
        Array[File] input_vcfs_indexes
        String prefix

        Boolean is_gvcf = false

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(input_vcfs, "GiB") * 2.5) + 10

    String gvcf_decorator = if is_gvcf then ".g" else ""

    command <<<
        set -euxo pipefail

        java -Xms2000m -Xmx2500m -jar /usr/picard/picard.jar \
          MergeVcfs \
          INPUT=~{sep=' INPUT=' input_vcfs} \
          OUTPUT=~{prefix}~{gvcf_decorator}.vcf.gz
    >>>

    output {
        File output_vcf = "~{prefix}~{gvcf_decorator}.vcf.gz"
        File output_vcf_index = "~{prefix}~{gvcf_decorator}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             3,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
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


task IndexFeatureFile {
    meta {
        description: "Create a Tribble index for a feature file using GATK.  Feature files are defined inside GATK and include VCF, BED, GTF, and other files."
    }

    input {
        File feature_file

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(feature_file, "GiB") * 2) + 10

    String fname = basename(feature_file)

    command <<<
        set -euxo pipefail

        mv ~{feature_file} ~{fname}
        gatk --java-options "-Xmx1500m" \
            IndexFeatureFile \
                -I ~{fname}
    >>>

    output {
        File index = "~{fname}.idx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.5.0.0"
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


task RevertBaseQualities {
    meta {
        description: "Replace base qualities in the bam file with those located in the `OQ` tag.  If  `ApplyBQSR` has not been run on the given bam file, no changes are made and the original file is returned."
    }

    input {
        File bam
        File? bai

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(bam, "GiB") * 4) + 10

    command <<<

        # Check if the input bam has been run through `ApplyBQSR`.
        # If not, we can just return the input bam.
        samtools view -H ~{bam} | grep '^@PG' > header.pg.txt

        grep -q 'ID:GATK ApplyBQSR' header.pg.txt > applybqsr.pg.txt
        rv=$?

        if [[ $rv -eq 0 ]] && grep -q '\-\-emit-original-quals' applybqsr.pg.txt ; then
            # OK - our data has had it's base quality scores recalibrated.
            # We must revert them:

            set -euxo pipefail

            gatk \
                RevertBaseQualityScores \
                    -I ~{bam} \
                    -O ~{prefix}.bam
        else
            # BQSR was not applied.  Just copy input -> output
            cp ~{bam} ~{prefix}.bam
            if [[ ! -e '~{bai}' ]] ; then
                samtools index ~{prefix}.bam
            else
                cp ~{bai} ~{prefix}.bam.bai
            fi
        fi
    >>>

    output {
        File bam_out = "~{prefix}.bam"
        File bai_out = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.5.0.0"
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
