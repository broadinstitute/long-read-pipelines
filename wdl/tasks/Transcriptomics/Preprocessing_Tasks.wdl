version 1.0

import "../../structs/Structs.wdl"

task DemuxMasSeqDataByIndex {

    meta {
        description : "This workflow will split MAS-seq data that is indexed with a 10bp sequence at the 3' end."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        array_bam : "Bam file containing annotated MAS-seq array reads that contain a 10bp index near the 3' end.."
    }

    input {
        File array_bam
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
        docker: "us.gcr.io/broad-dsp-lrma/lr-transcript_utils:0.0.14"
        memory: 4 + " GiB"
        disks: "local-disk " + disk_size + " HDD"
        bootDiskSizeGb: 25
        preemptible: 0
        cpu: 2
    }
}

task ConvertSplicedBamToGff {
    meta {
        description : "Convert a given splice aligned bam file into a gff file."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        bam : "Bamfile to be converted to gff."
    }

    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    String base_name = basename(bam, ".bam")

    Int disk_size_gb = 10 + 5*ceil(size(bam, "GB"))

    command <<<
        spliced_bam2gff -S -M ~{bam} > ~{base_name}.gff
    >>>

    output {
        File gff = "~{base_name}.gff"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-splicedbam2gff:0.0.1"
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

task GffCompare {
    meta {
        description : "Compare two GFF files"
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        gff_ref : "Gff file to be used as a reference."
        gff_query : "Gff file to be used as a query (compared against the gff_ref)."
        ref_fasta : "Reference fasta file."
        ref_fasta_index : "Reference fasta file index."
    }

    input {
        File gff_ref
        File gff_query
        File ref_fasta
        File? ref_fasta_index

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 2*ceil(size(gff_ref, "GB")) + 2*ceil(size(gff_query, "GB")) + 2*ceil(size(ref_fasta, "GB"))

    String query_base_name = basename(gff_query)

    command <<<
        # Because of how gffcompare works, we need to move the query file to our PWD:

        mv -v ~{gff_query} .

        time /gffcompare/gffcompare \
             -V \
             -r ~{gff_ref} \
             -s ~{ref_fasta} \
             ~{query_base_name} &> ~{prefix}.~{query_base_name}.gffcmp.log

        # Rename some output files so we can disambiguate them later:
        mv gffcmp.~{query_base_name}.refmap ~{prefix}.gffcmp.~{query_base_name}.refmap
        mv gffcmp.~{query_base_name}.tmap ~{prefix}.gffcmp.~{query_base_name}.tmap

        mv gffcmp.tracking ~{prefix}.~{query_base_name}.gffcmp.tracking
        mv gffcmp.loci ~{prefix}.~{query_base_name}.gffcmp.loci
        mv gffcmp.annotated.gtf ~{prefix}.~{query_base_name}.gffcmp.annotated.gtf
        mv gffcmp.stats ~{prefix}.~{query_base_name}.gffcmp.stats
    >>>

    output {
        File refmap        = "~{prefix}.gffcmp.~{query_base_name}.refmap"
        File tmap          = "~{prefix}.gffcmp.~{query_base_name}.tmap"

        File tracking      = "~{prefix}.~{query_base_name}.gffcmp.tracking"
        File loci          = "~{prefix}.~{query_base_name}.gffcmp.loci"
        File annotated_gtf = "~{prefix}.~{query_base_name}.gffcmp.annotated.gtf"
        File stats         = "~{prefix}.~{query_base_name}.gffcmp.stats"
        File log           = "~{prefix}.~{query_base_name}.gffcmp.log"
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

task RestoreOriginalReadNames {
    meta {
        description : "Copies the contents of the XM tag to the read name and sets the XM tag to the read name."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        bam : "Bam file in which to restore the original read names."
    }

    input {
        File bam
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 2*ceil(size(bam, "GB"))

    command <<<
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        time /python_scripts/restore_original_read_names.py ~{bam}

        # Rename some output files so we can disambiguate them later:
        mv out.original_read_names_restored.bam ~{prefix}.bam

        samtools index -@${np} ~{prefix}.bam
    >>>

    output {
        File bam_out = "~{prefix}.bam"
        File bai_out = "~{prefix}.bam.bai"
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

task CorrectUmisWithSetCover {
    meta {
        description : "Corrects the UMIs in the given reads using a set cover algorithm"
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    parameter_meta {
        bam : "Bam file for which to correct UMIs."
        prefix : "Prefix to assign to output files."
    }

    input {
        File bam
        String prefix = "out"

        Boolean is_extracted = true

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 10*ceil(size(bam, "GB"))

    String is_extracted_arg = if (is_extracted)  then "--pre-extracted" else ""

    command <<<
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        time /python_scripts/umi_correction.py  \
            --input_bam  ~{bam} \
            --output_bam  ~{prefix}.corrected_umis.unsorted.bam \
            --filtered_bam ~{prefix}.corrected_umis.failed_filters.bam \
            ~{is_extracted_arg} \
            --config /python_scripts/umi_correction.yaml

        samtools sort ~{prefix}.corrected_umis.unsorted.bam > ~{prefix}.corrected_umis.bam
        samtools index -@${np} ~{prefix}.corrected_umis.bam
    >>>

    output {
        File corrected_umi_reads = "~{prefix}.corrected_umis.bam"
        File corrected_umi_reads_index = "~{prefix}.corrected_umis.bam.bai"
        File uncorrected_umi_reads = "~{prefix}.corrected_umis.failed_filters.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             32,
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
