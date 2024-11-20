version 1.0

import "../../structs/Structs.wdl" as Structs
import "../Utility/SRUtils.wdl" as SRUTIL
import "../Utility/Utils.wdl" as Utils
import "../Utility/Finalize.wdl" as FF

workflow RemoveSingleOrganismContamination {
    meta {
        author: "Jonn Smith"
        description: "A workflow to remove contamination originating from a single organism from a dataset."
    }

    input {
        File? input_bam
        File? input_bai

        File? fq_end1
        File? fq_end2

        String SM
        String LB
        String platform = "illumina"

        String contaminant_ref_name
        File contaminant_ref_map_file

        String dir_prefix
        String? gcs_out_root_dir

        Boolean DEBUG_MODE = false
    }

    parameter_meta {
        input_bam:                "GCS path to unmapped bam"
        input_bai:                "GCS path to bai index for unmapped bam"

        fq_end1:            "GCS path to end1 of paired-end fastq"
        fq_end2:            "GCS path to end2 of paired-end fastq"

        SM:                 "the value to place in the BAM read group's SM field"
        LB:                 "the value to place in the BAM read group's LB (library) field"
        platform:                 "[default valued] the value to place in the BAM read group's PL (platform) field (default: illumina)"

        contaminant_ref_name:                 "Name of the contaminant genome to be used in output files."
        contaminant_ref_map_file:                 "Table indicating reference sequence and auxillary file locations."

        dir_prefix:         "directory prefix for output files"
        gcs_out_root_dir:    "GCS Bucket into which to finalize outputs.  If no bucket is given, outputs will not be finalized and instead will remain in their native execution location."

        DEBUG_MODE:         "[default valued] enables debugging tasks / subworkflows (default: false)"
    }

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_001_WdlExecutionStartTimestamp { input: }

    # Some basic error handling:
    if (!defined(input_bam) && (!defined(fq_end1) || !defined(fq_end2))) {
        call Utils.StopWorkflow as t_002_NoInputFileProvidedFailure {
            input: reason = "No input file has been provided!  You must provide either an input bam or input fastq1/fastq2 files."
        }
    }
    if (defined(input_bam) && (defined(fq_end1) || defined(fq_end2))) {
        call Utils.StopWorkflow as t_003_TooManyInputsProvidedFailure {
            input: reason = "Too many inputs provided!  You must provide EITHER an input bam OR input fastq1/fastq2 files."
        }
    }

    # Get ref info:
    Map[String, String] ref_map = read_map(contaminant_ref_map_file)

    if (defined(input_bam)) {
        # Convert the given bam to a uBAM (needed for previous aligned data):
        call SRUTIL.RevertSam as t_004_RevertSam {
            input:
                input_bam = select_first([input_bam]),
                prefix = SM + ".revertSam"
        }

        # Convert input SAM/BAM to FASTQ:
        call SRUTIL.BamToFq as t_005_Bam2Fastq {
            input:
                bam = t_004_RevertSam.bam,
                prefix = SM
        }
        call Utils.GetRawReadGroup as t_006_GetRawReadGroup { input: gcs_bam_path = select_first([input_bam]) }
    }

    File fq_e1 = select_first([fq_end1, t_005_Bam2Fastq.fq_end1])
    File fq_e2 = select_first([fq_end2, t_005_Bam2Fastq.fq_end2])
    
    # Align data to contaminant reference:
    call SRUTIL.Bowtie2 as t_007_AlignReads {
        input:
            fq_end1 = fq_e1,
            fq_end2 = fq_e2,

            ref_fasta = ref_map["fasta"],
            ref_fasta_index = ref_map["fai"],
            ref_bowtie_indices = [
                ref_map["1.bt2"],
                ref_map["2.bt2"],
                ref_map["3.bt2"],
                ref_map["4.bt2"],
                ref_map["rev.1.bt2"],
                ref_map["rev.2.bt2"]
            ],
            prefix = SM + ".contaminant_aligned." + contaminant_ref_name,
            rg_id = SM + "_" + LB,
            rg_pl = platform,
            rg_lb = LB,
            rg_sm = SM
    }

    call ExtractReadsWithSamtools as t_008_ExtractDecontaminatedReads {
        input:
            bam = t_007_AlignReads.bam,
            sam_flags = "256",
            extra_args = " -f 12 ",
            prefix = SM + ".decontaminated"
    }

    call ExtractReadsWithSamtools as t_009_ExtractContaminatedReads {
        input:
            bam = t_007_AlignReads.bam,
            sam_flags = "12",
            prefix = SM + ".contaminated_" + contaminant_ref_name + "_reads"
    }

    call SortBamWithoutIndexing as t_010_SortDecontaminatedReads {
        input:
            input_bam = t_008_ExtractDecontaminatedReads.output_bam,
            extra_args = " -n ",
            prefix = SM + ".decontaminated.sorted"
    }

    call Utils.SortSam as t_011_SortContaminatedReads {
        input:
            input_bam = t_009_ExtractContaminatedReads.output_bam,
            prefix = SM + ".contaminated_" + contaminant_ref_name + "_reads.sorted"
    }

    # Convert input SAM/BAM to FASTQ:
    call SRUTIL.BamToFq as t_012_CreateFastqFromDecontaminatedReads {
        input:
            bam = t_010_SortDecontaminatedReads.sorted_bam,
            prefix = SM + ".decontaminated"
    }

    ############################################
    #      _____ _             _ _
    #     |  ___(_)_ __   __ _| (_)_______
    #     | |_  | | '_ \ / _` | | |_  / _ \
    #     |  _| | | | | | (_| | | |/ /  __/
    #     |_|   |_|_| |_|\__,_|_|_/___\___|
    #
    ############################################

    if (defined(gcs_out_root_dir)) {
        # Chosen because it's a relatively small file.
        File keyfile = t_012_CreateFastqFromDecontaminatedReads.fq_unpaired

        # Create an outdir:
        String concrete_gcs_out_root_dir = select_first([gcs_out_root_dir])
        String outdir = if DEBUG_MODE then sub(concrete_gcs_out_root_dir, "/$", "") + "/RemoveSingleOrganismContamination/~{dir_prefix}/" + t_001_WdlExecutionStartTimestamp.timestamp_string else sub(concrete_gcs_out_root_dir, "/$", "") + "/RemoveSingleOrganismContamination/~{dir_prefix}"

        call FF.FinalizeToFile as t_013_FinalizeContaminatedBam { input: outdir = outdir, file = t_011_SortContaminatedReads.output_bam, keyfile = keyfile }
        call FF.FinalizeToFile as t_014_FinalizeContaminatedBamIndex { input: outdir = outdir, file = t_011_SortContaminatedReads.output_bam_index, keyfile = keyfile }
        call FF.FinalizeToFile as t_015_FinalizeDecontaminatedFq1 { input: outdir = outdir, file = t_012_CreateFastqFromDecontaminatedReads.fq_end1, keyfile = keyfile }
        call FF.FinalizeToFile as t_016_FinalizeDecontaminatedFq2 { input: outdir = outdir, file = t_012_CreateFastqFromDecontaminatedReads.fq_end2, keyfile = keyfile }
        call FF.FinalizeToFile as t_017_FinalizeDecontaminatedUnpaired { input: outdir = outdir, file = t_012_CreateFastqFromDecontaminatedReads.fq_unpaired, keyfile = keyfile }
    }

 
    ############################################
    #      ___        _               _
    #     / _ \ _   _| |_ _ __  _   _| |_
    #    | | | | | | | __| '_ \| | | | __|
    #    | |_| | |_| | |_| |_) | |_| | |_
    #     \___/ \__,_|\__| .__/ \__,_|\__|
    #                    |_|
    ############################################

    output {
        File contaminated_bam = select_first([t_013_FinalizeContaminatedBam.gcs_path, t_011_SortContaminatedReads.output_bam])
        File contaminated_bam_index = select_first([t_014_FinalizeContaminatedBamIndex.gcs_path, t_011_SortContaminatedReads.output_bam_index])

        File decontaminated_fq1 = select_first([t_015_FinalizeDecontaminatedFq1.gcs_path, t_012_CreateFastqFromDecontaminatedReads.fq_end1])
        File decontaminated_fq2 = select_first([t_016_FinalizeDecontaminatedFq2.gcs_path, t_012_CreateFastqFromDecontaminatedReads.fq_end2])
        File decontaminated_unpaired = select_first([t_017_FinalizeDecontaminatedUnpaired.gcs_path, t_012_CreateFastqFromDecontaminatedReads.fq_unpaired])
    }
}

task ExtractReadsWithSamtools {
    meta {
        description : "Filter reads based on sam flags.  Reads with ANY of the given flags will be removed from the given dataset.  Does not sort or index the results."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File bam
        String sam_flags

        String extra_args = ""

        String prefix = "filtered_reads"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:   "BAM file to be filtered."
        sam_flags: "Flags for which to remove reads.  Reads with ANY of the given flags will be removed from the given dataset."
        prefix : "[Optional] Prefix string to name the output file (Default: filtered_reads)."
    }

    Int disk_size = 20 + ceil(11 * size(bam, "GiB"))

    command <<<

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        samtools view -h -b -F ~{sam_flags} -@$np ~{extra_args} ~{bam} > ~{prefix}.bam
    >>>

    output {
        File output_bam = "~{prefix}.bam"
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

task SortBamWithoutIndexing {
    input {
        File input_bam
        String prefix = "sorted"

        String? extra_args = ""

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        input_bam: "input BAM"
        prefix:    "[default-valued] prefix for output BAM"
    }

    Int disk_size = 10 + 10*ceil(size(input_bam, "GB"))

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        samtools sort ~{extra_args} -@$num_core -o ~{prefix}.bam ~{input_bam}
    >>>

    output {
        File sorted_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
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
