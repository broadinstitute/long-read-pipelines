version 1.0

import "../../tasks/Z_One_Off_Analyses/BroadOnPremMalariaPipelineTasks.wdl" as BroadOnPremMalariaPipelineTasks
import "../../tasks/Utility/SRUtils.wdl" as SRUTIL
import "../../tasks/Utility/Utils.wdl" as Utils
import "../../structs/Structs.wdl"


workflow BroadOnPremMalariaPipeline_1_Alignment {
    meta {
        desciption: "Recreation of the first and second steps of the Broad OnPrem Malaria Pipeline.  Aligns reads and filters out human contamination, then calls variants and creates both a GVCF file and a recalibrated variants file."
    }
    input {
        String sample_name

        File? bam
        File? bai

        File? fq_end1
        File? fq_end2

        File ref_map_file
        String contaminant_ref_name
        File contaminant_ref_map_file

        File resource_vcf_7g8_gb4
        File resource_vcf_7g8_gb4_index
        File resource_vcf_hb3_dd2
        File resource_vcf_hb3_dd2_index
        File resource_vcf_3d7_hb3
        File resource_vcf_3d7_hb3_index
    }

    # Sanity checks:
    if ((!defined(bam) && !defined(bai)) && (!defined(fq_end1) && !defined(fq_end2))) {
        call Error as UserErrorMissingInputFiles {
            input:
                message = "Either 'bam' and 'bai' or 'fq_end1' and 'fq_end2' must be provided."
        }
    }

    # Get ref info:
    Map[String, String] ref_map = read_map(ref_map_file)
    Map[String, String] contaminant_ref_map = read_map(contaminant_ref_map_file)

    # Call our timestamp so we can store outputs without clobbering previous runs:
    call Utils.GetCurrentTimestampString as t_001_WdlExecutionStartTimestamp { input: }

    if (defined(bam)) {
        # Convert the given bam to a uBAM (needed for previous aligned data):
        call SRUTIL.RevertSam as t_002_RevertSam {
            input:
                input_bam = select_first([bam]),
                prefix = sample_name + ".revertSam"
        }

        # Convert input SAM/BAM to FASTQ:
        call SRUTIL.BamToFq as t_003_Bam2Fastq {
            input:
                bam = t_002_RevertSam.bam,
                prefix = sample_name
        }

        call Utils.GetRawReadGroup as t_004_GetRawReadGroup { input: gcs_bam_path = select_first([bam]) }
    }

    File fq_e1 = select_first([fq_end1, t_003_Bam2Fastq.fq_end1])
    File fq_e2 = select_first([fq_end2, t_003_Bam2Fastq.fq_end2])

    # 1 - Filter out human reads:
    call FilterOutTargetOrganismReads as t_005_FilterOutHumanReads {
        input:
            prefix = sample_name,
            fq_end1 = fq_e1,
            fq_end2 = fq_e2,
            contaminant_reference_fasta = contaminant_ref_map["fasta"],
            contaminant_reference_fai = contaminant_ref_map["fai"],
            contaminant_reference_bowtie_indices = [
                contaminant_ref_map["1.bt2"],
                contaminant_ref_map["2.bt2"],
                contaminant_ref_map["3.bt2"],
                contaminant_ref_map["4.bt2"],
                contaminant_ref_map["rev.1.bt2"],
                contaminant_ref_map["rev.2.bt2"]
            ],
            contaminant_organism_label = contaminant_ref_name
    }

    # 2 - Align to Plasmodium reference:
    String RG = "@RG\tID:FLOWCELL_~{sample_name}\tSM:~{sample_name}\tPL:ILLUMINA\tLB:LIB_~{sample_name}"

    # 2a - Align reads to reference with BWA-MEM2:
    call SRUTIL.BwaMem2 as t_006_AlignReads {
        input:
            fq_end1 = t_005_FilterOutHumanReads.fq1,
            fq_end2 = t_005_FilterOutHumanReads.fq2,
            ref_fasta = ref_map["fasta"],
            ref_fasta_index = ref_map["fai"],
            ref_dict = ref_map["dict"],
            ref_0123 = ref_map["0123"],
            ref_amb = ref_map["amb"],
            ref_ann = ref_map["ann"],
            ref_bwt = ref_map["bwt"],
            ref_pac = ref_map["pac"],
            mark_short_splits_as_secondary = true,
            read_group = RG,
            prefix = sample_name + ".aligned"
    }

    # 2b - Sort Bam:
    call Utils.SortSam as t_008_SortAlignedBam {
        input:
            input_bam = t_006_AlignReads.bam,
            prefix = sample_name + ".aligned.sorted"
    }

    # 3 - Mark duplicates:
    # TODO: Picard Version is not the same!  Does it matter?
    call SRUTIL.MarkDuplicates as t_009_MarkDuplicates {
        input:
            input_bam = t_008_SortAlignedBam.output_bam,
            prefix = sample_name + ".aligned.sorted.marked_duplicates"
    }

#    # 4 - Reorder bam is BROKEN.
#    # USING SORT SAM.
#    call ReorderSam as t_010_ReorderSam {
#        input:
#            input_bam = t_009_MarkDuplicates.bam,
#            reference_fasta = ref_map["fasta"],
#            reference_fai = ref_map["fai"],
#            reference_dict = ref_map["dict"],
#            prefix = sample_name + ".aligned.sorted.marked_duplicates.reordered"
#    }

    call Utils.SortSam as t_010_ReorderSam {
        input:
            input_bam = t_009_MarkDuplicates.bam,
            prefix = sample_name + ".aligned.sorted.marked_duplicates.reordered"
    }

    # 5 - Realign indels:
    call RealignIndels as t_011_RealignIndels {
        input:
            input_bam = t_010_ReorderSam.output_bam,
            input_bai = t_010_ReorderSam.output_bam_index,
            reference_fasta = ref_map["fasta"],
            reference_fai = ref_map["fai"],
            reference_dict = ref_map["dict"],
            prefix = sample_name + ".aligned.sorted.marked_duplicates.reordered.indels_realigned"
    }

    # 5a - Index Realign indels bam:
    call Utils.Index as t_012_IndexRealignIndels {
        input:
            bam = t_011_RealignIndels.indel_realigned_bam
    }

    # 6 - BQSR:
    call BQSR as t_012_BQSR {
        input:
            input_bam = t_011_RealignIndels.indel_realigned_bam,
            input_bai = t_012_IndexRealignIndels.bai,
            reference_fasta = ref_map["fasta"],
            reference_fai = ref_map["fai"],
            reference_dict = ref_map["dict"],
            known_sites = [
                          resource_vcf_7g8_gb4,
                          resource_vcf_hb3_dd2,
                          resource_vcf_3d7_hb3
                ],
            known_sites_indices = [
                          resource_vcf_7g8_gb4_index,
                          resource_vcf_hb3_dd2_index,
                          resource_vcf_3d7_hb3_index
            ],
            prefix = sample_name + ".aligned.sorted.marked_duplicates.reordered.indels_realigned.bqsr"
    }

    # 6a - Index bam:
    call Utils.Index as t_013_IndexBQSR {
        input:
            bam = t_012_BQSR.bqsr_bam
    }

    # 7 - call variants with haplotypecaller:
    # both VCF and GVCF mode:
    call HaplotypeCaller as t_013_HaplotypeCaller {
        input:
            input_bam = t_012_BQSR.bqsr_bam,
            input_bai = t_013_IndexBQSR.bai,
            reference_fasta = ref_map["fasta"],
            reference_fai = ref_map["fai"],
            reference_dict = ref_map["dict"],
            prefix = sample_name + ".raw"
    }
    call HaplotypeCaller as t_014_HaplotypeCallerGvcfMode {
        input:
            input_bam = t_012_BQSR.bqsr_bam,
            input_bai = t_013_IndexBQSR.bai,
            reference_fasta = ref_map["fasta"],
            reference_fai = ref_map["fai"],
            reference_dict = ref_map["dict"],
            prefix = sample_name,
            gvcf_mode = true
    }

    # 8 - Recalibrate variants:
    call BroadOnPremMalariaPipelineTasks.VariantRecalibrator as t_015_VariantRecalibrator {
        input:
            input_vcf = t_013_HaplotypeCaller.vcf,
            reference_fasta = ref_map["fasta"],
            reference_fai = ref_map["fai"],
            reference_dict = ref_map["dict"],
            resource_vcf_7g8_gb4 = resource_vcf_7g8_gb4,
            resource_vcf_hb3_dd2 = resource_vcf_hb3_dd2,
            resource_vcf_3d7_hb3 = resource_vcf_3d7_hb3,
            resource_vcf_3d7_hb3_index = resource_vcf_3d7_hb3_index,
            resource_vcf_hb3_dd2_index = resource_vcf_hb3_dd2_index,
            resource_vcf_7g8_gb4_index = resource_vcf_7g8_gb4_index,
            prefix = sample_name
    }

    # 9 - sort, compress, and index final outputs:
    call BroadOnPremMalariaPipelineTasks.SortCompressIndexVcf as t_016_SortCompressIndexRawVCF {
        input:
            input_vcf = t_013_HaplotypeCaller.vcf
    }

    call BroadOnPremMalariaPipelineTasks.SortCompressIndexVcf as t_017_SortCompressIndexVCF {
        input:
            input_vcf = t_015_VariantRecalibrator.vcf
    }

    # 9 - sort, compress, and index final outputs:
    call BroadOnPremMalariaPipelineTasks.SortCompressIndexVcf as t_018_SortCompressIndexGVCF {
        input:
            input_vcf = t_014_HaplotypeCallerGvcfMode.vcf
    }

    output {
        File raw_vcf = t_016_SortCompressIndexRawVCF.vcf
        File raw_vcf_index = t_016_SortCompressIndexRawVCF.vcf_index
        File recalibrated_vcf = t_017_SortCompressIndexVCF.vcf
        File recalibrated_vcf_index = t_017_SortCompressIndexVCF.vcf_index
        File gvcf = t_018_SortCompressIndexGVCF.vcf
        File gvcf_index = t_018_SortCompressIndexGVCF.vcf_index
    }
}

task FilterOutTargetOrganismReads {

    parameter_meta {
        prefix: "Prefix for output files."

        fq_end1: "FASTQ file containing end 1 of reads."
        fq_end2: "FASTQ file containing end 2 of reads."

        contaminant_reference_fasta: "Reference FASTA file for the genome of the organism to remove."
        contaminant_reference_fai: "Reference FASTA index file for the genome of the organism to remove."
        contaminant_reference_bowtie_indices: "Bowtie2 indices for the given target genome."

        contaminant_organism_label: "Label for the name of the target organism."

        runtime_attr_override: "Optional override for runtime attributes."
    }

    input {
        String prefix

        File fq_end1
        File fq_end2

        File contaminant_reference_fasta
        File contaminant_reference_fai
        Array[File] contaminant_reference_bowtie_indices

        String contaminant_organism_label = "contaminant_organism"

        RuntimeAttr? runtime_attr_override
    }

    Int bowtie_index_size = ceil(size(contaminant_reference_bowtie_indices, "GB"))
    Int disk_size = 10 + 10 * (bowtie_index_size + ceil(size([fq_end1, fq_end2, contaminant_reference_fasta, contaminant_reference_fai], "GB")))

    String contaminant_organism = sub(contaminant_organism_label, " ", "_")

    command <<<
        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        if [[ ${np} -gt 2 ]] ; then
            let np=${np}-1
        fi

        # Move the bowtie2 index files to the same directory as the reference fasta file:
        ref_dir=$(dirname ~{contaminant_reference_fasta})
        while read f ; do 
            indx_dirname=$( dirname ${f} )
            if [[ "${indx_dirname}" != "${ref_dir}" ]] ; then
                mv ${f} ${ref_dir}
            fi
        done < ~{write_lines(contaminant_reference_bowtie_indices)}

        # Get the basename of the bowtie2 index file so we can reference it in the bowtie2 command:
        bowtie2_index_basename=$(echo "~{contaminant_reference_bowtie_indices[0]}" | grep bt2 | head -n1 | sed 's@\.[a-zA-Z0-9]*\.bt2@@')

        echo "Aligning to contaminant genome:"
        bowtie2 -x ${bowtie2_index_basename} \
            --threads ${np} \
            -1 ~{fq_end1} \
            -2 ~{fq_end2} | \
        samtools view -@$((np-1)) -bS - > ~{prefix}.mapAndUnmapped.~{contaminant_organism}.bam

        echo "Filtering out ~{contaminant_organism} reads:"
        samtools view -@$((np-1)) -b -f 12 -F 256 ~{prefix}.mapAndUnmapped.~{contaminant_organism}.bam > ~{prefix}.unmapped.bam

        echo "Sorting unmapped reads:"
        samtools sort -@$((np-1)) -n ~{prefix}.unmapped.bam -o ~{prefix}.unmapped.sorted.bam

        echo "Converting back to FASTQ:"
        samtools fastq -@$((np-1)) -n ~{prefix}.unmapped.sorted.bam \
            -1 ~{prefix}.contaminationRemoved.1.fq \
            -2 ~{prefix}.contaminationRemoved.2.fq
    >>>

    output {
        File fq1 = "~{prefix}.contaminationRemoved.1.fq"
        File fq2 = "~{prefix}.contaminationRemoved.2.fq"
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


task ReorderSam {
    input {
        File input_bam

        String prefix

        File reference_fasta
        File reference_fai
        File reference_dict

        RuntimeAttr? runtime_attr_override
    }

    Int compression_level = 2

    Int disk_size = 1 + 4*ceil(size([input_bam, reference_fasta], "GB"))

    # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
    # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
    # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

    command <<<
        set -euxo pipefail

        tot_mem_mb=$(free -m | grep '^Mem' | awk '{print $2}')
        java_memory_size_mb=$((tot_mem_mb-5120))

        java -Dsamjdk.compression_level=~{compression_level} -Xms${java_memory_size_mb}m -jar /usr/picard/picard.jar \
            ReorderSam \
                INPUT=~{input_bam} \
                OUTPUT=~{prefix}.bam \
                REFERENCE_SEQUENCE=~{reference_fasta} \
                SEQUENCE_DICTIONARY=~{reference_dict} \

        samtools index -@${np} ~{prefix}.bam
    >>>

    output {
        File bam = "~{prefix}.bam"
        File bai = "~{prefix}.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
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

task RealignIndels {
    input {
        String prefix

        File input_bam
        File input_bai

        File reference_fasta
        File reference_fai
        File reference_dict

        RuntimeAttr? runtime_attr_override
    }

    Int compression_level = 2

    Int disk_size = 1 + 4*ceil(size([input_bam, reference_fasta], "GB"))

    # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
    # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
    # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

    command <<<
        ################################
        # Standard Preamble

        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        if [[ ${np} -gt 2 ]] ; then
            let np=${np}-1
        fi

        tot_mem_mb=$(free -m | grep '^Mem' | awk '{print $2}')

        ################################

        java_memory_size_mb=$((tot_mem_mb-5120))

        # Create regions for indel realignment:
        java -Xmx${java_memory_size_mb}M -jar /usr/GenomeAnalysisTK.jar \
            -T RealignerTargetCreator \
                -nct 1 \
                -nt ${np} \
                -R ~{reference_fasta} \
                -I ~{input_bam} \
                -o ~{prefix}.interval_list

        java -Xmx${java_memory_size_mb}M -jar /usr/GenomeAnalysisTK.jar \
            -T IndelRealigner \
                -nct 1 \
                -nt 1 \
                -R ~{reference_fasta} \
                -I ~{input_bam} \
                -targetIntervals ~{prefix}.interval_list \
                -o ~{prefix}.bam
    >>>

    output {
        File indel_realigned_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "broadinstitute/gatk3:3.5-0"
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


task BQSR {
    input {
        String prefix

        File input_bam
        File input_bai

        File reference_fasta
        File reference_fai
        File reference_dict

        Array[File] known_sites
        Array[File] known_sites_indices

        RuntimeAttr? runtime_attr_override
    }

    Int compression_level = 2

    Int known_sites_size = ceil(size(known_sites, "GB")) + ceil(size(known_sites_indices, "GB"))
    Int disk_size = 10 + known_sites_size + 20*ceil(size([input_bam, input_bai, reference_fasta, reference_fai, reference_dict], "GB"))

    # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
    # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
    # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

    command <<<
        ################################
        # Standard Preamble

        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        if [[ ${np} -gt 2 ]] ; then
            let np=${np}-1
        fi

        tot_mem_mb=$(free -m | grep '^Mem' | awk '{print $2}')

        ################################

        # Must copy the input bam index to the same directory as the input bam.
        # We have to do this because the GATK will look for the index in the same directory as the input bam
        # and it may not be there due to us indexing it in a separate task.
        mv ~{input_bai} $(dirname ~{input_bam})/

        java_memory_size_mb=$((tot_mem_mb-5120))

        java -Xmx${java_memory_size_mb}M -jar /usr/GenomeAnalysisTK.jar \
            -T BaseRecalibrator \
                -nct 8 \
                -nt 1 \
                -R ~{reference_fasta} \
                -I ~{input_bam} \
                -knownSites ~{sep=" -knownSites " known_sites} \
                -o ~{prefix}_recal_report.grp

        java -Xmx${java_memory_size_mb}M -jar /usr/GenomeAnalysisTK.jar \
            -T PrintReads \
                -nct 8 \
                -nt 1 \
                -R ~{reference_fasta} \
                -I ~{input_bam} \
                -BQSR ~{prefix}_recal_report.grp \
                -o ~{prefix}.bam
    >>>

    output {
        File bqsr_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "broadinstitute/gatk3:3.5-0"
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


task HaplotypeCaller {
    input {
        String prefix

        File input_bam
        File input_bai

        File reference_fasta
        File reference_fai
        File reference_dict

        Boolean gvcf_mode = false

        RuntimeAttr? runtime_attr_override
    }

    Int compression_level = 2

    Int disk_size = 1 + 4*ceil(size([input_bam, reference_fasta], "GB"))

    String gvcf_arg = if (gvcf_mode) then "-ERC GVCF" else ""
    String out_suffix = if (gvcf_mode) then ".g.vcf" else ".vcf"

    command <<<
        ################################
        # Standard Preamble

        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        if [[ ${np} -gt 2 ]] ; then
            let np=${np}-1
        fi

        tot_mem_mb=$(free -m | grep '^Mem' | awk '{print $2}')

        ################################

        # Must copy the input bam index to the same directory as the input bam.
        # We have to do this because the GATK will look for the index in the same directory as the input bam
        # and it may not be there due to us indexing it in a separate task.
        mv ~{input_bai} $(dirname ~{input_bam})/

        java_memory_size_mb=$((tot_mem_mb-5120))

        java -Xmx${java_memory_size_mb}M -jar /usr/GenomeAnalysisTK.jar \
            -T HaplotypeCaller \
                -nt 1 \
                -R ~{reference_fasta} \
                -I ~{input_bam} \
                ~{gvcf_arg} \
                -ploidy 2 \
                --interval_padding 100 \
                -o ~{prefix}.~{out_suffix} \
                -variant_index_type LINEAR \
                -variant_index_parameter 128000
    >>>

    output {
        File vcf = "~{prefix}.~{out_suffix}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "broadinstitute/gatk3:3.5-0"
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

task Error {

    input {
        String message
    }
    command {
        set -euxo pipefail
        echo ~{message}
        exit 1
    }
    runtime {
        docker: "ubuntu:22.04"
        memory: "512 MB"
        disks: "local-disk 10 HDD"
        bootDiskSizeGb: "10"
        preemptible: 1
        cpu: 1
    }
}
