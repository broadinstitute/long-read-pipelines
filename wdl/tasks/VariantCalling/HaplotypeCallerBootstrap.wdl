version 1.0

import "../../structs/Structs.wdl"
import "../Utility/Utils.wdl"
import "../Utility/SRUtils.wdl" as SRUTIL
import "../VariantCalling/SRJointGenotyping.wdl" as SRJOINT

workflow CallVariantsWithHaplotypeCaller {
    meta {
        author: "Jonn Smith, adapted by Raphael Brosula"
        description: "A workflow for calling small variants with GATK HaplotypeCaller from an Illumina BAM file, with optional specification of a SNP VCF."
    }

    input {
        File bam
        File bai

        String prefix
        String sample_id

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File? dbsnp_vcf

        Boolean call_vars_on_mitochondria = true

        Int ploidy = 2

        Float heterozygosity = 0.001
        Float heterozygosity_stdev = 0.01
        Float indel_heterozygosity = 0.000125

        Float max_reads_per_alignment_start
        Float max_num_haplotypes_in_population

        Boolean enable_pileup_mode = false

        String mito_contig = "chrM"
        Array[String] contigs_names_to_ignore = ["RANDOM_PLACEHOLDER_VALUE"]  ## Required for ignoring any filtering - this is kind of a hack - TODO: fix the task!

        RuntimeAttr? runtime_attr_override
    }

    # Scatter by chromosome:
    Array[String] use_filter = if (call_vars_on_mitochondria) then contigs_names_to_ignore else flatten([[mito_contig], contigs_names_to_ignore])
    call Utils.MakeChrIntervalList as SmallVariantsScatterPrep {
        input:
            ref_dict = ref_dict,
            filter = use_filter
    }

    # Call over the scattered intervals:
    scatter (c in SmallVariantsScatterPrep.chrs) {
        String contig_for_small_var = c[0]

        call HaplotypeCaller_GATK4_VCF as CallVariantsWithHC {
            input:
                input_bam = bam,
                input_bam_index = bai,
                prefix = prefix + "." + contig_for_small_var,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_fai,
                ref_dict = ref_dict,
                make_gvcf = true,
                make_bamout = true,
                enable_pileup_mode = enable_pileup_mode,
                single_interval = contig_for_small_var,
                contamination = 0,
                ploidy = ploidy,
                heterozygosity = heterozygosity,
                heterozygosity_stdev = heterozygosity_stdev,
                indel_heterozygosity = indel_heterozygosity,
                max_reads_per_alignment_start = max_reads_per_alignment_start,
                max_num_haplotypes_in_population = max_num_haplotypes_in_population,
                use_spanning_event_genotyping = true,
                runtime_attr_override = runtime_attr_override
        }
    }

    # Merge the output GVCFs:
    call SRUTIL.MergeVCFs as MergeGVCFs {
        input:
            input_vcfs = CallVariantsWithHC.output_vcf,
            input_vcfs_indexes = CallVariantsWithHC.output_vcf_index,
            prefix = prefix
    }

    # Merge the output BAMs:
    call MergeBamouts as MergeVariantCalledBamOuts {
        input:
            bams = CallVariantsWithHC.bamout,
            prefix = "~{prefix}.bamout"
    }

    # Index the Bamout:
    call Utils.Index as IndexBamout {
        input:
            bam = MergeVariantCalledBamOuts.output_bam
    }

   # Now reblock the GVCF to combine hom ref blocks and save $ / storage:
   call ReblockGVCF as ReblockHcGVCF {
        input:
            gvcf = MergeGVCFs.output_vcf,
            gvcf_index = MergeGVCFs.output_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_dict = ref_dict,
            prefix = prefix
    }

    # Collapse the GVCF into a regular VCF:
    call SRJOINT.GenotypeGVCFs as CollapseGVCFtoVCF {
        input:
            input_gvcf_data = ReblockHcGVCF.output_gvcf,
            input_gvcf_index = ReblockHcGVCF.output_gvcf_index,
            interval_list = SmallVariantsScatterPrep.interval_list,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_dict = ref_dict,
            prefix = prefix,
            dbsnp_vcf = dbsnp_vcf 
    }

    output {
        File output_gvcf = ReblockHcGVCF.output_gvcf
        File output_gvcf_index = ReblockHcGVCF.output_gvcf_index
        File output_vcf = CollapseGVCFtoVCF.output_vcf
        File output_vcf_index = CollapseGVCFtoVCF.output_vcf_index
        File bamout = MergeVariantCalledBamOuts.output_bam
        File bamout_index = IndexBamout.bai
    }
}

task HaplotypeCaller_GATK4_VCF {
    meta {
        author: "Jonn Smith"
        notes: "Adapted from the WARP pipeline found here: https://github.com/broadinstitute/warp.git"
    }

    input {
        File input_bam
        File input_bam_index

        String prefix

        File ref_dict
        File ref_fasta
        File ref_fasta_index

        Int ploidy = 2

        Float heterozygosity = 0.001
        Float heterozygosity_stdev = 0.01
        Float indel_heterozygosity = 0.000125

        Float max_reads_per_alignment_start = 50
        Float max_num_haplotypes_in_population = 128

        Boolean make_gvcf
        Boolean make_bamout

        String? single_interval
        File? interval_list
        Float? contamination

        Boolean use_spanning_event_genotyping = true

        Boolean enable_pileup_mode = false

        RuntimeAttr? runtime_attr_override
    }

    String output_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
    String output_file_name = prefix + output_suffix

    Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
    Int disk_size = 2*ceil(((size(input_bam, "GiB") + 30)) + ref_size) + 30

    String bamout_arg = if make_bamout then "-bamout ~{prefix}.bamout.bam" else ""

    String interval_arg = if (defined(interval_list) || defined(single_interval)) then " -L " else ""
    String interval_arg_value = if defined(interval_list) then select_first([interval_list]) else if defined(single_interval) then select_first([single_interval]) else ""

    parameter_meta {
        input_bam: { localization_optional: true }
    }

    command <<<
        set -euxo pipefail

        # We need at least 1 GB of available memory outside of the Java heap in order to execute native code, thus, limit
        # Java's memory by the total memory minus 1 GB. We need to compute the total memory as it might differ from
        # memory_size_gb because of Cromwell's retry with more memory feature.
        # Note: In the future this should be done using Cromwell's ${MEM_SIZE} and ${MEM_UNIT} environment variables,
        #       which do not rely on the output format of the `free` command.

        available_memory_mb=$(free -m | awk '/^Mem/ {print $2}')
        let java_memory_size_mb=available_memory_mb-1024
        echo Total available memory: ${available_memory_mb} MB >&2
        echo Memory reserved for Java: ${java_memory_size_mb} MB >&2

        gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
                -R ~{ref_fasta} \
                -I ~{input_bam} \
                ~{interval_arg}~{default="" interval_arg_value} \
                -O ~{output_file_name} \
                -contamination ~{default=0 contamination} \
                --sample-ploidy ~{ploidy} \
                --heterozygosity ~{heterozygosity} \
                --heterozygosity-stdev ~{heterozygosity_stdev} \
                --indel-heterozygosity ~{indel_heterozygosity} \
                --max-reads-per-alignment-start ~{max_reads_per_alignment_start} \
                --max-num-haplotypes-in-population ~{max_num_haplotypes_in_population} \
                --linked-de-bruijn-graph \
                ~{true="--pileup-detection --pileup-detection-enable-indel-pileup-calling" false="" enable_pileup_mode} \
                --annotate-with-num-discovered-alleles \
                -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
                ~{false="--disable-spanning-event-genotyping" true="" use_spanning_event_genotyping} \
                -G StandardAnnotation -G StandardHCAnnotation  \
                -A AssemblyComplexity \
                ~{true="-ERC GVCF" false="" make_gvcf} \
                --smith-waterman FASTEST_AVAILABLE \
                ~{bamout_arg} 

        # Removed for now because we need to qualify the pipeline with standard annotations first.
        # ~{true="-G AS_StandardAnnotation" false="" make_gvcf}

        # Cromwell doesn't like optional task outputs, so we have to touch this file.
        touch ~{prefix}.bamout.bam
    >>>

    output {
        File output_vcf = "~{output_file_name}"
        File output_vcf_index = "~{output_file_name}.tbi"
        File bamout = "~{prefix}.bamout.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
       cpu_cores:          2,
       mem_gb:             24,
       disk_gb:            disk_size,
       boot_disk_gb:       23,
       preemptible_tries:  1,
       max_retries:        1,
       docker:             "broadinstitute/gatk-nightly:2024-05-02-4.5.0.0-27-gec39c37d9-NIGHTLY-SNAPSHOT"
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


# This task is here because merging bamout files using Picard produces an error.
task MergeBamouts {

    input {
        Array[File] bams
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(bams, "GiB") * 2) + 10

    command <<<

        set -euxo pipefail

        # Make sure we use all our processors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        ithreads=${np}

        # If the number of processors = 1, then `let` will return 1 here:
        # So we need to turn off `set -e` for this command:
        set +e
        let mthreads=${np}-1
        set -e

        samtools merge -@${mthreads} ~{prefix}.bam ~{sep=" " bams}
        samtools index -@${ithreads} ~{prefix}.bam
        mv ~{prefix}.bam.bai ~{prefix}.bai
    >>>

    #########################
    RuntimeAttr default_attr = object {
       cpu_cores:          1,
       mem_gb:             4,
       disk_gb:            disk_size,
       boot_disk_gb:       10,
       preemptible_tries:  1,
       max_retries:        1,
       docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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
        File output_bam = "~{prefix}.bam"
        File output_bam_index = "~{prefix}.bai"
    }
}

task ReblockGVCF {

    input {
        File gvcf
        File gvcf_index

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String prefix

        Array[Int] gq_blocks = [20, 30, 40]

        Float? tree_score_cutoff

        Array[String]? annotations_to_keep

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil((size(gvcf, "GiB") * 4) + size(ref_fasta, "GiB") + size(ref_fasta_fai, "GiB") + size(ref_dict, "GiB") + 10)

    String annotations_to_keep_arg = if defined(annotations_to_keep) then "--annotations-to-keep" else ""

    command {
        set -euxo pipefail

        gatk --java-options "-Xms3000m -Xmx3000m" \
            ReblockGVCF \
                -R ~{ref_fasta} \
                -V ~{gvcf} \
                -do-qual-approx \
                -A AssemblyComplexity \
                --annotate-with-num-discovered-alleles \
                --floor-blocks \
                -GQB ~{sep=" -GQB " gq_blocks} \
                ~{"--tree-score-threshold-to-no-call " + tree_score_cutoff} \
                ~{annotations_to_keep_arg} ~{sep=" --annotations-to-keep " annotations_to_keep} \
                -O ~{prefix}.rb.g.vcf.gz
    }

    #########################
    RuntimeAttr default_attr = object {
       cpu_cores:          2,
       mem_gb:             4,
       disk_gb:            disk_size,
       boot_disk_gb:       15,
       preemptible_tries:  1,
       max_retries:        1,
       docker:             "broadinstitute/gatk-nightly:2024-05-02-4.5.0.0-27-gec39c37d9-NIGHTLY-SNAPSHOT"
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
        File output_gvcf = "~{prefix}.rb.g.vcf.gz"
        File output_gvcf_index = "~{prefix}.rb.g.vcf.gz.tbi"
    }
}
