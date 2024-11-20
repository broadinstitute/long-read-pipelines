version 1.0

import "../../structs/Structs.wdl"
import "../Utility/Utils.wdl"
import "../Utility/SRUtils.wdl" as SRUTIL
import "../VariantCalling/SRJointGenotyping.wdl" as SRJOINT

workflow CallVariantsWithHaplotypeCaller {
    meta {
        author: "Jonn Smith"
        description: "A workflow for calling small variants with GATK HaplotypeCaller from an Illumina BAM file using the methods laid out by Niare et al. (https://doi.org/10.1186/s12936-023-04632-0)."
    }

    parameter_meta {
        bam: "Input bam file containing reads from which to call variants."
        bai: "Input bam index for `bam`."

        prefix: "Prefix to use for output files."
        sample_id: "ID of the sample being called."

        call_vars_on_mitochondria:  "If true, will call variants on the mitochondrial contig."

        mito_contig: "Name of the mitochondrial contig."
        contigs_names_to_ignore:  "Array of names of contigs to ignore for the purposes of reporting variants."
    }

    input {
        File bam
        File bai

        String prefix
        String sample_id

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File genotype_gvcfs_intervals

        Boolean call_vars_on_mitochondria = false

        String mito_contig = "chrM"
        Array[String] contigs_names_to_ignore = ["RANDOM_PLACEHOLDER_VALUE"]  ## Required for ignoring any filtering - this is kind of a hack - TODO: fix the task!
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

        call HaplotypeCaller_NIARE_GATK4_VCF as CallVariantsWithHC {
            input:
                input_bam = bam,
                input_bam_index = bai,
                prefix = prefix + "." + contig_for_small_var,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_fai,
                ref_dict = ref_dict,
                interval_list = contig_for_small_var
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

    # Collapse the GVCF into a regular VCF:
    call SRJOINT.GenotypeGVCFs as CollapseGVCFtoVCF {
        input:
            input_gvcf_data = MergeGVCFs.output_vcf,
            input_gvcf_index = MergeGVCFs.output_vcf_index,
            interval_list = genotype_gvcfs_intervals,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            ref_dict = ref_dict,
            prefix = prefix,
    }

    output {
        File output_gvcf = MergeGVCFs.output_vcf
        File output_gvcf_index = MergeGVCFs.output_vcf_index
        File output_vcf = CollapseGVCFtoVCF.output_vcf
        File output_vcf_index = CollapseGVCFtoVCF.output_vcf_index
        File bamout = MergeVariantCalledBamOuts.output_bam
        File bamout_index = IndexBamout.bai
    }
}

task HaplotypeCaller_NIARE_GATK4_VCF {
    meta {
        author: "Jonn Smith"
        notes: "Call variants with GATK 4 in accordance with Niare et al. (https://doi.org/10.1186/s12936-023-04632-0).  Reproducing methods laid out by to perform testing."
    }

    input {
        File input_bam
        File input_bam_index

        String prefix

        File ref_dict
        File ref_fasta
        File ref_fasta_index

        String interval_list

        RuntimeAttr? runtime_attr_override
    }

    String output_file_name = prefix + ".g.vcf.gz"

    Int disk_size = 2*ceil(size([ref_fasta, ref_fasta_index, ref_dict, input_bam], "GiB") + 50)

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
        java_memory_size_mb=$((available_memory_mb-1024))
        echo Total available memory: ${available_memory_mb} MB >&2
        echo Memory reserved for Java: ${java_memory_size_mb} MB >&2

#        gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
#            HaplotypeCaller \
#                -R ~{ref_fasta} \
#                -I ~{input_bam} \
#                -ERC GVCF \
#                -ploidy 2 \                                    # Default: 2
#                --native-pair-hmm-threads 16 \                 # Default: 4
#                -O ~{output_file_name} \
#                --assembly-region-padding 100 \                # Default: 100
#                --max-num-haplotypes-in-population 128 \       # Default: 128
#                --kmer-size 10 \                               # Default: 10, 25
#                --kmer-size 25 \                               # Default: 10, 25
#                --min-dangling-branch-length 4 \               # default: 4
#                --heterozygosity 0.0029 \
#                --indel-heterozygosity 0.0017 \
#                --min-assembly-region-size 100 \               # default: 50
#                -L ~{interval_list} \
#                -mbq 5 \                                       # default 10
#                -DF MappingQualityReadFilter \
#                --base-quality-score-threshold 12 \            # 18
#                -bamout ~{prefix}.bamout.bam

        gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
                -R ~{ref_fasta} \
                -I ~{input_bam} \
                -ERC GVCF \
                -ploidy 2 \
                --native-pair-hmm-threads 16 \
                -O ~{output_file_name} \
                --assembly-region-padding 100 \
                --max-num-haplotypes-in-population 128 \
                --kmer-size 10 \
                --kmer-size 25 \
                --min-dangling-branch-length 4 \
                --heterozygosity 0.0029 \
                --indel-heterozygosity 0.0017 \
                --min-assembly-region-size 100 \
                -L ~{interval_list} \
                -mbq 5 \
                -DF MappingQualityReadFilter \
                --base-quality-score-threshold 12 \
                -bamout ~{prefix}.bamout.bam

    >>>

    output {
        File output_vcf = "~{output_file_name}"
        File output_vcf_index = "~{output_file_name}.tbi"
        File bamout = "~{prefix}.bamout.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
       cpu_cores:          4,
       mem_gb:             32,
       disk_gb:            disk_size,
       boot_disk_gb:       25,
       preemptible_tries:  1,
       max_retries:        1,
       docker:             "us.gcr.io/broad-dsp-lrma/sr-malaria-niare-pipeline:0.0.1"
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

    meta {
        author: "Jonn Smith"
        notes: "Merge separate bamouts into a single file using `samtools merge`."
    }

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
        mthreads=$((np-1))
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
       boot_disk_gb:       25,
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

task GenomicsDbImport {
    meta {
        author: "Jonn Smith"
        notes: "Import variants into Genomics DB in accordance with Niare et al. (https://doi.org/10.1186/s12936-023-04632-0).  Reproducing methods laid out by to perform testing."
    }

    input {
        File sample_name_map

        File interval_list

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String prefix

        Int batch_size = 100

        RuntimeAttr? runtime_attr_override
    }

    Int ref_size = ceil(size([sample_name_map, interval_list, ref_fasta, ref_fasta_fai, ref_dict], "GB"))

    Int disk_size = 1 + 4*ref_size

    command <<<
        set -euxo pipefail
        # We need at least 1 GB of available memory outside of the Java heap in order to execute native code, thus, limit
        # Java's memory by the total memory minus 1 GB. We need to compute the total memory as it might differ from
        # memory_size_gb because of Cromwell's retry with more memory feature.
        # Note: In the future this should be done using Cromwell's ${MEM_SIZE} and ${MEM_UNIT} environment variables,
        #       which do not rely on the output format of the `free` command.

        available_memory_mb=$(free -m | awk '/^Mem/ {print $2}')
        java_memory_size_mb=$((available_memory_mb-1024))
        echo Total available memory: ${available_memory_mb} MB >&2
        echo Memory reserved for Java: ${java_memory_size_mb} MB >&2

        gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            GenomicsDBImport \
                --sample-name-map ~{sample_name_map} \
                --genomicsdb-workspace-path ~{prefix}.genomicsDB \
                --batch-size ~{batch_size} \
                -L ~{interval_list} \
                --genomicsdb-segment-size 8048576 \
                --genomicsdb-vcf-buffer-size 160384

        tar -cf ~{prefix}.genomicsDB.tar ~{prefix}.genomicsDB
    >>>

    output {
        File output_genomicsdb = "~{prefix}.genomicsDB.tar"
    }

    #########################
        RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/sr-malaria-niare-pipeline:0.0.1"
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

task GenotypeGVCFs {
    meta {
        author: "Jonn Smith"
        notes: "Genotype GVCFs with GATK 4 in accordance with Niare et al. (https://doi.org/10.1186/s12936-023-04632-0).  Reproducing methods laid out by to perform testing.  Minor adaptations due to infrastructure."
    }

    input {
        File input_gvcf_data
        File? input_gvcf_index  # Required if passing a VCF file.

        File interval_list

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String prefix

        Int batch_size = 100

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        input_gvcf_data: { help: "Either a single GVCF file or a GenomicsDB Tar file." }
        interval_list: {
            localization_optional: true
        }
    }

    Int ref_size = ceil(size([input_gvcf_data, input_gvcf_index, ref_fasta, ref_fasta_fai, ref_dict, interval_list], "GB"))

    Int disk_size = 1 + 4*ref_size

    command <<<
        set -euxo pipefail
        # We need at least 1 GB of available memory outside of the Java heap in order to execute native code, thus, limit
        # Java's memory by the total memory minus 1 GB. We need to compute the total memory as it might differ from
        # memory_size_gb because of Cromwell's retry with more memory feature.
        # Note: In the future this should be done using Cromwell's ${MEM_SIZE} and ${MEM_UNIT} environment variables,
        #       which do not rely on the output format of the `free` command.

        # We must determine if our input variants are in a genomicsdb file or in a VCF.
        # The easiest way is to see if the input is a .tar file:

        is_genomics_db=true
        filename=$(basename -- "~{input_gvcf_data}")
        extension="${filename##*.}"
        if [[ "${extension}" != "tar" ]] ; then
            is_genomics_db=false
        fi

        if $is_genomics_db ; then
            tar -xf ~{input_gvcf_data}
            INPUT_FILE="gendb://$(basename ~{input_gvcf_data} .tar)"
        else
            INPUT_FILE=~{input_gvcf_data}
        fi

        available_memory_mb=$(free -m | awk '/^Mem/ {print $2}')
        java_memory_size_mb=$((available_memory_mb-1024))
        echo Total available memory: ${available_memory_mb} MB >&2
        echo Memory reserved for Java: ${java_memory_size_mb} MB >&2

        gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            GenotypeGVCFs \
                --genomicsdb-use-bcf-codec true \
                -R ~{ref_fasta} \
                -V ${INPUT_FILE} \
                -L ~{interval_list} \
                --max-genotype-count 1024 \
                -O ~{prefix}.vcf.gz \
                -stand-call-conf 30

    >>>

    output {
        File output_vcf = "~{prefix}.vcf.gz"
        File output_vcf_index = "~{prefix}.vcf.gz.tbi"
    }

    #########################
        RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/sr-malaria-niare-pipeline:0.0.1"
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

task NormalizeVcfSplittingMultiallelics {
    meta {
        author: "Jonn Smith"
        notes: "Normalize multiallelic variants in accordance with Niare et al. (https://doi.org/10.1186/s12936-023-04632-0).  Reproducing methods laid out by to perform testing.  Minor adaptations due to infrastructure."
    }

    input {
        File input_vcf
        File input_vcf_index

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String prefix

        Int batch_size = 100

        RuntimeAttr? runtime_attr_override
    }

    Int ref_size = ceil(size([input_vcf, input_vcf_index, ref_fasta, ref_fasta_fai, ref_dict], "GB"))

    Int disk_size = 1 + 4*ref_size

    command <<<
        set -euxo pipefail
        # We need at least 1 GB of available memory outside of the Java heap in order to execute native code, thus, limit
        # Java's memory by the total memory minus 1 GB. We need to compute the total memory as it might differ from
        # memory_size_gb because of Cromwell's retry with more memory feature.
        # Note: In the future this should be done using Cromwell's ${MEM_SIZE} and ${MEM_UNIT} environment variables,
        #       which do not rely on the output format of the `free` command.


      bcftools norm -m-any ~{input_vcf} | \
        bcftools norm --check-ref -w -f ~{ref_fasta} | \
        bcftools annotate \
            -Ob \
            -x 'ID' \
            -I +'%CHROM:%POS:%POS:%REF:%ALT' | \
        bcftools view -i 'AC>0' -Oz -o ~{prefix}.vcf.gz

      tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File output_vcf = "~{prefix}.vcf.gz"
        File output_vcf_index = "~{prefix}.vcf.gz.tbi"
    }

    #########################
        RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/sr-malaria-niare-pipeline:0.0.1"
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

task VariantRecalibratorIndel {
    meta {
        author: "Jonn Smith"
        notes: "Run VariantRecalibrator on INDELs in accordance with Niare et al. (https://doi.org/10.1186/s12936-023-04632-0).  Reproducing methods laid out by to perform testing.  Minor adaptations due to infrastructure."
    }

    input {
        File input_vcf
        File input_vcf_index

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File sites_only_vcf
        File sites_only_vcf_index

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int ref_size = ceil(size([input_vcf, input_vcf_index, ref_fasta, ref_fasta_fai, ref_dict, sites_only_vcf, sites_only_vcf_index], "GB"))

    Int disk_size = 1 + 4*ref_size

    command <<<
        set -euxo pipefail
        # We need at least 1 GB of available memory outside of the Java heap in order to execute native code, thus, limit
        # Java's memory by the total memory minus 1 GB. We need to compute the total memory as it might differ from
        # memory_size_gb because of Cromwell's retry with more memory feature.
        # Note: In the future this should be done using Cromwell's ${MEM_SIZE} and ${MEM_UNIT} environment variables,
        #       which do not rely on the output format of the `free` command.

        available_memory_mb=$(free -m | awk '/^Mem/ {print $2}')
        java_memory_size_mb=$((available_memory_mb-1024))
        echo Total available memory: ${available_memory_mb} MB >&2
        echo Memory reserved for Java: ${java_memory_size_mb} MB >&2

        gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            VariantRecalibrator \
                -R ~{ref_fasta} \
                -V ~{input_vcf} \
                --trust-all-polymorphic \
                -an QD -an DP -an FS -an SOR -an MQ \
                -mode INDEL \
                --max-gaussians 4 \
                -resource:Brown,known=true,training=true,truth=true,prior=15.0 ~{sites_only_vcf} \
                -O ~{prefix}.indel_recal \
                --output-model ~{prefix}.indel.model.report \
                --tranches-file  ~{prefix}.raw.indel.tranches \
#                --rscript-file  ~{prefix}.raw.indel.plots.R

    >>>

    output {
        File recalibration = "~{prefix}.indel_recal"
        File recalibration_index = "~{prefix}.indel_recal.idx"
        File tranches = "~{prefix}.raw.indel.tranches"
        File model_report = "~{prefix}.indel.model.report"
    }

    #########################
        RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/sr-malaria-niare-pipeline:0.0.1"
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

task VariantRecalibratorSnp {
    meta {
        author: "Jonn Smith"
        notes: "Run VariantRecalibrator on SNPs in accordance with Niare et al. (https://doi.org/10.1186/s12936-023-04632-0).  Reproducing methods laid out by to perform testing.  Minor adaptations due to infrastructure."
    }

    input {
        File input_vcf
        File input_vcf_index

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File sites_only_vcf
        File sites_only_vcf_index

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int ref_size = ceil(size([input_vcf, input_vcf_index, ref_fasta, ref_fasta_fai, ref_dict, sites_only_vcf, sites_only_vcf_index], "GB"))

    Int disk_size = 1 + 4*ref_size

    command <<<
        set -euxo pipefail
        # We need at least 1 GB of available memory outside of the Java heap in order to execute native code, thus, limit
        # Java's memory by the total memory minus 1 GB. We need to compute the total memory as it might differ from
        # memory_size_gb because of Cromwell's retry with more memory feature.
        # Note: In the future this should be done using Cromwell's ${MEM_SIZE} and ${MEM_UNIT} environment variables,
        #       which do not rely on the output format of the `free` command.

        available_memory_mb=$(free -m | awk '/^Mem/ {print $2}')
        java_memory_size_mb=$((available_memory_mb-1024))
        echo Total available memory: ${available_memory_mb} MB >&2
        echo Memory reserved for Java: ${java_memory_size_mb} MB >&2

        gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            VariantRecalibrator \
                -R ~{ref_fasta} \
                -V ~{input_vcf} \
                --trust-all-polymorphic \
                -an QD -an DP -an FS -an SOR -an MQ \
                -mode SNP \
                --max-gaussians 4 \
                -resource:Brown,known=true,training=true,truth=true,prior=15.0 ~{sites_only_vcf} \
                -O ~{prefix}.snp_recal \
                --tranches-file  ~{prefix}.raw.snp.tranches \
                --output-model ~{prefix}.snp.model.report \
#                --rscript-file  ~{prefix}.raw.snp.plots.R

    >>>

    output {
        File recalibration = "~{prefix}.snp_recal"
        File recalibration_index = "~{prefix}.snp_recal.idx"
        File tranches = "~{prefix}.raw.snp.tranches"
        File model_report = "~{prefix}.snp.model.report"
    }

    #########################
        RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/sr-malaria-niare-pipeline:0.0.1"
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

task ApplyVqsrIndel {
    meta {
        author: "Jonn Smith"
        notes: "Run GATK 4 ApplyVqsr for INDELs in accordance with Niare et al. (https://doi.org/10.1186/s12936-023-04632-0).  Reproducing methods laid out by to perform testing.  Minor adaptations due to infrastructure."
    }

    input {
        File input_vcf
        File input_vcf_index

        File recal_file
        File recal_file_index
        File recal_tranches

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int ref_size = ceil(size(input_vcf, "GB") + size(recal_file, "GB") + size(recal_tranches, "GB"))

    Int disk_size = 1 + 4*ref_size

    command <<<
        set -euxo pipefail
        # We need at least 1 GB of available memory outside of the Java heap in order to execute native code, thus, limit
        # Java's memory by the total memory minus 1 GB. We need to compute the total memory as it might differ from
        # memory_size_gb because of Cromwell's retry with more memory feature.
        # Note: In the future this should be done using Cromwell's ${MEM_SIZE} and ${MEM_UNIT} environment variables,
        #       which do not rely on the output format of the `free` command.

        available_memory_mb=$(free -m | awk '/^Mem/ {print $2}')
        java_memory_size_mb=$((available_memory_mb-1024))
        echo Total available memory: ${available_memory_mb} MB >&2
        echo Memory reserved for Java: ${java_memory_size_mb} MB >&2

        gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
           ApplyVQSR \
            -V ~{input_vcf} \
            --recal-file ~{recal_file} \
            --tranches-file  ~{recal_tranches} \
            --create-output-variant-index true \
            --lod-score-cutoff -2.0 \
            --exclude-filtered false \
            -mode INDEL \
            -O ~{prefix}.indel_recal.vcf.gz

    >>>

    output {
        File output_vcf = "~{prefix}.indel_recal.vcf.gz"
        File output_vcf_index = "~{prefix}.indel_recal.vcf.gz.tbi"
    }

    #########################
        RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/sr-malaria-niare-pipeline:0.0.1"
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

task ApplyVqsrSnp {
    meta {
        author: "Jonn Smith"
        notes: "Run GATK 4 ApplyVQSR for SNPs in accordance with Niare et al. (https://doi.org/10.1186/s12936-023-04632-0).  Reproducing methods laid out by to perform testing.  Minor adaptations due to infrastructure."
    }

    input {
        File input_vcf
        File input_vcf_index

        File recal_file
        File recal_file_index
        File recal_tranches

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int ref_size = ceil(size(input_vcf, "GB") + size(recal_file, "GB") + size(recal_tranches, "GB"))

    Int disk_size = 1 + 4*ref_size

    command <<<
        set -euxo pipefail
        # We need at least 1 GB of available memory outside of the Java heap in order to execute native code, thus, limit
        # Java's memory by the total memory minus 1 GB. We need to compute the total memory as it might differ from
        # memory_size_gb because of Cromwell's retry with more memory feature.
        # Note: In the future this should be done using Cromwell's ${MEM_SIZE} and ${MEM_UNIT} environment variables,
        #       which do not rely on the output format of the `free` command.

        available_memory_mb=$(free -m | awk '/^Mem/ {print $2}')
        java_memory_size_mb=$((available_memory_mb-1024))
        echo Total available memory: ${available_memory_mb} MB >&2
        echo Memory reserved for Java: ${java_memory_size_mb} MB >&2

        gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
           ApplyVQSR \
            -V ~{input_vcf} \
            --recal-file ~{recal_file} \
            --tranches-file  ~{recal_tranches} \
            --create-output-variant-index true \
            --lod-score-cutoff 0.0 \
            --exclude-filtered false \
            -mode SNP \
            -O ~{prefix}.snp_recal.vcf.gz

    >>>

    output {
        File output_vcf = "~{prefix}.snp_recal.vcf.gz"
        File output_vcf_index = "~{prefix}.snp_recal.vcf.gz.tbi"
    }

    #########################
        RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/sr-malaria-niare-pipeline:0.0.1"
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

task MergeMultiAllelicSitesPostRecalibration {
    meta {
        author: "Jonn Smith"
        notes: "Merge multi-allelic sites in accordance with Niare et al. (https://doi.org/10.1186/s12936-023-04632-0).  Reproducing methods laid out by to perform testing.  Minor adaptations due to infrastructure."
    }

    input {
        File input_vcf
        File input_vcf_index

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int ref_size = ceil(size([input_vcf, input_vcf_index, ref_fasta, ref_fasta_fai, ref_dict], "GB"))

    Int disk_size = 1 + 4*ref_size

    command <<<
        set -euxo pipefail
        # We need at least 1 GB of available memory outside of the Java heap in order to execute native code, thus, limit
        # Java's memory by the total memory minus 1 GB. We need to compute the total memory as it might differ from
        # memory_size_gb because of Cromwell's retry with more memory feature.
        # Note: In the future this should be done using Cromwell's ${MEM_SIZE} and ${MEM_UNIT} environment variables,
        #       which do not rely on the output format of the `free` command.

        # Make sure we use all our processors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')

        available_memory_mb=$(free -m | awk '/^Mem/ {print $2}')
        java_memory_size_mb=$((available_memory_mb-1024))
        echo Total available memory: ${available_memory_mb} MB >&2
        echo Memory reserved for Java: ${java_memory_size_mb} MB >&2

        gatk --java-options "-Xmx${java_memory_size_mb}m -Xms${java_memory_size_mb}m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            SelectVariants \
                -R ~{ref_fasta} \
                -V ~{input_vcf} \
                -O ~{prefix}.raw.recal.pass.vcf.gz \
                --exclude-filtered true

        bcftools norm \
            -m+any ~{prefix}.raw.recal.pass.vcf.gz \
            --check-ref -f ~{ref_fasta} \
            -Oz \
            -o ~{prefix}.pass.merged.vcf.gz \
            --threads ${np}

        tabix -p vcf ~{prefix}.pass.merged.vcf.gz

    >>>

    output {
        File output_vcf = "~{prefix}.pass.merged.vcf.gz"
        File output_vcf_index = "~{prefix}.pass.merged.vcf.gz.tbi"
    }

    #########################
        RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/sr-malaria-niare-pipeline:0.0.1"
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

