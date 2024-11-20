version 1.0

##########################################################################################
# Calls (small) variants on an input LR BAM with GATK4 HaplotypeCaller pipeline.

# This is essentially copying "dsde_pipelines_tasks/VariantCalling.wdl"
# with some customization to fit the process described in
# "https://github.com/PacificBiosciences/hg002-ccs/"
##########################################################################################

import "dsde_pipelines_tasks/GermlineVariantDiscovery.wdl" as Calling
import "dsde_pipelines_tasks/Qc.wdl" as QC
import "dsde_pipelines_tasks/Utilities.wdl" as DSDEPipelinesUtils

import "../../tasks/Utility/Utils.wdl" as Utils

workflow GATKBestPraciceForLR {
    input {
        File calling_interval_list
        Int calling_intervals_scatter_count
        Int calling_intervals_break_bands_at_multiples_of

        Float? contamination

        File input_bam

        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File par_regions_bed

        Boolean? sample_is_female

        String output_prefix

        Boolean run_qc_on_variants
        Boolean make_bamout = false
        File? var_calling_metrics_eval_interval_list

        File dbsnp_vcf
        File dbsnp_vcf_index

        Int agg_preemptible_tries = 1 # can not be optional as the tasks and sub-workflows copied over has this non-optional
    }

    parameter_meta {
        calling_interval_list :                         "an interval list file holding intervals on which to call variants"
        calling_intervals_scatter_count :               "a run-time performance tunning integer indicating how many scatters to have on the intervals"
        calling_intervals_break_bands_at_multiples_of : "if > 0, reference bands will be broken up at genomic positions that are multiples of this number"

        contamination :    "[optional] a floating number in range [0,1] indicating fraction of contamination in sequencing data (for all samples) to aggressively remove"
        input_bam :        "BAM on which to call variants"
        sample_is_female : "boolean indicating if the sample is a female or not"

        ref_fasta :       "reference to which the BAM was aligned to"
        ref_fasta_index : "index accompanying the reference"
        ref_dict :        "dictionary accompanying the reference"
        par_regions_bed : "BED file holding PAR region of the reference genome"

        output_prefix : "prefix to output files"

        run_qc_on_variants : "boolean indicating to run QC on the output variants or not"
        dbsnp_vcf :          "DBSNP vcf for quality control purposes"
        dbsnp_vcf_index :    "index accompanying the DBSNP vcf"
        var_calling_metrics_eval_interval_list : "[optional] an interval_list file telling us where to run variant QC"

        make_bamout: "[default-valued] For CNNScoreVariants to run with a 2D model, a bamout must be created by HaplotypeCaller. The bamout is a bam containing information on how HaplotypeCaller remapped reads while it was calling variants. See https://gatkforums.broadinstitute.org/gatk/discussion/5484/howto-generate-a-bamout-file-showing-how-haplotypecaller-has-remapped-sequence-reads for more details."
    }

    ###########################################################################
    # Break the calling interval_list into sub-intervals
    # Perform variant calling on the sub-intervals, and then gather the results
    call DSDEPipelinesUtils.ScatterIntervalList as ScatterIntervalList {
        input:
            interval_list = calling_interval_list,
            scatter_count = calling_intervals_scatter_count,
            break_bands_at_multiples_of = calling_intervals_break_bands_at_multiples_of
    }

    ###########################################################################
    # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
    # If we take the number we are scattering by and reduce by 20 we will have enough disk space
    # to account for the fact that the data is quite uneven across the shards.
    Int potential_hc_divisor = ScatterIntervalList.interval_count - 20
    Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

    String dbsnp_vcf_gspath = dbsnp_vcf

    # Call variants in parallel over WGS calling intervals
    scatter (scattered_interval_list in ScatterIntervalList.out) {

        # Generate GVCF by interval
        call Calling.HaplotypeCaller_GATK4_VCF as HaplotypeCallerGATK4 {
            input:
                contamination = contamination,
                input_bam = input_bam,
                interval_list = scattered_interval_list,
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                par_regions_bed = par_regions_bed,
                hc_scatter = hc_divisor,
                make_gvcf = true, # always produce gVCF
                make_bamout = make_bamout,
                preemptible_tries = agg_preemptible_tries,

                sample_gender_known_as_female = sample_is_female,

                vcf_basename = output_prefix + ".scatter.hc",
        }

        call GenotypeGVCFs as Genotyping {
            input:
                input_gvcf = HaplotypeCallerGATK4.output_vcf,
                input_gvcf_index = HaplotypeCallerGATK4.output_vcf_index,
                interval = scattered_interval_list,
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,

                dbsnp_vcf_gspath = dbsnp_vcf_gspath
        }

        # If bamout files were created, we need to sort and gather them into one bamout
        if (make_bamout) {
            call Utils.SortSam as SortBamout {
                input:
                    input_bam = HaplotypeCallerGATK4.bamout,
                    prefix = output_prefix,
            }
        }

        File gvcfs_to_merge = HaplotypeCallerGATK4.output_vcf
        File gvcf_indices_to_merge = HaplotypeCallerGATK4.output_vcf_index

        File vcfs_to_merge = Genotyping.output_vcf
        File vcf_indices_to_merge = Genotyping.output_vcf_index
    }

    ###########################################################################
    # Combine by-interval (g)VCFs into a single sample (g)VCF file
    call Calling.MergeVCFs as MergeGVCFs {
        input:
            input_vcfs = gvcfs_to_merge,
            input_vcfs_indexes = gvcf_indices_to_merge,
            output_vcf_name = output_prefix + ".haplotypecaller.g.vcf.gz",
            preemptible_tries = agg_preemptible_tries
    }

    call Calling.MergeVCFs as MergeVCFs {
        input:
            input_vcfs = vcfs_to_merge,
            input_vcfs_indexes = vcf_indices_to_merge,
            output_vcf_name = output_prefix + ".vcf.gz",
            preemptible_tries = agg_preemptible_tries
    }

    ###########################################################################
    if (make_bamout) {
        call MergeBamouts {
            input:
                bams = select_all(SortBamout.output_bam),
                output_base_name = output_prefix
        }
    }

    ###########################################################################
    if (run_qc_on_variants) {
        # Validate the (g)VCF output of HaplotypeCaller
        call QC.ValidateVCF as ValidateGVCF {
            input:
                input_vcf = MergeGVCFs.output_vcf,
                input_vcf_index = MergeGVCFs.output_vcf_index,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                calling_interval_list = calling_interval_list,
                is_gvcf = true,
                dbsnp_vcf = dbsnp_vcf,
                dbsnp_vcf_index = dbsnp_vcf_index,
                gatk_docker = "us.gcr.io/broad-gatk/gatk:latest",
                preemptible_tries = agg_preemptible_tries
        }
        call QC.ValidateVCF as ValidateVCF {
            input:
                input_vcf = MergeVCFs.output_vcf,
                input_vcf_index = MergeVCFs.output_vcf_index,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                ref_dict = ref_dict,
                calling_interval_list = calling_interval_list,
                is_gvcf = false,
                dbsnp_vcf = dbsnp_vcf,
                dbsnp_vcf_index = dbsnp_vcf_index,
                gatk_docker = "us.gcr.io/broad-gatk/gatk:latest",
                preemptible_tries = agg_preemptible_tries
        }

        # QC the (g)VCF
        call QC.CollectVariantCallingMetrics as CollectVariantCallingMetricsGVCF {
            input:
                input_vcf = MergeGVCFs.output_vcf,
                input_vcf_index = MergeGVCFs.output_vcf_index,
                metrics_basename = output_prefix,
                ref_dict = ref_dict,
                is_gvcf = true,
                dbsnp_vcf = dbsnp_vcf,
                dbsnp_vcf_index = dbsnp_vcf_index,
                evaluation_interval_list = select_first([var_calling_metrics_eval_interval_list, ScatterIntervalList.out[2]]),
                preemptible_tries = agg_preemptible_tries
        }
        call QC.CollectVariantCallingMetrics as CollectVariantCallingMetrics {
            input:
                input_vcf = MergeVCFs.output_vcf,
                input_vcf_index = MergeVCFs.output_vcf_index,
                metrics_basename = output_prefix,
                ref_dict = ref_dict,
                is_gvcf = false,
                dbsnp_vcf = dbsnp_vcf,
                dbsnp_vcf_index = dbsnp_vcf_index,
                evaluation_interval_list = select_first([var_calling_metrics_eval_interval_list, ScatterIntervalList.out[2]]),
                preemptible_tries = agg_preemptible_tries
        }
    }

    ###########################################################################
    call PostProcess as PostProcess {
        input:
            input_vcf = MergeVCFs.output_vcf,
            input_vcf_index = MergeVCFs.output_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict
    }

    output {
        File out_pp_vcf = PostProcess.output_vcf
        File out_pp_vcf_index = PostProcess.output_vcf_index

        File output_vcf = MergeVCFs.output_vcf
        File output_vcf_index = MergeVCFs.output_vcf_index

        File output_gvcf = MergeGVCFs.output_vcf
        File output_gvcf_index = MergeGVCFs.output_vcf_index

        File? vcf_summary_metrics = CollectVariantCallingMetrics.summary_metrics
        File? vcf_detail_metrics = CollectVariantCallingMetrics.detail_metrics

        File? gvcf_summary_metrics = CollectVariantCallingMetricsGVCF.summary_metrics
        File? gvcf_detail_metrics = CollectVariantCallingMetricsGVCF.detail_metrics

        File? bamout = MergeBamouts.output_bam
        File? bamout_index = MergeBamouts.output_bam_index
    }
}

# This task is here because merging bamout files using Picard produces an error.
task MergeBamouts {
    input {
        Array[File] bams
        String output_base_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(bams, "GiB") * 2) + 10

    command {
        samtools merge ~{output_base_name}.bam ~{sep=" " bams}
        samtools index ~{output_base_name}.bam
        mv ~{output_base_name}.bam.bai ~{output_base_name}.bai
    }

    output {
        File output_bam = "~{output_base_name}.bam"
        File output_bam_index = "~{output_base_name}.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsde-methods/samtools-cloud:v1.clean"
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

# Copied from
# https://github.com/broadinstitute/dsde-pipelines/blob/1275cdb34f6589ea51d54daa33b3d948ac9da316/genomes_in_the_cloud/joint_genotyping_workflow/JointGenotypingWf.wdl
# then did some modification to suite
# https://github.com/PacificBiosciences/hg002-ccs/blob/master/smallvariants/gatk_genotypeGVCFs.sh
# Note that the original task depends on GenomicsDB for large scale joint genotyping, we don't have that concern here for a single sample pipeline.
# The original task was using "String" for what should be "File"s, again because of the large cohort JG hence using the NIO cloud feature.
# Here we revert it back to File to take advantage of call-caching considering that we don't have a large cohort.
task GenotypeGVCFs {
    input {
        File input_gvcf
        File input_gvcf_index

        File interval

        File ref_dict
        File ref_fasta
        File ref_fasta_index

        String dbsnp_vcf_gspath # using NIO to cut localization time because it is expected to be large

        Float std_call_conf = 2.0 # this default value follows the CCS paper custom script value

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(input_gvcf, "GiB") + size(ref_fasta, "GiB")) + 10
    String out_vcf = basename(input_gvcf, ".g.vcf.gz") + ".vcf.gz"

    command <<<
        set -euo pipefail

        # compared to the original, we take out "--only-output-calls-starting-in-intervals"
        # because we don't have to deal with a large (close to TB) VCF
        gatk --java-options -Xms5g GenotypeGVCFs \
            -V ~{input_gvcf} \
            -R ~{ref_fasta} \
            -O ~{out_vcf} \
            -D ~{dbsnp_vcf_gspath} \
            -L ~{interval} \
            --annotation-group StandardAnnotation \
            --annotation-group AS_StandardAnnotation \
            --annotation-group StandardHCAnnotation \
            --standard-min-confidence-threshold-for-calling ~{std_call_conf}
    >>>

    output {
        File output_vcf = "~{out_vcf}"
        File output_vcf_index = "~{out_vcf}.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             7,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-gatk/gatk:latest"
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

# Postprocessing the VCFs following the CCS paper
task PostProcess {
    input {
        File input_vcf
        File input_vcf_index

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(input_vcf, "GiB")) + 10

    String out_vcf = "/cromwell_root/" + basename(input_vcf, ".vcf.gz") + ".haplotypecaller.vcf.gz"

    command <<<
        set -euo pipefail

        cd /opt/ && \
            bash postprocess.sh ~{input_vcf} ~{out_vcf} ~{ref_fasta} && \
            cd -
    >>>

    output {
        File output_vcf = "~{out_vcf}"
        File output_vcf_index = "~{out_vcf}.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             7,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gatk:0.1.1"
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
