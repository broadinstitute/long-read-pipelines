version 1.0

##################################################
# This is essentially copying "dsde_pipelines_tasks/VariantCalling.wdl"
# with some customization to fit the process described in 
# "https://github.com/PacificBiosciences/hg002-ccs/"
##################################################

import "dsde_pipelines_tasks/GermlineVariantDiscovery.wdl" as Calling
import "dsde_pipelines_tasks/Qc.wdl" as QC
import "dsde_pipelines_tasks/Utilities.wdl" as DSDEPipelinesUtils

import "Utils.wdl" as Utils

workflow GATKBestPraciceForLR {

    input {
      File calling_interval_list
      Int haplotype_scatter_count
      Int break_bands_at_multiples_of

      Float? contamination

      File input_bam
      File ref_fasta
      File ref_fasta_index
      File ref_dict

      Boolean sample_is_female

      String gatk4_docker_tag

      String base_file_name
      String final_vcf_base_name

      Boolean run_qc_on_variants
      Boolean make_gvcf = true
      Boolean make_bamout = false

      Int agg_preemptible_tries = 1 # can not be optional as the tasks and sub-workflows copied over has this non-optional
    }

    parameter_meta {
      make_bamout: "For CNNScoreVariants to run with a 2D model, a bamout must be created by HaplotypeCaller. The bamout is a bam containing information on how HaplotypeCaller remapped reads while it was calling variants. See https://gatkforums.broadinstitute.org/gatk/discussion/5484/howto-generate-a-bamout-file-showing-how-haplotypecaller-has-remapped-sequence-reads for more details."
    }

    ###########################################################################
    # Break the calling interval_list into sub-intervals
    # Perform variant calling on the sub-intervals, and then gather the results
    call DSDEPipelinesUtils.ScatterIntervalList as ScatterIntervalList {
      input:
        interval_list = calling_interval_list,
        scatter_count = haplotype_scatter_count,
        break_bands_at_multiples_of = break_bands_at_multiples_of
    }

    ###########################################################################
    # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
    # If we take the number we are scattering by and reduce by 20 we will have enough disk space
    # to account for the fact that the data is quite uneven across the shards.
    Int potential_hc_divisor = ScatterIntervalList.interval_count - 20
    Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

    # Call variants in parallel over WGS calling intervals
    scatter (scattered_interval_list in ScatterIntervalList.out) {

        # Generate GVCF by interval
        call Calling.HaplotypeCaller_GATK4_VCF as HaplotypeCallerGATK4 {
          input:
            contamination = contamination,
            input_bam = input_bam,
            interval_list = scattered_interval_list,
            vcf_basename = base_file_name,
            ref_dict = ref_dict,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            hc_scatter = hc_divisor,
            make_gvcf = make_gvcf,
            make_bamout = make_bamout,
            preemptible_tries = agg_preemptible_tries,

            sample_is_female = sample_is_female,

            gatk4_docker_tag = gatk4_docker_tag
         }

        # If bamout files were created, we need to sort and gather them into one bamout
        if (make_bamout) {
          call Utils.SortSam as SortBamout {
            input:
              input_bam = HaplotypeCallerGATK4.bamout,
              output_bam_basename = final_vcf_base_name,
              compression_level = 2
          }
        }

      File vcfs_to_merge = HaplotypeCallerGATK4.output_vcf
      File vcf_indices_to_merge = HaplotypeCallerGATK4.output_vcf_index
    }

    ###########################################################################
    # Combine by-interval (g)VCFs into a single sample (g)VCF file
    String merge_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
    call Calling.MergeVCFs as MergeVCFs {
      input:
        input_vcfs = vcfs_to_merge,
        input_vcfs_indexes = vcf_indices_to_merge,
        output_vcf_name = final_vcf_base_name + merge_suffix,
        preemptible_tries = agg_preemptible_tries
    }

    ###########################################################################
    if (make_bamout) {
      call MergeBamouts {
        input:
          bams = select_all(SortBamout.output_bam),
          output_base_name = final_vcf_base_name
      }
    }

    ###########################################################################
    if (run_qc_on_variants) {
        # Validate the (g)VCF output of HaplotypeCaller
        call QC.ValidateVCF as ValidateVCF {
          input:
            input_vcf = MergeVCFs.output_vcf,
            input_vcf_index = MergeVCFs.output_vcf_index,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            calling_interval_list = calling_interval_list,
            is_gvcf = make_gvcf,
            preemptible_tries = agg_preemptible_tries
        }

        # QC the (g)VCF
        call QC.CollectVariantCallingMetrics as CollectVariantCallingMetrics {
          input:
            input_vcf = MergeVCFs.output_vcf,
            input_vcf_index = MergeVCFs.output_vcf_index,
            metrics_basename = final_vcf_base_name,
            ref_dict = ref_dict,
            is_gvcf = make_gvcf,
            preemptible_tries = agg_preemptible_tries
        }
    }

    output {
      File output_vcf = MergeVCFs.output_vcf
      File output_vcf_index = MergeVCFs.output_vcf_index
      File? vcf_summary_metrics = CollectVariantCallingMetrics.summary_metrics
      File? vcf_detail_metrics = CollectVariantCallingMetrics.detail_metrics
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
        disk_gb:            "local-disk ~{disk_size} HDD",
        boot_disk_gb:       10,
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
