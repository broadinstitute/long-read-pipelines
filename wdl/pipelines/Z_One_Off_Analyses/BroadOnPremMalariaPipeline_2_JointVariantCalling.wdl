version 1.0

import "../../tasks/Utility/Utils.wdl" as Utils
import "../../structs/Structs.wdl"

import "../../tasks/Z_One_Off_Analyses/BroadOnPremMalariaPipelineTasks.wdl" as BroadOnPremMalariaPipelineTasks

workflow BroadOnPremMalariaPipeline_2_JointVariantCalling {
    meta {
        desciption: "Recreation of the second step of the Broad OnPrem Malaria Pipeline.  Joint calls and recalibrates/filters variants."
    }
    input {
        String sample_name

        Array[File] vcf_files
        Array[File] vcf_index_files

        File ref_map_file

        File resource_vcf_7g8_gb4
        File resource_vcf_hb3_dd2
        File resource_vcf_3d7_hb3
    }

    # Get ref info:
    Map[String, String] ref_map = read_map(ref_map_file)

    call GenotypeGVCFs as t_001_GenotypeGVCFs {
        input:
            prefix = sample_name,
            input_vcfs = vcf_files,
            input_vcf_indices = vcf_index_files,
            reference_fasta = ref_map["fasta"],
            reference_fai = ref_map["fai"],
            reference_dict = ref_map["dict"]
    }

    call BroadOnPremMalariaPipelineTasks.VariantRecalibrator as t_002_VariantRecalibrator {
        input:
            prefix = sample_name,
            input_vcf = t_001_GenotypeGVCFs.vcf,
            reference_fasta = ref_map["fasta"],
            reference_fai = ref_map["fai"],
            reference_dict = ref_map["dict"],
            resource_vcf_7g8_gb4 = resource_vcf_7g8_gb4,
            resource_vcf_hb3_dd2 = resource_vcf_hb3_dd2,
            resource_vcf_3d7_hb3 = resource_vcf_3d7_hb3
    }

    # 9 - sort, compress, and index final outputs:
    call BroadOnPremMalariaPipelineTasks.SortCompressIndexVcf as t_003_SortCompressIndexVcf {
        input:
            input_vcf = t_002_VariantRecalibrator.vcf,
    }

    output {
        File joint_vcf = t_003_SortCompressIndexVcf.vcf
        File joint_vcf_index = t_003_SortCompressIndexVcf.vcf_index
    }
}

task GenotypeGVCFs {
    input {
        String prefix

        Array[File] input_vcfs
        Array[File] input_vcf_indices

        File reference_fasta
        File reference_fai
        File reference_dict

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 10*ceil(
        size(input_vcfs, "GB")
        + size(input_vcf_indices, "GB")
        + size(reference_fasta, "GB")
        + size(reference_fai, "GB")
        + size(reference_dict, "GB")
    )

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

        java -Xmx${java_memory_size_mb}M -jar /usr/GenomeAnalysisTK.jar \
            -T GenotypeGVCFs \
                -R ~{reference_fasta} \
                -o ~{prefix}_CombinedGVCFs.vcf \
                -V ~{sep=" -V " input_vcfs}
    >>>

    output {
        File vcf = "~{prefix}_CombinedGVCFs.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size, # This uses the variable calculated above
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
