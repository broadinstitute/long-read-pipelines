version 1.0

import "../../structs/Structs.wdl"
import "../Utility/VariantUtils.wdl"
import "../Visualization/VisualizeResourceUsage.wdl"

workflow Run {
    meta {
        desciption:
        "For running DeepVariant on CCS/ONT WGS data."
    }
    parameter_meta {
        how_to_shard_wg_for_calling: "An array of the BAM's shard; each element is assumed to be a tuple of (ID for the shard, (BAM of the shard, BAI of the shard))"
        model_for_dv_andor_pepper: "Model string to be used on DV or the PEPPER-Margin-DeepVariant toolchain. Please refer to their github pages for accepted values."
        prefix: "Prefix for output files"
    }

    input {
        Array[Pair[String, Pair[File, File]]] how_to_shard_wg_for_calling
        String prefix
        String model_for_dv_andor_pepper

        # reference info
        Map[String, String] ref_map

        # optimizations
        Int dv_threads
        Int dv_memory
        Boolean use_gpu
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }

    output {
        File g_vcf = MergeDeepVariantGVCFs.vcf
        File g_tbi = MergeDeepVariantGVCFs.tbi
        File vcf = MergeDeepVariantVCFs.vcf
        File tbi = MergeDeepVariantVCFs.tbi

        Array[File] nongpu_resource_usage_logs = nongpu_resource_usage_log
        Array[File] nongpu_resource_usage_visual = VisualizeDVRegularResoureUsage.plot_pdf
    }

    ####################################################################################################################################

    #############
    # per shard DV
    #############
    scatter (triplet in how_to_shard_wg_for_calling) {
        String shard_id = triplet.left
        File shard_bam = triplet.right.left
        File shard_bai = triplet.right.right

        if (!use_gpu) {
            call DV as DeepV {
                input:
                    bam           = shard_bam,
                    bai           = shard_bai,
                    ref_fasta     = ref_map['fasta'],
                    ref_fasta_fai = ref_map['fai'],

                    model_type = model_for_dv_andor_pepper,

                    threads = select_first([dv_threads]),
                    memory  = select_first([dv_memory]),
                    zones = zones
            }
        }

        if (use_gpu) {
            call DV_gpu as DeepV_G {
                input:
                    bam           = shard_bam,
                    bai           = shard_bai,
                    ref_fasta     = ref_map['fasta'],
                    ref_fasta_fai = ref_map['fai'],

                    model_type = model_for_dv_andor_pepper,

                    threads = select_first([dv_threads]),
                    memory  = select_first([dv_memory]),
                    zones = zones
            }
        }

        File dv_vcf = select_first([DeepV.VCF, DeepV_G.VCF])
        File dv_gvcf = select_first([DeepV.gVCF, DeepV_G.gVCF])
        File nongpu_resource_usage_log = select_first([DeepV.resouce_monitor_log, DeepV_G.resouce_monitor_log])

        call VisualizeResourceUsage.SimpleRscript as VisualizeDVRegularResoureUsage {
            input:
            resource_log = nongpu_resource_usage_log,
            output_pdf_name = "~{prefix}.deepvariant.regular-resources-usage.~{triplet.left}.pdf",
            plot_title = "DeepVariant, on input ~{prefix}, at locus ~{triplet.left}"
        }
    }

    #############
    # merge
    #############
    String dv_prefix = prefix + ".deepvariant"

    call VariantUtils.MergeAndSortVCFs as MergeDeepVariantGVCFs {
        input:
            vcfs     = dv_gvcf,
            prefix   = dv_prefix + ".g",
            ref_fasta_fai = ref_map['fai']
    }

    call VariantUtils.MergeAndSortVCFs as MergeDeepVariantVCFs {
        input:
            vcfs     = dv_vcf,
            prefix   = dv_prefix,
            ref_fasta_fai = ref_map['fai']
    }
}

task DV {

    parameter_meta {
        model_type: "which DV pre-trained model to use. Must be one of [PACBIO, ONT_R104] (or anything later supported after DV's 1.5.0 release)."
    }
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        String model_type

        Int threads
        Int memory
        String zones

        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, ".bam") + ".deepvariant"
    String output_root = "/cromwell_root/dv_output"

    Int bam_sz = ceil(size(bam, "GB"))
    Boolean is_big_bam = bam_sz > 100
    Int inflation_factor = if (is_big_bam) then 10 else 5
    Int minimal_disk = 50
	Int disk_size = if inflation_factor * bam_sz > minimal_disk then inflation_factor * bam_sz else minimal_disk

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        mkdir -p "~{output_root}"

        export MONITOR_MOUNT_POINT="/cromwell_root/"
        bash /opt/vm_local_monitoring_script.sh &> resources.log &
        job_id=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')

        /opt/deepvariant/bin/run_deepvariant \
            --model_type=~{model_type} \
            --ref=~{ref_fasta} \
            --reads=~{bam} \
            --output_vcf="~{output_root}/~{prefix}.vcf.gz" \
            --output_gvcf="~{output_root}/~{prefix}.g.vcf.gz" \
            --num_shards="${num_core}" || cat resources.log
        if ps -p "${job_id}" > /dev/null; then kill "${job_id}"; fi

        find "~{output_root}/" -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g' \
            > "~{output_root}/dir_structure.txt"
    >>>

    output {

        File resouce_monitor_log = "resources.log"
        File? gpu_monitor_log = "gpu.usages.log"
        File output_dir_structure = "~{output_root}/dir_structure.txt"

        File VCF        = "~{output_root}/~{prefix}.vcf.gz"
        File VCF_tbi    = "~{output_root}/~{prefix}.vcf.gz.tbi"

        File gVCF       = "~{output_root}/~{prefix}.g.vcf.gz"
        File gVCF_tbi   = "~{output_root}/~{prefix}.g.vcf.gz.tbi"

        File visual_report_html = "~{output_root}/~{prefix}.visual_report.html"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          threads,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-deepvariant:1.5.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task DV_gpu {
    parameter_meta {
        model_type: "which DV pre-trained model to use. Must be one of [PACBIO, ONT_R104] (or anything later supported after DV's 1.5.0 release)."
    }

    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        String model_type

        Int threads
        Int memory
        String zones

        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, ".bam") + ".deepvariant"
    String output_root = "/cromwell_root/dv_output"

    Int bam_sz = ceil(size(bam, "GB"))
    Boolean is_big_bam = bam_sz > 100
    Int inflation_factor = if (is_big_bam) then 10 else 5
    Int minimal_disk = 100
	Int disk_size = if inflation_factor * bam_sz > minimal_disk then inflation_factor * bam_sz else minimal_disk

    Int max_cpu = 12
    Int use_this_cpu = if threads > max_cpu then max_cpu else threads
    Int max_memory = 64
    Int use_this_memory = if memory > max_memory then max_memory else memory

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        mkdir -p "~{output_root}"

        export MONITOR_MOUNT_POINT="/cromwell_root/"
        bash vm_local_monitoring_script.sh &> resources.log &
        job_id=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')
        gpustat -a -i 1 &> gpu.usages.log &
        gpu_tracking_job_id=$(ps -aux | grep -F 'gpustat' | head -1 | awk '{print $2}')

        /opt/deepvariant/bin/run_deepvariant \
            --model_type=~{model_type} \
            --ref=~{ref_fasta} \
            --reads=~{bam} \
            --output_vcf="~{output_root}/~{prefix}.vcf.gz" \
            --output_gvcf="~{output_root}/~{prefix}.g.vcf.gz" \
            --num_shards="${num_core}" || cat resources.log
        if ps -p "${job_id}" > /dev/null; then kill "${job_id}"; fi
        if ps -p "${gpu_tracking_job_id}" > /dev/null; then kill "${gpu_tracking_job_id}"; fi

        find "~{output_root}/" -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g' \
            > "~{output_root}/dir_structure.txt"
    >>>

    output {

        File resouce_monitor_log = "resources.log"
        File? gpu_monitor_log = "gpu.usages.log"
        File output_dir_structure = "~{output_root}/dir_structure.txt"

        File VCF        = "~{output_root}/~{prefix}.vcf.gz"
        File VCF_tbi    = "~{output_root}/~{prefix}.vcf.gz.tbi"

        File gVCF       = "~{output_root}/~{prefix}.g.vcf.gz"
        File gVCF_tbi   = "~{output_root}/~{prefix}.g.vcf.gz.tbi"

        File visual_report_html = "~{output_root}/~{prefix}.visual_report.html"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          use_this_cpu,
        mem_gb:             use_this_memory,
        disk_gb:            disk_size,
        boot_disk_gb:       30,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-deepvariant:1.5.0-gpu"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])

        gpuType: "nvidia-tesla-v100"
        gpuCount: 1
    }
}
