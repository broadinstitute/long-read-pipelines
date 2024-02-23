version 1.0

import "../../structs/Structs.wdl"
import "../../structs/ReferenceMetadata.wdl"
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
        haploid_contigs: "Optimization since DV 1.6 to improve calling on haploid contigs (e.g. human allosomes); see DV github page for more info."
        par_regions_bed: "The pseudoautosomal (PAR) regions of the human allosomes in BED format."
    }

    input {
        Array[Pair[String, Pair[File, File]]] how_to_shard_wg_for_calling
        String prefix
        String model_for_dv_andor_pepper
        String? haploid_contigs
        File? par_regions_bed

        # reference info
        HumanReferenceBundle ref_bundle

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

        Array[File] nongpu_resource_usage_logs = select_all(nongpu_resource_usage_log)
        Array[File] nongpu_resource_usage_visual = select_all(VisualizeDVRegularResoureUsage.plot_pdf)

        Array[File] native_visual_report_htmls = select_all(dv_self_visual_html)
    }

    ####################################################################################################################################

    #############
    # per shard DV
    #############
    scatter (triplet in how_to_shard_wg_for_calling) {
        String shard_id = triplet.left
        File shard_bam = triplet.right.left
        File shard_bai = triplet.right.right

        if (shard_id != "alts") {
        if (!use_gpu) {
            Boolean is_t2t_chrX = length(how_to_shard_wg_for_calling) < 20 && shard_id == "18_X"
            Boolean is_t2t_shard8 = length(how_to_shard_wg_for_calling) < 20 && shard_id == "9_15"
            Boolean is_38_chr1  = length(how_to_shard_wg_for_calling) > 20 && shard_id == "1-p"
            Boolean is_38_shard3  = length(how_to_shard_wg_for_calling) > 20 && shard_id == "11-q_17-q"

            Boolean pay_fee_and_go = is_t2t_chrX || is_38_chr1 || is_t2t_shard8 || is_38_shard3
            Int use_this_memory = if ( pay_fee_and_go ) then 48 else dv_memory
            # RuntimeAttr preemption_override = {"preemptible_tries": if ( pay_fee_and_go ) then 0 else 1}
            call DV as DeepV {
                input:
                    bam           = shard_bam,
                    bai           = shard_bai,
                    ref_fasta     = ref_bundle.fasta,
                    ref_fasta_fai = ref_bundle.fai,

                    haploid_contigs = haploid_contigs,
                    par_regions_bed = par_regions_bed,

                    model_type = model_for_dv_andor_pepper,

                    threads = select_first([dv_threads]),
                    memory  = use_this_memory, # select_first([dv_memory]),
                    zones = zones,
                    # runtime_attr_override = preemption_override
            }
        }

        if (use_gpu) {
            call DV_gpu as DeepV_G {
                input:
                    bam           = shard_bam,
                    bai           = shard_bai,
                    ref_fasta     = ref_bundle.fasta,
                    ref_fasta_fai = ref_bundle.fai,

                    haploid_contigs = haploid_contigs,
                    par_regions_bed = par_regions_bed,

                    model_type = model_for_dv_andor_pepper,

                    threads = select_first([dv_threads]),
                    memory  = select_first([dv_memory]),
                    zones = zones
            }
        }

        File dv_vcf = select_first([DeepV.VCF, DeepV_G.VCF])
        File dv_gvcf = select_first([DeepV.gVCF, DeepV_G.gVCF])
        File nongpu_resource_usage_log = select_first([DeepV.resouce_monitor_log, DeepV_G.resouce_monitor_log])
        File dv_self_visual_html = select_first([DeepV.visual_report_html, DeepV_G.visual_report_html])

        call VisualizeResourceUsage.SimpleRscript as VisualizeDVRegularResoureUsage {
            input:
            resource_log = nongpu_resource_usage_log,
            output_pdf_name = "~{prefix}.deepvariant.regular-resources-usage.~{triplet.left}.pdf",
            plot_title = "DeepVariant, on input ~{prefix}, at locus ~{triplet.left}"
        }
        }
    }

    #############
    # merge
    #############
    String dv_prefix = prefix + ".deepvariant"

    call VariantUtils.MergeAndSortVCFs as MergeDeepVariantGVCFs {
        input:
            vcfs     = select_all(dv_gvcf),
            prefix   = dv_prefix + ".g",
            ref_fasta_fai = ref_bundle.fai
    }

    call VariantUtils.MergeAndSortVCFs as MergeDeepVariantVCFs {
        input:
            vcfs     = select_all(dv_vcf),
            prefix   = dv_prefix,
            ref_fasta_fai = ref_bundle.fai
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

        String? haploid_contigs
        File? par_regions_bed

        String model_type

        Int threads
        Int memory
        String zones

        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, ".bam") + ".deepvariant"
    String output_root = "/cromwell_root/dv_output"

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
            ~{true='--haploid_contigs ' false='' defined(haploid_contigs)} ~{select_first([haploid_contigs, ""])} \
            ~{true='--par_regions_bed ' false='' defined(par_regions_bed)} ~{select_first([par_regions_bed, ""])} \
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
    Int bam_sz = ceil(size(bam, "GB"))
    Boolean is_big_bam = bam_sz > 100
    Int inflation_factor = if (is_big_bam) then 10 else 5
    Int minimal_disk = 20
	Int disk_size = if inflation_factor * bam_sz > minimal_disk then inflation_factor * bam_sz else minimal_disk

    RuntimeAttr default_attr = object {
        cpu_cores:          threads,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-deepvariant:1.6.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones: zones
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

        String? haploid_contigs
        File? par_regions_bed

        String model_type

        Int threads
        Int memory
        String zones

        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, ".bam") + ".deepvariant"
    String output_root = "/cromwell_root/dv_output"

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
            ~{true='--haploid_contigs ' false='' defined(haploid_contigs)} ~{select_first([haploid_contigs, ""])} \
            ~{true='--par_regions_bed ' false='' defined(par_regions_bed)} ~{select_first([par_regions_bed, ""])} \
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

    Int bam_sz = ceil(size(bam, "GB"))
    Boolean is_big_bam = bam_sz > 100
    Int inflation_factor = if (is_big_bam) then 10 else 5
    Int minimal_disk = 100
	Int disk_size = if inflation_factor * bam_sz > minimal_disk then inflation_factor * bam_sz else minimal_disk

    Int max_cpu = 12
    Int use_this_cpu = if threads > max_cpu then max_cpu else threads
    Int max_memory = 64
    Int use_this_memory = if memory > max_memory then max_memory else memory

    RuntimeAttr default_attr = object {
        cpu_cores:          use_this_cpu,
        mem_gb:             use_this_memory,
        disk_gb:            disk_size,
        boot_disk_gb:       30,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-deepvariant:1.6.0-gpu"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])

        zones: zones
        gpuType: "nvidia-tesla-v100"
        gpuCount: 1
    }
}
