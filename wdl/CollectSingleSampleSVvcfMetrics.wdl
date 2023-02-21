version 1.0

import "tasks/Structs.wdl"
import "tasks/VariantUtils.wdl"

workflow CollectSingleSampleSVvcfMetrics {
    meta {
        description: "Collect single sample SV VCF metrics. Only a few selected callers are supported."
    }
    input {
        File vcf
        File tbi
        String caller
        Array[String] variant_types
        File? baseline_vcf

        File ref_map_file
        File contig_list
        String sv_pipeline_base_docker
    }

    parameter_meta {
        vcf: "The single sample SV VCF, gzipped."
        tbi: "The accompanying index file."
        caller: "The caller used to make the VCF. Currently supports [pbsv, sniffles2, pav]"
        variant_types: "A list of SVTYPE's expected to appear in the VCF."
        baseline_vcf: "TO BE FILLED"
        ref_fai: ""
        contig_list: "TO BE FILLED"
        sv_pipeline_base_docker: "Docker image for 'sv-pipeline-base' in the GATK-SV pipeline"

        sv_vcf_metrics: "The metrics in a Map format. We recommend changing the name of this out on Terra to match the caller."
    }

    Map[String, String] ref_map = read_map(ref_map_file)
    call VariantUtils.GetVCFSampleName {
        input: fingerprint_vcf = vcf
    }
    call Standardize {
        input:
            vcf = vcf, tbi = tbi, ref_fai = ref_map['fai'], caller = caller, prefix = GetVCFSampleName.sample_name
    }

    String dummy_prefix = "placeholder"
    call VCFMetrics {
        input:
            vcf = Standardize.standardized_vcf, baseline_vcf = baseline_vcf,
            samples = [GetVCFSampleName.sample_name], variant_types = variant_types,
            prefix = dummy_prefix, contig_list = contig_list, sv_pipeline_base_docker = sv_pipeline_base_docker
    }
    call CleanMap {
        input: metrics_tsv = VCFMetrics.out, prefix_used = dummy_prefix
    }

    output {
        Map[String, Int] sv_vcf_metrics = CleanMap.metrics
    }
}

task Standardize {
    meta {
        description: "Standarize/format a selected number of caller's SV VCF, so that metrics can be collected."
    }
    input {
        File vcf
        File tbi
        File ref_fai

        String caller
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 8*ceil(size([vcf, tbi, ref_fai], "GB")) + 1

    command <<<
        set -euxo pipefail

        svtk standardize \
            --include-reference-sites \
            --contigs ~{ref_fai} \
            --prefix ~{prefix} ~{vcf} - ~{caller} | \
            bcftools sort /dev/stdin -o ~{prefix}.truvari.std.vcf.gz -O z

        tabix ~{prefix}.truvari.std.vcf.gz
    >>>

    output {
        File standardized_vcf = "~{prefix}.truvari.std.vcf.gz"
        File standardized_tbi = "~{prefix}.truvari.std.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             24,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-truvari:shuang-sp"
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

task VCFMetrics {
    meta {
        description: "Collect metrics from a pre-formatted SV VCF. Source: https://github.com/broadinstitute/gatk-sv/blob/43c4e3f2e2b3e85e065269b9932bc144c67cded2/wdl/TestUtils.wdl#L87"
    }
    input {
        File vcf
        File? baseline_vcf
        Array[String] samples
        Array[String] variant_types
        String prefix
        File contig_list
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    File samples_list = write_lines(samples)

    RuntimeAttr runtime_default = object {
        mem_gb: 3.75,
        disk_gb: 10,
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    output {
        File out = "~{prefix}.vcf.tsv"
    }
    command <<<

        svtest vcf \
         ~{vcf} \
         ~{contig_list} \
         ~{samples_list} \
         ~{sep=',' variant_types} \
         ~{prefix} \
         ~{if defined(baseline_vcf) then "--baseline-vcf " + baseline_vcf else ""} \
         > ~{prefix}.vcf.tsv

    >>>
    runtime {
        memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_base_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
}

task CleanMap {
    meta {
        description: "A specific task to cleanup the TSV output by the GATK_SV task"
    }
    input {
        File metrics_tsv
        String prefix_used
    }

    command <<<
        set -eux

        sed "s#~{prefix_used}_vcf_##" ~{metrics_tsv} > cleaned.2col.tsv
        cat cleaned.2col.tsv

        awk -F '\t' '{print $1}' cleaned.2col.tsv > names.txt
        awk -F '\t' '{print $2}' cleaned.2col.tsv > values.txt
    >>>

    output {
        Array[String] attrs  = read_lines("names.txt")
        Array[Int]    values = read_lines("values.txt")
        Map[String, Int] metrics = read_map("cleaned.2col.tsv")
    }
    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
