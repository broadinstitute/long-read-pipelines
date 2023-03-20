version 1.0

import "../../../structs/Structs.wdl"

import "../../../tasks/Utility/Utils.wdl"

import "../../../tasks/Utility/Finalize.wdl" as FF

workflow HiFiCNV {
    meta {
        desciption:
        "Runs the PacBio HiFiCNV tool on a single (human) HiFi bam."
    }
    input {
        File bam
        File bai
        File ref_map_file
        File exclude_bed
        File sex_specific_cn
        String gcs_out_root_dir
    }
    parameter_meta {
        exclude_bed:      "BED holding regions that are known to cause artifacts during HiFiCNV data processing (e.g. centromeres)."
        sex_specific_cn:  "Sex-specific files annotating the PAR regions on and expected copy numbers of sex chromosomes."
        ref_map_file:     "table indicating reference sequence and auxillary file locations"
        gcs_out_root_dir: "GCS bucket to store the reads, variants, and metrics files"
    }

    String workflow_name = 'HiFiCNV'
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_name}"

    Map[String, String] ref_map = read_map(ref_map_file)

    call Utils.InferSampleName { input: bam = bam, bai = bai}
    call PacBioHiFiCNV { input:
        bam = bam, bai = bai,
        sample_name = InferSampleName.sample_name,
        output_prefix = InferSampleName.sample_name,
        ref_fasta = ref_map['fasta'],
        ref_fasta_fai = ref_map['fai'],
        exclude_bed = exclude_bed,
        sex_specific_cn = sex_specific_cn
    }

    call FF.FinalizeToFile as FinalizeLog { input: outdir = outdir + '/~{InferSampleName.sample_name}', file = PacBioHiFiCNV.log }
    call FF.FinalizeToFile as FinalizeVCF { input: outdir = outdir + '/~{InferSampleName.sample_name}', file = PacBioHiFiCNV.vcf }
    call FF.FinalizeToFile as FinalizeBedGraph { input: outdir = outdir + '/~{InferSampleName.sample_name}', file = PacBioHiFiCNV.bedgraph }
    call FF.FinalizeToFile as FinalizeBigWig { input: outdir = outdir + '/~{InferSampleName.sample_name}', file = PacBioHiFiCNV.depth_bw }

    output {
        Map[String, String] hificnv_outs = {'vcf': FinalizeVCF.gcs_path,
                                            'bedgraph': FinalizeLog.gcs_path,
                                            'log': FinalizeBedGraph.gcs_path,
                                            'depth_bw': FinalizeBigWig.gcs_path
                                            }
    }
}

task PacBioHiFiCNV {
    input {
        File bam
        File bai
        String sample_name
        String output_prefix
        File ref_fasta
        File ref_fasta_fai

        File exclude_bed
        File sex_specific_cn
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eux

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        hificnv \
            --bam ~{bam} \
            --ref ~{ref_fasta} \
            --exclude ~{exclude_bed} \
            --expected-cn ~{sex_specific_cn} \
            --threads "${num_core}" \
            --output-prefix ~{output_prefix}

        tree
    >>>
    output {
        File vcf = "~{output_prefix}.${sample_name}.vcf.gz"
        File bedgraph = "~{output_prefix}.${sample_name}.copynum.bedgraph"
        File log = "~{output_prefix}.log"
        File depth_bw = "~{output_prefix}.${sample_name}.depth.bw"
    }

    #########################
    Int min_disk = 40
    Float disk_multiplier = 1
    Int disk_size = ceil(disk_multiplier * size(bam, "GiB")) + 20
    Int use_this_disk_sz = if (min_disk>disk_size) then min_disk else disk_size

    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             6,
        disk_gb:            use_this_disk_sz,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/hificnv:1.0.0"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
