version 1.0

import "../../../structs/Structs.wdl"

import "../../../tasks/Utility/Finalize.wdl" as FF

workflow HiFiCNV {
    meta {
        description:
        "Runs the PacBio HiFiCNV tool on a single (human) HiFi bam."
    }
    parameter_meta {
        exclude_bed:      "BED holding regions that are known to cause artifacts during HiFiCNV data processing (e.g. centromeres)."
        sex_specific_cn:  "Sex-specific files annotating the PAR regions on and expected copy numbers of sex chromosomes."
        ref_map_file:     "table indicating reference sequence and auxiliary file locations"
        gcs_out_root_dir: "GCS bucket to store the variants, and metrics files"
    }
    input {
        File bam
        File bai

        File ref_map_file
        File exclude_bed
        File sex_specific_cn

        String gcs_out_root_dir
    }
    output {
        String vcf      = FinalizeVCF.gcs_path
        String vcf_tbi  = FinalizeVcfIndex.gcs_path
        String bedgraph = FinalizeBedGraph.gcs_path
        String depth_bw = FinalizeBigWig.gcs_path
        String log      = FinalizeLog.gcs_path
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    # Use the BAM filename as the hificnv output prefix. The sample name is
    # inferred inside the PacBioHiFiCNV task from the *localized* BAM (see there),
    # replacing the separate InferSampleName task that streamed the BAM header
    # directly from GCS and failed on requester-pays buckets.
    String prefix = basename(bam, ".bam")

    call PacBioHiFiCNV { input:
        bam = bam, bai = bai,
        output_prefix = prefix,
        ref_fasta = ref_map['fasta'],
        ref_fasta_fai = ref_map['fai'],
        exclude_bed = exclude_bed,
        sex_specific_cn = sex_specific_cn
    }

    String workflow_name = 'HiFiCNV'
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_name}/~{prefix}"
    call FF.FinalizeToFile as FinalizeLog      { input: outdir = outdir, file = PacBioHiFiCNV.log }
    call FF.FinalizeToFile as FinalizeVCF      { input: outdir = outdir, file = PacBioHiFiCNV.vcf }
    call FF.FinalizeToFile as FinalizeVcfIndex { input: outdir = outdir, file = PacBioHiFiCNV.vcf_tbi }
    call FF.FinalizeToFile as FinalizeBedGraph { input: outdir = outdir, file = PacBioHiFiCNV.bedgraph }
    call FF.FinalizeToFile as FinalizeBigWig   { input: outdir = outdir, file = PacBioHiFiCNV.depth_bw }
}

task PacBioHiFiCNV {
    input {
        File bam
        File bai
        String output_prefix
        File ref_fasta
        File ref_fasta_fai

        File exclude_bed
        File sex_specific_cn

        RuntimeAttr? runtime_attr_override
    }

    output {
        # The command renames hificnv's per-sample outputs
        # (<output_prefix>.<sample>.<ext>) to deterministic <output_prefix>.<ext>.
        # Output paths MUST depend only on output_prefix (a task input): the GCP
        # Batch backend builds the delocalization file list BEFORE the command
        # runs, so a name derived from the runtime-inferred sample (e.g. via
        # read_string) resolves to empty and delocalization of the "required"
        # outputs fails the job (with a misleading exit-code-0 message).
        File vcf      = "~{output_prefix}.vcf.gz"
        File vcf_tbi  = "~{output_prefix}.vcf.gz.tbi"
        File bedgraph = "~{output_prefix}.copynum.bedgraph"
        File depth_bw = "~{output_prefix}.depth.bw"
        File log      = "~{output_prefix}.log"
    }

    command <<<
    set -eux

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        # re-generate the bai to ensure it's not corrupted, as corrupted bai can cause hificnv to fail without a clear error message
        rm ~{bai} \
        && \
        samtools index ~{bam}

        # Infer the sample name from the (already-localized) BAM @RG SM tag -- the
        # same value hificnv uses to name its outputs. Doing this on the local
        # file avoids the direct-from-GCS header streaming that made the separate
        # InferSampleName task fail on requester-pays buckets. Same checks as that
        # task: fail on missing / multiple / unnamedsample.
        samtools view -H ~{bam} | grep '^@RG' | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g' | sort | uniq > sample_name.txt
        if [[ $(wc -l < sample_name.txt) -lt 1 ]]; then echo "No @RG SM sample name found!" && exit 1; fi
        if [[ $(wc -l < sample_name.txt) -gt 1 ]]; then echo "Multiple sample names found!" && exit 1; fi
        if grep -iq "unnamedsample" sample_name.txt; then echo "Sample name found to be unnamedsample!" && exit 1; fi
        sample_name=$(cat sample_name.txt)

        hificnv \
            --bam ~{bam} \
            --ref ~{ref_fasta} \
            --exclude ~{exclude_bed} \
            --expected-cn ~{sex_specific_cn} \
            --threads "${num_core}" \
            --output-prefix ~{output_prefix}

        # hificnv names the per-sample outputs <output_prefix>.<sample>.<ext>.
        # Rename them to deterministic <output_prefix>.<ext> so the declared
        # outputs depend only on output_prefix, not on the runtime-inferred
        # sample name (see the output block for why that matters on GCP Batch).
        # (The log is already named <output_prefix>.log, so it needs no rename.)
        mv ~{output_prefix}.${sample_name}.vcf.gz           ~{output_prefix}.vcf.gz
        mv ~{output_prefix}.${sample_name}.copynum.bedgraph ~{output_prefix}.copynum.bedgraph
        mv ~{output_prefix}.${sample_name}.depth.bw         ~{output_prefix}.depth.bw

        # hificnv does not emit a VCF index; create a tabix index so the .vcf.gz
        # can be random-accessed downstream.
        tabix -p vcf ~{output_prefix}.vcf.gz

        tree
    >>>

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
        docker:             "us.gcr.io/broad-dsp-lrma/hificnv:1.0.1"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])

        noAddress: true
    }
}
