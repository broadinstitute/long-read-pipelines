version 1.0

import "tasks/utils/AlignAndCheckFingerprintCCS.wdl" as major

import "tasks/Utils.wdl"
import "tasks/PBUtils.wdl"
import "tasks/utils/BAMutils.wdl"
import "tasks/utils/GeneralUtils.wdl" as GU
import "tasks/AlignedMetrics.wdl"

import "tasks/Finalize.wdl" as FF

workflow ProcessOnInstrumentDemuxedChunk {

    meta {
        desciption: "!!! WARN: THIS IS PROJECT-CENTER SPECIFIC !!! Given an on-instrument demultiplexed hifi_reads.bam, perform alignment and QC check."
    }

    input {
        File uBAM

        String smrtcell_id
        String barcode_name

        String biosample_id
        String sample_id
        String downstream_sample_id

        String fingerprint_store
        String sample_id_at_store
        Boolean turn_off_fingperprint_check = false

        File ref_map_file

        String gcs_out_root_dir
    }

    ###################################################################################
    # prep work

    # where to store final results
    String workflow_name = "ProcessOnInstrumentDemuxedChunk"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name
    String outdir_aln     = outdir + '/alignments'
    String outdir_metrics = outdir + '/metrics'

    ###################################################################################
    # generate PBI
    # update header so as to be conformant with other input formats
    call Utils.ComputeAllowedLocalSSD {input: intended_gb = ceil(2 * size(uBAM, "GiB"))}
    call ReheaderInputUnalignedBAM as RH {input: uBAM = uBAM, sample_name = downstream_sample_id, ssd_cnt = ComputeAllowedLocalSSD.numb_of_local_ssd}
    call PBUtils.PBIndex as Index {input: bam = RH.reheadered_bam}

    call BAMutils.GetReadGroupInfo as RG {input: uBAM = RH.reheadered_bam, keys = ['SM', 'LB', 'PU']}

    # major work
    call major.AlignAndCheckFingerprintCCS {
        input:
            uBAM = uBAM,
            uPBI = Index.pbi,
            bam_sample_name = RG.read_group_info['SM'],
            library = RG.read_group_info['LB'],

            turn_off_fingperprint_check = turn_off_fingperprint_check,
            fp_store = fingerprint_store,
            sample_id_at_store = sample_id_at_store,
            ref_map_file = ref_map_file
    }

    call AlignedMetrics.MosDepthWGS { input: bam = AlignAndCheckFingerprintCCS.aligned_bam, bai = AlignAndCheckFingerprintCCS.aligned_bai}
    ###################################################################################
    # finalize
    String movie_name = RG.read_group_info['PU']
    String bc = barcode_name
    String bc_specific_aln_out    = outdir + '/alignments/' + smrtcell_id + '/' + bc
    String bc_specific_metric_out = outdir + "/metrics/"    + smrtcell_id + '/' + bc

    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = bc_specific_aln_out, file = AlignAndCheckFingerprintCCS.aligned_bam, name = movie_name + '.' + bc + '.bam' }
    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = bc_specific_aln_out, file = AlignAndCheckFingerprintCCS.aligned_bai, name = movie_name + '.' + bc + '.bai' }
    call FF.FinalizeToFile as FinalizeAlignedPbi { input: outdir = bc_specific_aln_out, file = AlignAndCheckFingerprintCCS.aligned_pbi, name = movie_name + '.' + bc + '.pbi' }

    call FF.FinalizeToFile as FinalizeAlnMetrics { input: outdir = bc_specific_metric_out, file = AlignAndCheckFingerprintCCS.alignment_metrics_tar_gz }

    if (! turn_off_fingperprint_check) {
        call FF.FinalizeToFile as FinalizeFPDetails  { input: outdir = bc_specific_metric_out, file = select_first([AlignAndCheckFingerprintCCS.fingerprint_detail_tar_gz]) }
    }

    ###################################################################################

    call GU.GetTodayDate as today {}

    output {
        String? fingerprint_check_results = AlignAndCheckFingerprintCCS.fp_status
        Float? fingerprint_check_LOD = AlignAndCheckFingerprintCCS.fp_lod_expected_sample
        File? fingerprint_check_tar_gz = FinalizeFPDetails.gcs_path

        File aligned_bam = FinalizeAlignedBam.gcs_path
        File aligned_bai = FinalizeAlignedBai.gcs_path
        File aligned_pbi = FinalizeAlignedPbi.gcs_path
        Float wgs_cov = MosDepthWGS.wgs_cov
        File alignment_metrics_tar_gz = FinalizeAlnMetrics.gcs_path
        Map[String, Float] alignment_metrics = AlignAndCheckFingerprintCCS.alignment_metrics
        String movie = movie_name

        String last_postprocessing_date = today.yyyy_mm_dd
    }
}

task ReheaderInputUnalignedBAM {
    meta {
        desciption: "Reheader the input unaligned BAM's SM entry in the readgroup line, so that all pipelines follow the same convention."
    }
    input {
        File uBAM
        String sample_name

        Int ssd_cnt
    }

    String prefix = basename(uBAM, ".bam")
    Int ssd_sz = ssd_cnt * 375

    command <<<
        set -eux

        samtools view --no-PG -H ~{uBAM} > header.txt
        awk '$1 ~ /^@RG/' header.txt > rg_line.txt

        # fix SM:
        if ! grep -qF "SM:" rg_line.txt; then
            sed -i "s/$/SM:tbd/" rg_line.txt
        fi
        awk -v sm="~{sample_name}" -F '\t' 'BEGIN {OFS="\t"} { for (i=1; i<=NF; ++i) { if ($i ~ "SM:") $i="SM:"sm } print}' \
            rg_line.txt \
            > fixed_rg_line.txt
        cat fixed_rg_line.txt
        sed -n '/@RG/q;p' header.txt > first_half.txt
        sed -n '/@RG/,$p' header.txt | sed '1d' > second_half.txt

        cat first_half.txt fixed_rg_line.txt second_half.txt > fixed_header.txt

        date
        samtools reheader fixed_header.txt ~{uBAM} > "~{prefix}.RH.bam"
        date
    >>>

    output {
        File reheadered_bam = "~{prefix}.RH.bam"
    }

    runtime {
        cpu:            2
        memory:         "8 GiB"
        disks:          "local-disk ~{ssd_sz} LOCAL"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}
