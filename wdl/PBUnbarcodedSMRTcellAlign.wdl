version 1.0

import "tasks/utils/AlignAndCheckFingerprintCCS.wdl" as major

import "tasks/Utils.wdl"
import "tasks/PBUtils.wdl"
import "tasks/utils/BAMutils.wdl"
import "tasks/utils/GeneralUtils.wdl" as GU
import "tasks/AlignedMetrics.wdl"

import "tasks/Finalize.wdl" as FF

workflow PBUnbarcodedSMRTcellAlign {
    meta {
        description: "Performs alignment of analysis-ready (CCSed, CLR) BAM."
    }
    input {
        File uBAM
        File uPBI
        String movie

        String sample
        String participant_name

        String? fingerprint_store
        String? sample_id_at_store
        Boolean turn_off_fingperprint_check = false
        File ref_map_file

        Int? num_shards

        String gcs_out_root_dir
    }

    String workflow_name = "PBUnbarcodedSMRTcellAlign"

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{participant_name}/~{workflow_name}"

    String output_prefix = sample + "." + movie

    ###################################################################################
    if (! turn_off_fingperprint_check) {
        if (!defined(fingerprint_store) || !(defined(sample_id_at_store))) {
            call Utils.StopWorkflow {input: reason = "Fingerprint information not provided."}
        }
    }

    # update header so as to be conformant with other input formats
    call Utils.ComputeAllowedLocalSSD {input: intended_gb = ceil(2 * size(uBAM, "GiB"))}
    call ReheaderInputUnalignedBAM as RH {input: uBAM = uBAM, sample_name = sample, ssd_cnt = ComputeAllowedLocalSSD.numb_of_local_ssd}
    call PBUtils.PBIndex as Index {input: bam = RH.reheadered_bam}

    call BAMutils.GetReadGroupInfo as RG {input: uBAM = RH.reheadered_bam, keys = ['SM', 'LB']}

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
    ##########
    # store the results into designated bucket
    ##########
    String cdir = outdir

    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = outdir, file = AlignAndCheckFingerprintCCS.aligned_bam, name = output_prefix + '.bam' }
    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = outdir, file = AlignAndCheckFingerprintCCS.aligned_bai, name = output_prefix + '.bam.bai' }
    call FF.FinalizeToFile as FinalizeAlignedPbi { input: outdir = outdir, file = AlignAndCheckFingerprintCCS.aligned_pbi, name = output_prefix + '.bam.pbi' }

    call FF.FinalizeToFile as FinalizeAlnMetrics { input: outdir = outdir, file = AlignAndCheckFingerprintCCS.alignment_metrics_tar_gz, name = output_prefix + '.aln_metrics.tar.gz' }

    if (! turn_off_fingperprint_check) {
        call FF.FinalizeToFile as FinalizeFPDetails  { input: outdir = outdir, file = select_first([AlignAndCheckFingerprintCCS.fingerprint_detail_tar_gz]), name = output_prefix + '.fp_details.tar.gz' }
    }

    call GU.GetTodayDate as today {}

    output {
        File aligned_bam = FinalizeAlignedBam.gcs_path
        File aligned_bai = FinalizeAlignedBai.gcs_path
        File aligned_pbi = FinalizeAlignedPbi.gcs_path

        String? fingerprint_check_results = AlignAndCheckFingerprintCCS.fp_status
        Float? fingerprint_check_LOD = AlignAndCheckFingerprintCCS.fp_lod_expected_sample
        File? fingerprint_check_tar_gz = FinalizeFPDetails.gcs_path

        Float wgs_cov = MosDepthWGS.wgs_cov
        File alignment_metrics_tar_gz = FinalizeAlnMetrics.gcs_path

        String last_alignment_date = today.yyyy_mm_dd
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
            sed -i "s/$/\tSM:tbd/" rg_line.txt
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

