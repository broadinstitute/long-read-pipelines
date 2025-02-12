version 1.0

import "tasks/Utility/Utils.wdl" as DeprecatedUtils
import "../tasks/Alignment/AlignAndCheckFingerprintCCS.wdl" as major
import "../tasks/Utility/Utils.wdl"
import "../tasks/Utility/GeneralUtils.wdl" as GU
import "../tasks/Utility/Finalize.wdl" as FF

workflow PostprocessCCSedDemultiplexedSMRTCell {
    meta {
        desciption:
        "A workflow for postprocessing files output by PreprocessBarcodedCCSedSMRTCell. Note: the various IDs are assumed to be in-phase with the barcode names."
    }

    input {
        String smrtcell_id
        String movie_name

        File demuxed_barcode_2_folder

        String barcode_names
        String biosample_ids
        String sample_ids
        String downstream_sample_ids

        String fingerprint_store
        String sample_ids_at_store
        Boolean turn_off_fingperprint_check = false

        File ref_map_file

        String gcs_out_root_dir
    }

    parameter_meta {
        # inputs
        demuxed_barcode_2_folder: "A 2-col TSV file, where the 1st col holds the barcodes expected in the SMRTCell, and 2nd col holds the corresponding folder hosting the demultiplexed cells."

        barcode_names:         "Names of barcode for each sample on-board this SMRTCell, delimited by comma. Note: the various IDs are assumed to be in-phase with the barcode names."
        biosample_ids:         "Names of samples at biology level, e.g. maps to a human, delimited by comma"
        sample_ids:            "Aliquot IDs, delimited by comma"
        downstream_sample_ids: "These sample names will be encoded into the SM field of the read group line of the BAM header, and used by downstream applications and users. Note: it's assumed the order here matches the order given in barcode_names."

        fingerprint_store:   "Bucket name and prefix (gs://...) storing the fingerprinting VCFs"
        sample_ids_at_store: "A comma delimited string holding samples ID at fingerprint store"
        turn_off_fingperprint_check: "Please turn of fingerprint check if the reference is not GRCh38."

        # outputs
        fingerprint_check_results: "Fingerprint check results summary for each barcode/sample, a string delimited by comma. In phase with input barcodes/samples. Each unit can only be one of [PASS, FAIL, BORDERLINE]."
        fingerprint_check_LODs: "Fingerprint check LODs for each bacode/sample, delimited by comma. In phase with input barcodes/samples."

        barcode_2_aln_dir: "A 2-col TSV file, where the 1st col holds barcode names, and 2nd col holds the corresponding dir holding alignment files."

        barcode_2_aln_metric_tar_gz: "A 2-col TSV file, where the 1st col holds barcode names, and 2nd col holds the corresponding alignment metric files (tag.gz-ed)."
        barcode_2_fpcheck_tar_gz: "A 2-col TSV file, where the 1st col holds barcode names, and 2nd col holds the corresponding fingerprint check files (tar.gz-ed)."
    }

    ###################################################################################
    # prep work

    # where to store final results
    String workflow_name = "PostprocessCCSedDemultiplexedSMRTCell"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name
    String outdir_aln     = outdir + '/alignments'
    String outdir_metrics = outdir + '/metrics'

    # associate barcode to various forms of sample ids
    Array[Pair[String, String]] barcode_2_folder = read_map(demuxed_barcode_2_folder)  # type coercion Map to Array[Pair] may only work in WDL 1.0
    call DeprecatedUtils.SplitDelimitedString as get_barcodes {input: s = barcode_names, separate = ','}
    call DeprecatedUtils.SplitDelimitedString as get_biosamples {input: s = biosample_ids, separate = ','}
    call DeprecatedUtils.SplitDelimitedString as get_aliquots {input: s = sample_ids, separate = ','}
    call DeprecatedUtils.SplitDelimitedString as get_ds_ids {input: s = downstream_sample_ids, separate = ','}
    call DeprecatedUtils.SplitDelimitedString as get_fp_ids {input: s = sample_ids_at_store, separate = ','}
    if (length(barcode_2_folder) != length(get_fp_ids.arr)) {
        call Utils.StopWorkflow as unmatched_barcodes_and_samples {
            input: reason = "Length of barcode names array and sample ids array don't match."
        }
    }
    call ConstructBarcodeToIDs {
        input:
            barcode_names = get_barcodes.arr,
            biosample_ids = get_biosamples.arr,
            sample_ids = get_aliquots.arr,
            downstream_sample_ids = get_ds_ids.arr,
            sample_ids_at_store = get_fp_ids.arr
    }

    ###################################################################################
    # major deal
    scatter (bc_n_dir in barcode_2_folder) { # each barcode, i.e. each sample loaded on the SMRTCell

        String bc = bc_n_dir.left

        call GetDemxedFilePaths {input: demux_dir = bc_n_dir.right}

        call GetReadGroupInfo as RG {input: uBAM = GetDemxedFilePaths.bam_path, keys = ['SM', 'LB']}

        call major.AlignAndCheckFingerprintCCS {
            input:
                uBAM = GetDemxedFilePaths.bam_path,
                uPBI = GetDemxedFilePaths.pbi_path,
                bam_sample_name = RG.read_group_info['SM'],
                library = RG.read_group_info['LB'],

                turn_off_fingperprint_check = turn_off_fingperprint_check,
                fp_store = fingerprint_store,
                sample_id_at_store = ConstructBarcodeToIDs.barcode_2_fp_id[bc],
                ref_map_file = ref_map_file
        }

        ###################################################################################
        # finalize each barcode
        String bc_specific_aln_out    = outdir + '/alignments/' + smrtcell_id + '/' + bc
        String bc_specific_metric_out = outdir + "/metrics/"    + smrtcell_id + '/' + bc

        call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = bc_specific_aln_out, file = AlignAndCheckFingerprintCCS.aligned_bam, name = movie_name + '.' + bc + '.bam' }
        call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = bc_specific_aln_out, file = AlignAndCheckFingerprintCCS.aligned_bai, name = movie_name + '.' + bc + '.bai' }
        call FF.FinalizeToFile as FinalizeAlignedPbi { input: outdir = bc_specific_aln_out, file = AlignAndCheckFingerprintCCS.aligned_pbi, name = movie_name + '.' + bc + '.pbi' }

        call FF.FinalizeToFile as FinalizeAlnMetrics { input: outdir = bc_specific_metric_out, file = AlignAndCheckFingerprintCCS.alignment_metrics_tar_gz }

        if (! turn_off_fingperprint_check) {
            call FF.FinalizeToFile as FinalizeFPDetails  { input: outdir = bc_specific_metric_out, file = select_first([AlignAndCheckFingerprintCCS.fingerprint_detail_tar_gz]) }
        }
    }

    ###################################################################################
    # For Terra data tables
    call LocateBarcodeSpecificFoldersOrFiles as bc_2_aln_dir {input: barcodes_list = bc, finalized_dir_or_file_for_each_barcode = bc_specific_aln_out}
    call FF.FinalizeToFile as finalize_bc_2_aln_dir_tsv {input: file = bc_2_aln_dir.barcode_2_gs_path, outdir = outdir_aln + '/' + smrtcell_id }

    call LocateBarcodeSpecificFoldersOrFiles as bc_2_aln_metric_dir {input: barcodes_list = bc, finalized_dir_or_file_for_each_barcode = FinalizeAlnMetrics.gcs_path}
    call FF.FinalizeToFile as finalize_bc_2_aln_metrics_tsv {input: file = bc_2_aln_metric_dir.barcode_2_gs_path, outdir = outdir_metrics + '/' + smrtcell_id }
    if (! turn_off_fingperprint_check) {
        call LocateBarcodeSpecificFoldersOrFiles as bc_2_fp_details_dir {input: barcodes_list = bc, finalized_dir_or_file_for_each_barcode = select_all(FinalizeFPDetails.gcs_path)}
        call FF.FinalizeToFile as finalize_bc_2_fp_details_tsv {input: file = bc_2_fp_details_dir.barcode_2_gs_path, outdir = outdir_metrics + '/' + smrtcell_id }
    }

    Array[String?] fp_res  = AlignAndCheckFingerprintCCS.fp_status
    Array[Float?]  fp_lods = AlignAndCheckFingerprintCCS.fp_lod_expected_sample

    call GU.GetTodayDate as today {}

    output {
        File barcode_2_aln_dir = finalize_bc_2_aln_dir_tsv.gcs_path

        Array[String?] fingerprint_check_results = fp_res
        Array[Float?] fingerprint_check_LODs = fp_lods

        File barcode_2_aln_metric_tar_gz = finalize_bc_2_aln_metrics_tsv.gcs_path
        File? barcode_2_fpcheck_tar_gz = finalize_bc_2_fp_details_tsv.gcs_path

        String last_postprocessing_date = today.yyyy_mm_dd
    }
}

###################################################################################
# todo: these files are barcode-related, gather them in a single WDL

task GetDemxedFilePaths {
    meta {
        desciption:
        "Given a dir hoding demultiplexed files (i.e. bam, pbi, and xml), get their cloud paths."
    }

    input {
        String demux_dir
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        gsutil ls ~{demux_dir} > files.list
        cat files.list

        grep "consensusreadset.xml" files.list > xml.gspath
        grep ".bam" files.list | grep -vF '.pbi' > bam.gspath
        grep ".pbi" files.list > pbi.gspath
    >>>

    output {
        String xml_path = read_string("xml.gspath")
        String bam_path = read_string("bam.gspath")
        String pbi_path = read_string("pbi.gspath")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task ConstructBarcodeToIDs {
    meta {
        desciption:
        "Construct maps from barcode to various forms of sample ids. IMPORTANT: the input arrays are all assumed to be in-phase."
    }

    input {
        Array[String] barcode_names
        Array[String] biosample_ids
        Array[String] sample_ids
        Array[String] downstream_sample_ids
        Array[String] sample_ids_at_store
    }

    command <<<
        set -eux

        paste ~{write_lines(barcode_names)} \
              ~{write_lines(biosample_ids)} \
              > "bc2bio.tsv"

        paste ~{write_lines(barcode_names)} \
              ~{write_lines(sample_ids)} \
              > "bc2aliquot.tsv"

        paste ~{write_lines(barcode_names)} \
              ~{write_lines(downstream_sample_ids)} \
              > "bc2downstream.tsv"

        paste ~{write_lines(barcode_names)} \
              ~{write_lines(sample_ids_at_store)} \
              > "bc2fp.tsv"
    >>>

    output {
        Map[String, String] barcode_2_bio_id = read_map("bc2bio.tsv")
        Map[String, String] barcode_2_aliquot_id = read_map("bc2aliquot.tsv")
        Map[String, String] barcode_2_downstream_id = read_map("bc2downstream.tsv")
        Map[String, String] barcode_2_fp_id = read_map("bc2fp.tsv")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task LocateBarcodeSpecificFoldersOrFiles {
    meta {
        desciption: "Generates a 2-col TSV where the 1st col is the barcode name, and 2nd is the GCS \'folder\' or file for that barcode."
    }
    input {
        Array[String] barcodes_list
        Array[String] finalized_dir_or_file_for_each_barcode
    }

    command <<<
        set -eux
        paste <(cat ~{write_lines(barcodes_list)}) \
              <(cat ~{write_lines(finalized_dir_or_file_for_each_barcode)}) \
              > barcode_2_gs_path.tsv
    >>>

    output {
        File barcode_2_gs_path = "barcode_2_gs_path.tsv"
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task GetReadGroupInfo {
    meta {
        desciption:
        "Get some read group information Given a single-readgroup BAM. Will fail if the information isn't present."
    }

    parameter_meta {
        uBAM: "The input BAM file."
        keys: "A list of requested fields in the RG line, e.g. ID, SM, LB."
    }

    input {
        String uBAM  # not using file as call-caching brings not much benefit

        Array[String] keys
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{uBAM} | grep "^@RG" | tr '\t' '\n' > rh_header.txt

        for attribute in ~{sep=' ' keys}; do
            value=$(grep "^${attribute}" rh_header.txt | awk -F ':' '{print $2}')
            echo -e "${attribute}\t${value}" >> "result.txt"
        done
    >>>

    output {
        Map[String, String] read_group_info = read_map("result.txt")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}
