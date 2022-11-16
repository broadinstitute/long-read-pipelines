version 1.0

import "tasks/utils/BAMutils.wdl"
import "tasks/Utils.wdl"
import "tasks/utils/GeneralUtils.wdl" as GU
import "tasks/Finalize.wdl" as FF

workflow HiFiBamToFastQ {
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

        String gcs_out_root_dir_ref_free
    }

    ###################################################################################
    # prep work

    # where to store final results
    String workflow_name = "PostprocessCCSedDemultiplexedSMRTCell"
    String outdir = sub(gcs_out_root_dir_ref_free, "/$", "") + "/" + workflow_name
    String outdir_ref_free = outdir + '/RefFree'
    # String outdir_aln     = outdir + '/alignments'
    # String outdir_metrics = outdir + '/metrics'

    # associate barcode to various forms of sample ids
    Array[Pair[String, String]] barcode_2_folder = read_map(demuxed_barcode_2_folder)  # type coercion Map to Array[Pair] may only work in WDL 1.0
    call Utils.SplitDelimitedString as get_barcodes {input: s = barcode_names, sep = ','}
    call Utils.SplitDelimitedString as get_biosamples {input: s = biosample_ids, sep = ','}
    call Utils.SplitDelimitedString as get_aliquots {input: s = sample_ids, sep = ','}
    call Utils.SplitDelimitedString as get_ds_ids {input: s = downstream_sample_ids, sep = ','}
    call Utils.SplitDelimitedString as get_fp_ids {input: s = sample_ids_at_store, sep = ','}
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

        call InputBarcodeWasNotDetected {input: bacode_and_demux_dir = bc_n_dir}
        if (!InputBarcodeWasNotDetected.res) {

            call GetDemxedFilePaths {input: demux_dir = bc_n_dir.right}

            # call BAMutils.GetReadGroupInfo as RG {input: uBAM = GetDemxedFilePaths.bam_path, keys = ['SM', 'LB']}

            call Utils.BamToFastq {input: bam = GetDemxedFilePaths.bam_path, prefix = "does_not_matter"}

            ###################################################################################
            # finalize each barcode
            # String bc_specific_aln_out    = outdir + '/alignments/' + smrtcell_id + '/' + bc
            # String bc_specific_metric_out = outdir + "/metrics/"    + smrtcell_id + '/' + bc

            # call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = bc_specific_aln_out, file = AlignAndCheckFingerprintCCS.aligned_bam, name = movie_name + '.' + bc + '.bam' }
            # call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = bc_specific_aln_out, file = AlignAndCheckFingerprintCCS.aligned_bai, name = movie_name + '.' + bc + '.bai' }
            # call FF.FinalizeToFile as FinalizeAlignedPbi { input: outdir = bc_specific_aln_out, file = AlignAndCheckFingerprintCCS.aligned_pbi, name = movie_name + '.' + bc + '.pbi' }


            String bc_specific_fastq_out = outdir_ref_free + '/' + smrtcell_id + '/' + bc
            call FF.FinalizeToFile as FinalizeFQ { input: outdir = bc_specific_fastq_out, file = BamToFastq.reads_fq, name =  movie_name + '.' + bc + '.hifi.fq.gz' }
        }

        String blah = select_first([FinalizeFQ.gcs_path, 'NA'])
    }

    ###################################################################################
    # For Terra data tables
    call LocateBarcodeSpecificFoldersOrFiles as bc_2_hifi_fq {input: barcodes_list = bc, finalized_dir_or_file_for_each_barcode = blah}
    call FF.FinalizeToFile as finalize_bc_2_hifi_fq_tsv {input: file = bc_2_hifi_fq.barcode_2_gs_path, outdir = outdir_ref_free + '/' + smrtcell_id }

    output {
        File barcode_2_hifi_fq = finalize_bc_2_hifi_fq_tsv.gcs_path
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

task InputBarcodeWasNotDetected {
    input {
        Pair[String, String] bacode_and_demux_dir
    }

    command <<<
        set -eux
        if [[ ~{bacode_and_demux_dir.right} == "NotDetected" ]]; then echo "true"; else echo "false"; fi
    >>>

    output {
        Boolean res = read_boolean(stdout())
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
