version 1.0

import "tasks/CCSLima.wdl"
import "tasks/SMRTtools.wdl"
import "tasks/metrics/CollectSMRTCellUnalignedMetrics.wdl" as uBAMCustomMetrics
import "tasks/PBUtils.wdl" as PB
import "tasks/Finalize.wdl" as FF
import "tasks/Utils.wdl"
import "tasks/utils/GeneralUtils.wdl" as GU

workflow PreprocessDemultiplexedCCSedSMRTCell {
    meta {
        desciption:
        "A workflow for preprocessing on-board demultiplxed SMRTCell."
    }
    input {
        String smrtcell_data_dir
        String library
        String barcode_names
        String downstream_sample_ids

        File barcodes_fasta
        String barcode_design
        File? lima_options_file

        String gcs_out_root_dir
    }

    parameter_meta {
        # inputs
        smrtcell_data_dir: "Bucket location of the 'folder' that contains the CCSed BAM, PBI, XMLs and other companion files"
        barcodes_fasta:    "FASTA containing the barcode names and their sequences."
        barcode_design:    "Design employed when barcoding the DNA molecules (https://lima.how/barcode-design.html)."
        lima_options_file: "A 2-col TSV file with lima option long name in the 1st col, and value in the 2nd col. For flag-like options, the 2nd column should be set to 'true'. Omit if default values are desired."
        library:           "Formatted as run_id|well"

        barcode_names:         "Names of barcode for each sample on-board this SMRTCell"
        downstream_sample_ids: "These sample names will be encoded into the SM field of the read group line of the BAM header, and used by downstream applications and users. Note: it's assumed the order here matches the order given in barcode_names."
        demultiplex_with_hifi_reads_only: "If true, will filter to Hifi reads first (if it hasn't been done on-instrument) before demultiplexing."

        # outputs
        hifi_bam:  "Available if the input bam is CCS-ed only, but not Hifi filtered"
        unbarcoded_files: "GCS path to the dir holding unbarcoded BAM and other auxiliary files."
        all_barcode_2_folder: "A 2-col TSV file where the 1st col holds lima-discovered barcode names and 2nd col holds GCS path to the dir holding corresponding demultiplexed BAM and auxiliary files."
        selected_barcode_2_folder: "Same format as all_barcode_2_folder, except here the 1st column contains only the barcodes given in barcode_names."

        smrtcell_pbi_summary:  "A custom summary (2-col TSV) of the PBI accompanying the raw input bam."
        runqc_reports_tar_gz:   "tar.gz file of results from running SMRTtools runqc-reports"
        barcode_report_tar_gz:  "tar.gz file of results from running SMRTtools barcode_report"
        all_lima_metris_tar_gz: "tar.gz file of metric files produced by lima itself, and custom reports from parsing lima.report."
    }

    String workflow_name = "PreprocessDemultiplexedCCSedSMRTCell"

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name + "/" + movie
    String outdir_metrics = outdir + "/metrics"

    call Utils.SplitDelimitedString as get_barcodes {input: s = barcode_names, sep = ','}
    call Utils.SplitDelimitedString as get_sample_ids {input: s = downstream_sample_ids, sep = ','}
    if (length(get_barcodes.arr) != length(get_sample_ids.arr)) {
        call Utils.StopWorkflow as unmatched_barcodes_and_samples {
            input: reason = "Length of barcode names array and sample ids array don't match."
        }
    }
    call Utils.ConstructMap {input: keys = get_barcodes.arr, values = get_sample_ids.arr}
    Map[String, String] barcode_2_sample = ConstructMap.converted

    call DistinguishFiles as prep { input: smrtcell_data_dir = smrtcell_data_dir }

    String movie = basename(basename(basename(prep.bam_path, ".reads.bam"), ".hifi_reads.bam"), ".unbarcoded")

    ###################################################################################
    call PB.SummarizeCCSReport { input: report = prep.ccs_report_txt_path }
    call uBAMCustomMetrics.CollectSMRTCellUnalignedMetrics {input: smrtcell_pbi = prep.pbi_path}
    call FF.FinalizeToFile as FinalizePBISummary {
        input:
            file = CollectSMRTCellUnalignedMetrics.pbi_summary,
            outdir = outdir_metrics
    }

    if (! prep.is_hifi) {

        call Utils.StopWorkflow {input: reason = "On-board demultiplexed BAMs are expected to be HiFi bams. That assumption is broken here."}

        call PB.ExtractHifiReads {
            input:
                bam         = prep.bam_path,
                sample_name = downstream_sample_ids,
                library     = library,
                prefix      = "~{movie}.hifi_reads"
        }
        call PB.PBIndex { input: bam = ExtractHifiReads.hifi_bam }
        call FF.FinalizeToFile as FinalizeHifiBam {input: file = ExtractHifiReads.hifi_bam, outdir = outdir}
        call FF.FinalizeToFile as FinalizeHifiPbi {input: file = PBIndex.pbi, outdir = outdir}
    }

    ###################################################################################
    # runqc-reports
    call SMRTtools.RunQCreports {
        input: smrtcell_data_dir = smrtcell_data_dir, movie_name = movie
    }
    call FF.FinalizeToFile as FinalizeRunQCreports {
        input:
            file = RunQCreports.runqc_reports_tar_gz,
            outdir = outdir_metrics
    }

    ###################################################################################
    # # demultiplexing
    # if(defined(lima_options_file)) {
    #     call CCSLima.ParseLimaOptionsFile { input: lima_options_file = select_first([lima_options_file]) }
    # }
    # Array[String] custom_options = select_first([ParseLimaOptionsFile.parsed_options, [""]])

    # # demultiplex CCS-ed BAM
    # File bam_to_demultiplex = if (prep.is_hifi || ! demultiplex_with_hifi_reads_only) then prep.bam_path else select_first([ExtractHifiReads.hifi_bam])
    # call CCSLima.CCSLima as Lima {
    #     input:
    #         bam = bam_to_demultiplex,
    #         movie = movie,
    #         barcode_file = barcodes_fasta,
    #         barcode_design = barcode_design,
    #         custom_options = custom_options
    # }

    # # un-barcoded
    # call FF.FinalizeToDir as FinalizeUnbarcoded {
    #     input:
    #         files = Lima.unbarcoded_files,
    #         outdir = outdir + "/unbarcoded"
    # }

    call GatherAutoDetectedBarcodes { input: smrtcell_data_dir = smrtcell_data_dir}

    call AllExpectedBarcodesFound {input: all_autodetected_barcodes=GatherAutoDetectedBarcodes.barcodes, expected_barcode_names=get_barcodes.arr}
    if (!AllExpectedBarcodesFound.indeed) {
        call Utils.StopWorkflow as missing_barcodes {input: reason = "Not all barcodes expected by user are detected by lima."}
    }

    call GatherUnexpectedBarcodes {input: expected_barcodes = get_barcodes.arr, auto_detected_barcodes = GatherAutoDetectedBarcodes.barcodes}

    call SelectExpectedBarcodesOnly {
        input: expected_barcode_names = get_barcodes.arr, all_barcode_2_dir = GatherAutoDetectedBarcodes.barcode_2_dir_tsv
    }
    call FF.FinalizeToFile as FinalizeAllBarcodeMetaFiles      {input: file = GatherAutoDetectedBarcodes.barcode_2_dir_tsv, outdir = outdir, name = "all_barcode_2_folder.tsv"}
    call FF.FinalizeToFile as FinalizeSelectedBarcodeMetaFiles {input: file = SelectExpectedBarcodesOnly.barcode_2_dir,     outdir = outdir, name = "selected_barcode_2_folder.tsv"}

    # # each barcode in the cell
    # scatter (idx in range(length(Lima.barcodes_list))) {

    #     String bc = Lima.barcodes_list[idx]
    #     call AutodetectedBarcodeIsInExpected as check_bc { input: autodetected_barcode = bc, expected_barcode_names = get_barcodes.arr }
    #     if (check_bc.is_in) {  # lima autodetected barcodes is assumed to be a superset of actually used barcodes for samples on this Cell
    #         call Utils.FixSampleName {
    #             input:
    #                 bam = Lima.demux_bams[idx],
    #                 sample_name = barcode_2_sample[bc]
    #         }
    #     }

    #     File bam_2_finalize = select_first([FixSampleName.reheadered_bam, Lima.demux_bams[idx]])
    #     call FF.FinalizeToDir as FinalizeDemultiplexed {
    #         input:
    #             files = [bam_2_finalize, Lima.demux_pbis[idx], Lima.demux_xmls[idx]],
    #             outdir = outdir + "/" + Lima.barcodes_list[idx]
    #     }
    # }
    # call CCSLima.LocateBarcodeSpecificFolders {
    #     input: barcodes_list = Lima.barcodes_list, finalized_dir_for_each_barcode = FinalizeDemultiplexed.gcs_dir
    # }
    # call SelectExpectedBarcodesOnly {
    #     input: expected_barcode_names = get_barcodes.arr, all_barcode_2_dir = LocateBarcodeSpecificFolders.barcode_2_dir
    # }
    # call FF.FinalizeToFile as FinalizeAllBarcodeMetaFiles      {input: file = LocateBarcodeSpecificFolders.barcode_2_dir, outdir = outdir, name = "all_barcode_2_folder.tsv"}
    # call FF.FinalizeToFile as FinalizeSelectedBarcodeMetaFiles {input: file = SelectExpectedBarcodesOnly.barcode_2_dir,   outdir = outdir, name = "selected_barcode_2_folder.tsv"}

    # metrics output by lima itself
    call GatherLimaFiles {input: smrtcell_data_dir = smrtcell_data_dir, movie = movie}
    if (defined(GatherLimaFiles.lima_json_path_file)) {
        String lima_json_path = read_string(select_first([GatherLimaFiles.lima_json_path_file]))
    }
    if (defined(GatherLimaFiles.lima_log_path_file)) {
        String lima_log_path = read_string(select_first([GatherLimaFiles.lima_log_path_file]))
    }
    if (defined(GatherLimaFiles.lima_guess_path_file)) {
        String lima_guess_path = read_string(select_first([GatherLimaFiles.lima_guess_path_file]))
    }
    call CCSLima.GatherLimaAndCustomDemultiplexingMetrics {
        input:
            lima_summary = GatherLimaFiles.lima_summary_path,
            lima_counts = GatherLimaFiles.lima_counts_path,
            lima_report = GatherLimaFiles.lima_report_path,
            lima_xml = GatherLimaFiles.lima_xml_path,
            lima_guess = lima_guess_path,
            lima_json = lima_json_path,
            barcodes_fasta = barcodes_fasta
    }
    call FF.FinalizeToFile as FinalizeLimaMetrics {
        input:
            file = GatherLimaAndCustomDemultiplexingMetrics.all_metrics_gz,
            outdir = outdir_metrics
    }

    # metrics on demultiplexing output by barcode_report
    call SMRTtools.BarcodeReportWithOnBoardDemuxedFolder {
        input:
            movie_name = movie,
            smrtcell_data_dir = smrtcell_data_dir
    }
    call FF.FinalizeToFile as FinalizeBarcodeReports {
        input: file = BarcodeReportWithOnBoardDemuxedFolder.barcode_report_tar_gz, outdir = outdir_metrics
    }

    call GU.GetTodayDate as today {}

    output {
        # File? hifi_bam = FinalizeHifiBam.gcs_path
        # File? hifi_pbi = FinalizeHifiPbi.gcs_path

        # String unbarcoded_files = FinalizeUnbarcoded.gcs_dir
        File      all_barcode_2_folder = FinalizeAllBarcodeMetaFiles.gcs_path
        File selected_barcode_2_folder = FinalizeSelectedBarcodeMetaFiles.gcs_path

        File barcode_report_tar_gz = FinalizeBarcodeReports.gcs_path
        File all_lima_metris_tar_gz = FinalizeLimaMetrics.gcs_path

        File smrtcell_pbi_summary = FinalizePBISummary.gcs_path
        File runqc_reports_tar_gz = FinalizeRunQCreports.gcs_path

        Float zmws_ccs_input = SummarizeCCSReport.zmws_input
        Float zmws_ccs_pass_filters = SummarizeCCSReport.zmws_pass_filters
        Float zmws_ccs_fail_filters = SummarizeCCSReport.zmws_fail_filters
        Float zmws_ccs_shortcut_filters = SummarizeCCSReport.zmws_shortcut_filters
        Float zmws_ccs_pass_filters_pct = SummarizeCCSReport.zmws_pass_filters_pct
        Float zmws_ccs_fail_filters_pct = SummarizeCCSReport.zmws_fail_filters_pct
        Float zmws_ccs_shortcut_filters_pct = SummarizeCCSReport.zmws_shortcut_filters_pct

        String last_preprocessing_date = today.yyyy_mm_dd
    }
}

task DistinguishFiles {
    meta {
        desciption: "Given a path to an CCSed SMRTCell's data folder mirrored onto the bucket, separate files for different downstream processing"
    }

    input {
        String smrtcell_data_dir
    }

    command <<<
        set -eux

        gsutil ls ~{smrtcell_data_dir} > files.list

        grep -E "consensusreadset.xml$" files.list > main.xml.txt
        grep -E "sts.xml$" files.list > sts.xml.txt
        grep -E "reads.bam$" files.list > bam.txt
        grep -E ".pbi$" files.list > pbi.txt

        grep -E "ccs_reports.txt$" files.list > ccs_report.txt

        if grep -qF 'hifi' bam.txt; then echo "true" > is_hifi.txt; else echo "false" > is_hifi.txt; fi
    >>>

    output {
        String main_xml_path = read_string("main.xml.txt")
        String bam_path = read_string("bam.txt")
        String pbi_path = read_string("pbi.txt")
        String ccs_report_txt_path = read_string("ccs_report.txt")

        Boolean is_hifi = read_boolean("is_hifi.txt")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task AutodetectedBarcodeIsInExpected {
    input {
        String autodetected_barcode
        Array[String] expected_barcode_names
    }

    command <<<
        set -eux

        if grep -q "^~{autodetected_barcode}$" ~{write_lines(expected_barcode_names)}; then
            echo "true" > res.txt
        else
            echo "false" > res.txt
        fi
    >>>

    output {
        Boolean is_in = read_boolean("res.txt")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task AllExpectedBarcodesFound {
    meta {
        desciption:
        "Check that lima detected barcode names is a superset of user expected barcode names."
    }
    input {
        Array[String] all_autodetected_barcodes
        Array[String] expected_barcode_names
    }

    command <<<
        set -eux

        cp ~{write_lines(all_autodetected_barcodes)} auto.txt
        cp ~{write_lines(expected_barcode_names)} expected.txt

        cat auto.txt
        cat expected.txt

        touch should_be_empty.txt
        comm -13 <(sort auto.txt | uniq) <(sort expected.txt | uniq) > should_be_empty.txt
        if [[ -s should_be_empty.txt ]]; then echo "false" > "res.txt"; else echo "true" > "res.txt"; fi
    >>>

    output {
        Boolean indeed = read_boolean("res.txt")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task SelectExpectedBarcodesOnly {
    meta {
        desciption:
        "Select only barcodes that are expected by user"
    }

    input {
        Array[String] expected_barcode_names
        File all_barcode_2_dir
    }

    command <<<
        set -eux
        touch "selected_barcode_2_folder.tsv"
        for bc in ~{sep=' ' expected_barcode_names}; do
            grep -w "^${bc}" ~{all_barcode_2_dir} >> "selected_barcode_2_folder.tsv"
        done
    >>>

    output {
        File barcode_2_dir = "selected_barcode_2_folder.tsv"
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

# todo: this assumes summetry barcode design
task GatherAutoDetectedBarcodes {
    input {
        String smrtcell_data_dir
    }

    String formatted_gcs_dir = sub(smrtcell_data_dir, "/$", "") + "/"

    command <<<
        set -eux

        gsutil ls ~{formatted_gcs_dir} > everthing.list
        cat everthing.list
        n=$(echo ~{formatted_gcs_dir} | awk -F '/' '{print NF}')
        awk -F '/' -v parent_ncol=$n '{if(NF-parent_ncol==1) print $(NF-1)}' everthing.list | \
            grep "^bc" | \
            awk -F '--' '{print $1}' \
            > auto_detected_barcodes.txt
        cat auto_detected_barcodes.txt

        while IFS= read -r bc
        do
            location=$(grep -F "${bc}" everthing.list)
            echo -e "${bc}\t${location}" >> barcode_2_dir.tsv
        done < auto_detected_barcodes.txt
        cat barcode_2_dir.tsv
    >>>

    output {
        Array[String] barcodes = read_lines("auto_detected_barcodes.txt")
        File barcode_2_dir_tsv = "barcode_2_dir.tsv"
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task GatherLimaFiles {
    input {
        String smrtcell_data_dir
        String movie
    }

    String formatted_gcs_dir = sub(smrtcell_data_dir, "/$", "") + "/"

    command <<<
        set -eux
        # touch lima_files.path
        # p=$(grep -F 'consensusreadset.xml' everthing.list | grep -vF 'unbarcoded')
        # echo -e "${xml}\t${p}" >> lima_files.path
        # grep -F 'lima' everthing.list > everthing.lima.list
        # p=$(grep -F 'counts' everthing.lima.list)
        # echo -e "${counts}\t${p}" >> lima_files.path
        # p=$(grep -F 'report' everthing.lima.list)
        # echo -e "${report}\t${p}" >> lima_files.path
        # p=$(grep -F 'summary' everthing.lima.list)
        # echo -e "${summary}\t${p}" >> lima_files.path
        # if grep -qF 'log' everthing.lima.list; then
        #     p=$(grep -F 'log' everthing.lima.list)
        #     echo -e "${log}\t${p}" >> lima_files.path
        # fi
        # if grep -qF 'guess' everthing.lima.list; then
        #     p=$(grep -F 'guess' everthing.lima.list | head -n 1)
        #     echo -e "${guess}\t${p}" >> lima_files.path
        # fi
        # if grep -qF "~{movie}.json" everthing.list; then
        #     p=$(grep -F "~{movie}.json" everthing.list)
        #     echo -e "${josn}\t${p}" >> lima_files.path
        # fi

        gsutil ls ~{formatted_gcs_dir} > everthing.list
        # file sure to exist
        grep -F 'lima' everthing.list > everthing.lima.list
        grep -F 'summary' everthing.lima.list > "lima_summary.path"
        grep -F 'counts' everthing.lima.list > "lima_counts.path"
        grep -F 'report' everthing.lima.list > "lima_report.path"

        grep -F 'consensusreadset.xml' everthing.list | grep -vF 'unbarcoded' > "lima_xml.path"

        # files exist sometimes
        if grep -qF 'log' everthing.lima.list; then
            grep -F 'log' everthing.lima.list > "lima_log.path"
        fi
        if grep -qF 'guess' everthing.lima.list; then
            grep -F 'guess' everthing.lima.list | head -n 1 > "lima_guess.path"
        fi
        if grep -qF "~{movie}.json" everthing.list; then
            grep -qF "~{movie}.json" everthing.list > "lima_json.path"
        fi
    >>>

    output {
        String lima_summary_path = read_string("lima_summary.path")
        String lima_counts_path = read_string("lima_counts.path")
        String lima_report_path = read_string("lima_report.path")
        String lima_xml_path = read_string("lima_xml.path")

        File? lima_log_path_file = "lima_log.path"
        File? lima_guess_path_file = "lima_guess.path"
        File? lima_json_path_file = "lima_json.path"
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task GatherUnexpectedBarcodes {
    input {
        Array[String] expected_barcodes
        Array[String] auto_detected_barcodes
    }

    command <<<
        set -eu

        comm -13 <(sort ~{write_lines(expected_barcodes)}) \
                 <(sort ~{write_lines(auto_detected_barcodes)}) \
            > unexpected_barcodes.txt
    >>>

    output {
        Array[String] unexpected_barcodes = read_lines("unexpected_barcodes.txt")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
