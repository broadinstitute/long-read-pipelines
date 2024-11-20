version 1.0

import "tasks/CCSLima.wdl"
import "tasks/SMRTtools.wdl"
import "tasks/Utility/Utils.wdl" as DeprecatedUtils
import "../tasks/QC/CollectSMRTCellUnalignedMetrics.wdl" as uBAMCustomMetrics
import "../tasks/Utility//PBUtils.wdl" as PB
import "../tasks/Utility//Finalize.wdl" as FF
import "../tasks/Utility/Utils.wdl"
import "../tasks/Utility/GeneralUtils.wdl" as GU

workflow PreprocessBarcodedCCSedSMRTCell {
    meta {
        desciption:
        "A workflow for preprocessing barcoded (potentially multiplexed) SMRTCell. The cell is assumed to be CCS-ed on-instrument, and the whole data folder is mirrored onto a cloud bucket."
    }
    input {
        String smrtcell_data_dir
        String library
        String barcode_names
        String downstream_sample_ids

        Boolean demultiplex_with_hifi_reads_only

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

    call DistinguishFiles as prep { input: smrtcell_data_dir = smrtcell_data_dir }
    String movie = basename(basename(prep.bam_path, ".reads.bam"), ".hifi_reads.bam")
    String workflow_name = "PreprocessBarcodedCCSedSMRTCell"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name + "/" + movie
    String outdir_metrics = outdir + "/metrics"

    call DeprecatedUtils.SplitDelimitedString as get_barcodes {input: s = barcode_names, separate = ','}
    call DeprecatedUtils.SplitDelimitedString as get_sample_ids {input: s = downstream_sample_ids, separate = ','}
    if (length(get_barcodes.arr) != length(get_sample_ids.arr)) {
        call Utils.StopWorkflow as unmatched_barcodes_and_samples {
            input: reason = "Length of barcode names array and sample ids array don't match."
        }
    }
    call ConstructMap {input: keys = get_barcodes.arr, values = get_sample_ids.arr}
    Map[String, String] barcode_2_sample = ConstructMap.converted

    ###################################################################################
    call PB.SummarizeCCSReport { input: report = prep.ccs_report_txt_path }
    call uBAMCustomMetrics.CollectSMRTCellUnalignedMetrics {input: smrtcell_pbi = prep.pbi_path}
    call FF.FinalizeToFile as FinalizePBISummary {
        input:
            file = CollectSMRTCellUnalignedMetrics.pbi_summary,
            outdir = outdir_metrics
    }

    if (! prep.is_hifi) {

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
    # demultiplexing
    if(defined(lima_options_file)) {
        call CCSLima.ParseLimaOptionsFile { input: lima_options_file = select_first([lima_options_file]) }
    }
    Array[String] custom_options = select_first([ParseLimaOptionsFile.parsed_options, [""]])

    # demultiplex CCS-ed BAM
    File bam_to_demultiplex = if (prep.is_hifi || ! demultiplex_with_hifi_reads_only) then prep.bam_path else select_first([ExtractHifiReads.hifi_bam])
    call CCSLima.CCSLima as Lima {
        input:
            bam = bam_to_demultiplex,
            movie = movie,
            barcode_file = barcodes_fasta,
            barcode_design = barcode_design,
            custom_options = custom_options
    }

    # un-barcoded
    call FF.FinalizeToDir as FinalizeUnbarcoded {
        input:
            files = Lima.unbarcoded_files,
            outdir = outdir + "/unbarcoded"
    }

    call AllExpectedBarcodesFound {input: all_autodetected_barcodes=Lima.barcodes_list, expected_barcode_names=get_barcodes.arr}
    if (!AllExpectedBarcodesFound.indeed) {
        call Utils.StopWorkflow as missing_barcodes {input: reason = "Not all barcodes expected by user are detected by lima."}
    }

    # each barcode in the cell
    scatter (idx in range(length(Lima.barcodes_list))) {

        String bc = Lima.barcodes_list[idx]
        call AutodetectedBarcodeIsInExpected as check_bc { input: autodetected_barcode = bc, expected_barcode_names = get_barcodes.arr }
        if (check_bc.is_in) {  # lima autodetected barcodes is assumed to be a superset of actually used barcodes for samples on this Cell
            call FixSampleName {
                input:
                    bam = Lima.demux_bams[idx],
                    sample_name = barcode_2_sample[bc]
            }
        }

        File bam_2_finalize = select_first([FixSampleName.reheadered_bam, Lima.demux_bams[idx]])
        call FF.FinalizeToDir as FinalizeDemultiplexed {
            input:
                files = [bam_2_finalize, Lima.demux_pbis[idx], Lima.demux_xmls[idx]],
                outdir = outdir + "/" + Lima.barcodes_list[idx]
        }
    }
    call CCSLima.LocateBarcodeSpecificFolders {
        input: barcodes_list = Lima.barcodes_list, finalized_dir_for_each_barcode = FinalizeDemultiplexed.gcs_dir
    }
    call SelectExpectedBarcodesOnly {
        input: expected_barcode_names = get_barcodes.arr, all_barcode_2_dir = LocateBarcodeSpecificFolders.barcode_2_dir
    }
    call FF.FinalizeToFile as FinalizeAllBarcodeMetaFiles      {input: file = LocateBarcodeSpecificFolders.barcode_2_dir, outdir = outdir, name = "all_barcode_2_folder.tsv"}
    call FF.FinalizeToFile as FinalizeSelectedBarcodeMetaFiles {input: file = SelectExpectedBarcodesOnly.barcode_2_dir,   outdir = outdir, name = "selected_barcode_2_folder.tsv"}

    # metrics output by lima itself
    call CCSLima.GatherLimaAndCustomDemultiplexingMetrics {
        input:
            lima_counts = Lima.lima_counts,
            lima_guess = Lima.lima_guess,
            lima_report = Lima.lima_report,
            lima_summary = Lima.lima_summary,
            lima_json = Lima.lima_json,
            lima_xml = Lima.lima_xml,
            barcodes_fasta = barcodes_fasta
    }
    call FF.FinalizeToFile as FinalizeLimaMetrics {
        input:
            file = GatherLimaAndCustomDemultiplexingMetrics.all_metrics_gz,
            outdir = outdir_metrics
    }

    # metrics on demultiplexing output by barcode_report
    call SMRTtools.BarcodeReport {
        input:
            movie_name = movie,
            all_lima_files = flatten([select_all([Lima.lima_counts, Lima.lima_guess, Lima.lima_report, Lima.lima_summary, Lima.lima_json, Lima.lima_xml]), Lima.unbarcoded_files, Lima.demux_bams, Lima.demux_pbis, Lima.demux_xmls]),
            barcodes_fasta = barcodes_fasta
    }
    call FF.FinalizeToFile as FinalizeBarcodeReports {
        input: file = BarcodeReport.barcode_report_tar_gz, outdir = outdir_metrics
    }

    call GU.GetTodayDate as today {}

    output {
        File? hifi_bam = FinalizeHifiBam.gcs_path
        File? hifi_pbi = FinalizeHifiPbi.gcs_path
        String unbarcoded_files = FinalizeUnbarcoded.gcs_dir
        File      all_barcode_2_folder = FinalizeAllBarcodeMetaFiles.gcs_path
        File selected_barcode_2_folder = FinalizeSelectedBarcodeMetaFiles.gcs_path

        File smrtcell_pbi_summary = FinalizePBISummary.gcs_path
        File runqc_reports_tar_gz = FinalizeRunQCreports.gcs_path
        File barcode_report_tar_gz = FinalizeBarcodeReports.gcs_path
        File all_lima_metris_tar_gz = FinalizeLimaMetrics.gcs_path

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

task FixSampleName {

    meta {
        desciption:
        "This fixes the sample name of a demultiplexed BAM"
    }

    parameter_meta {
        bam: {
            localization_optional: true,
            description: "BAM file"
        }
        sample_name: "sample name"
    }

    input {
        File bam
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, ".bam")
    Int disk_size = 3*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools view --no-PG -H ~{bam} > header.txt
        awk '$1 ~ /^@RG/' header.txt > rg_line.txt
        if ! grep -qF "SM:" rg_line.txt; then
            sed -i "s/$/SM:tbd/" rg_line.txt
        fi
        awk -v lib="~{sample_name}" -F '\t' 'BEGIN {OFS="\t"} { for (i=1; i<=NF; ++i) { if ($i ~ "SM:") $i="SM:"lib } print}' \
            rg_line.txt \
            > fixed_rg_line.txt

        sed -n '/@RG/q;p' header.txt > first_half.txt
        sed -n '/@RG/,$p' header.txt | sed '1d' > second_half.txt

        cat first_half.txt fixed_rg_line.txt second_half.txt > fixed_header.txt
        cat fixed_header.txt

        mv ~{bam} old.bam
        date
        samtools reheader fixed_header.txt old.bam > ~{prefix}.bam
        date
    >>>

    output {
        File reheadered_bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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

task ConstructMap {

    meta {
        desciption:
        "Use only when the keys are guaranteed to be unique and the two arrays are of the same length."
    }

    parameter_meta {
        keys: "The keys of the map"
        values: "The values of the map"
    }

    input {
        Array[String] keys
        Array[String] values
    }
    command <<<
        set -eux
        paste ~{write_lines(keys)} ~{write_lines(values)} > converted.tsv
        cat converted.tsv
    >>>

    output {
        Map[String, String] converted = read_map("converted.tsv")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
