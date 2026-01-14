version 1.0

import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/FastqUtils.wdl"

import "../../../tasks/Visualization/NanoPlot.wdl" as QC0
import "../../TechAgnostic/Utility/FASTQstats.wdl" as QC1
import "../../TechAgnostic/Utility/DystPeaker.wdl" as QC2

import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../TechAgnostic/Utility/SaveFilesToDestination.wdl" as SAVE

workflow Work {

    meta {
        description:
        "A workflow that merges multiple HiFi fastq files for a sample into a single FASTQ."
    }

    parameter_meta {
        # input
        prefix:
        "prefix for output files"

        short_reads_threshold:
        "A threshold below which the reads will be classified as short, used in read-length related metrics collection workflow."

        run_nanoplot:
        "if true, will kick off Nanoplot to collect metrics on the uBAM; this isn't necessary if you intend to process the reads through our alignment-based pipelines since Nanoplot is run there automatically and produces a superset of metrics"

        gcs_out_root_dir:
        "GCS bucket to store the merged FASTQ and metrics files"

        # output
        nanoplot_u_summ:
        "A few metrics output by Nanoplot on the uBAM"
        seqkit_stats:
        "A few metrics output by seqkit stats"

        read_len_summaries:
        "A few metrics summarizing the read length distribution"
        read_len_peaks:
        "Peaks of the read length distribution (heruistic)"
        read_len_deciles:
        "Deciles of the read length distribution"

        metrics_files:
        "A map where keys are summary-names and values are paths to files generated from the various QC/metrics tasks"
    }

    input {
        Array[File] ccs_fqs

        String prefix

        Int short_reads_threshold

        Boolean run_nanoplot = false

        String disk_type = "SSD"
        String gcs_out_root_dir
    }

    output {
        String last_processing_date = today.yyyy_mm_dd

        ########################################
        File merged_fq = finalized_merged_fq_path

        ########################################
        Map[String, Float] seqkit_stats = FASTQstats.stats

        # read length metrics
        # File read_len_hist = DystPeaker.read_len_hist
        Array[Int] read_len_peaks = DystPeaker.read_len_peaks
        Array[Int] read_len_deciles = DystPeaker.read_len_deciles
        Map[String, String] read_len_summaries = DystPeaker.read_len_summaries

        Map[String, Float]? nanoplot_u_summ = NanoPlotFromFastqs.stats_map

        # file-based outputs all packed into a finalization map
        Map[String, String] metrics_files = metrics_files_out
    }

    String outdir = gcs_out_root_dir

    #########################################################################################
    # aggregation
    if (length(ccs_fqs) > 1) {
        call FastqUtils.MergeFastqs as MergeAllFastqs { input: fastqs = ccs_fqs, prefix = prefix, disk_type = disk_type}
    }
    File ccs_fq  = select_first([ MergeAllFastqs.merged_fastq, ccs_fqs[0] ])

    # save files
    String dummy = basename(ccs_fq)
    String dummy_b = sub(dummy, ".gz$", "")
    String sample_outdir = sub(outdir, "/+$", "") + '/sample_fastq'
    if (dummy != dummy_b) {
        call FF.FinalizeToFile as FinalizeMergedFQ { input: outdir = sample_outdir, file = ccs_fq, name = prefix + ".fq.gz" }
    }
    if (dummy == dummy_b) {
        call FF.CompressAndFinalize as CompressAndFinalizeMergedFQ { input: outdir = sample_outdir, file = ccs_fq, name = prefix + ".fq.gz" }
    }
    String finalized_merged_fq_path = select_first([FinalizeMergedFQ.gcs_path, CompressAndFinalizeMergedFQ.gcs_path])

    ###########################################################
    # more QCs and metrics
    ##########
    if (run_nanoplot) {
        call QC0.NanoPlotFromFastqs { input: fastqs = [ccs_fq] , is_ont_rich_fastq = false }
        FinalizationManifestLine a = object
                                    {files_to_save: flatten([[NanoPlotFromFastqs.stats], NanoPlotFromFastqs.plots]),
                                    is_singleton_file: false,
                                    destination: sub(outdir, "/+$", "") + "/nanoplot",
                                    output_attribute_name: "nanoplot"}
    }
    ##########
    call QC1.FASTQstats { input: reads=ccs_fq, file_type='FASTQ' }
    ##########
    call QC2.DystPeaker { input:
        input_file=ccs_fq, input_is_bam=false,
        id=prefix,
        short_reads_threshold=short_reads_threshold,
        gcs_out_root_dir=outdir
    }
    Map[String, String] read_len_out = {"read_len_hist": DystPeaker.read_len_hist,
                                        "raw_rl_file": DystPeaker.read_len_summaries['raw_rl_file']}
    ##########
    if (run_nanoplot) {
        call SAVE.SaveFilestoDestination as SaveRest { input:
            instructions = select_all([a]),
            already_finalized = [read_len_out]
        }
    }
    Map[String, String] metrics_files_out = select_first([SaveRest.result, read_len_out])

    ###########################################################
    call GU.GetTodayDate as today {}
}
