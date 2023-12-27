version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/BAMutils.wdl"     as BU
import "../../../tasks/Utility/GeneralUtils.wdl" as GU

import "../../../tasks/Visualization/NanoPlot.wdl"    as QC0
import "../../TechAgnostic/Utility/CountTheBeans.wdl" as QC1
import "../../TechAgnostic/Utility/DystPeaker.wdl"    as QC2
import "../../TechAgnostic/Utility/FASTQstats.wdl"    as QC3

import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../TechAgnostic/Utility/SaveFilesToDestination.wdl" as SAVE

workflow ProcessOnInstrumentDemuxedChunkRefFree {

    meta {
        desciption: "Given an on-instrument demultiplexed hifi_reads.bam, perform ref-independent prep work, and collect a few metrics"
    }

    parameter_meta {
        input_hifi_fq:
        "The preexisting Hifi FASTQ; if provided the metrics calculation is faster; if not, you most likely want to set convert_to_fq to true"

        readgroup_id:
        "Unique ID for a readgroup, for example (D|E)A[0-9]{6}-<barcode>"

        short_reads_threshold:
        "a length threshold below which reads are classified as short"

        bam_descriptor:
        "a one-word description of the purpose of the BAM (e.g. 'a_particular_seq_center'; used for saving the reads that don't have any MM/ML tags; doesn't need to be single-file specific)"

        run_nanoplot:
        "if true, will kick off Nanoplot to collect metrics on the uBAM; this isn't necessary if you also run the alignment version to process the uBAM as Nanoplot is run there automatically and produces a superset of metrics"

        convert_to_fq:
        "if true, input HiFi uBAM will be converted to FASTQ"
        hifi_fq:
        "available only if convert_to_fq is true."

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

        methyl_tag_simple_stats:
        "Simple stats on the reads with & without SAM methylation tags (MM/ML)."

        uBAM_metrics_files:
        "A map where keys are summary-names and values are paths to files generated from the various QC/metrics tasks"
    }

    input {
        File uBAM
        String readgroup_id
        File? input_hifi_fq

        String bam_descriptor
        Int short_reads_threshold

        Boolean run_nanoplot = false
        Boolean convert_to_fq = false
        Boolean i_donot_want_fq = false

        String disk_type = "SSD"

        String gcs_out_root_dir
    }

    output {
        String? movie = movie_name

        File? converted_hifi_fq = FinalizeFQ.gcs_path

        String last_processing_date = today.yyyy_mm_dd

        # todo merge these two together
        Map[String, Float] seqkit_stats = select_first([SeqKitOnBAM.stats, SeqKitOnFASTQ.stats])
        Map[String, Float]? nanoplot_u_summ = NanoPlotFromUBam.stats_map

        Map[String, String] read_len_summaries = DystPeaker.read_len_summaries
        Array[Int] read_len_peaks = DystPeaker.read_len_peaks
        Array[Int] read_len_deciles = DystPeaker.read_len_deciles

        # methylation call rate stats
        Map[String, String] methyl_tag_simple_stats = NoMissingBeans.methyl_tag_simple_stats

        # file-based outputs all packed into a finalization map
        Map[String, String] uBAM_metrics_files = files_out
    }

    ###################################################################################
    # prep work

    # arg validation
    if (defined(input_hifi_fq) == convert_to_fq) {
        if (!i_donot_want_fq) {
            call Utils.StopWorkflow as MutExProvided { input:
                reason = "You most likely want to either convert to FASTQ here, or have a pre-existing FASTQ, but not (miss) both."
            }
        }
    }

    # where to store final results
    String workflow_name = "ProcessOnInstrumentDemuxedChunkRefFree"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name

    String outdir_ref_free = outdir + '/RefFree'  # ideally, this isn't necessary, but it's a relic we'd not touch right now (2023-12-23)
    String bc_specific_out = outdir_ref_free + '/' + readgroup_id
    String bc_specific_metrics_out = bc_specific_out + "/metrics"

    ###################################################################################
    # convert to FASTQ if there isn't already one
    if (!defined(input_hifi_fq) && !i_donot_want_fq) {
        call BU.GetReadGroupInfo as RG { input: bam = uBAM, keys = ['PU']}
        String movie_name = RG.read_group_info['PU']

        ###################################################################################
        call BU.BamToFastq { input: bam = uBAM, prefix = "does_not_matter"}
        call FF.FinalizeToFile as FinalizeFQ { input:
            outdir = bc_specific_out,
            file = BamToFastq.reads_fq,
            name = readgroup_id + '.hifi.fq.gz'
        }
    }

    ###################################################################################
    # more QCs and metrics

    ##########
    if (run_nanoplot) {
        call QC0.NanoPlotFromUBam { input: uBAM = uBAM }
        FinalizationManifestLine a = object
                                    {files_to_save: flatten([[NanoPlotFromUBam.stats], NanoPlotFromUBam.plots]),
                                    is_singleton_file: false,
                                    destination: bc_specific_metrics_out + "/nanoplot",
                                    output_attribute_name: "nanoplot"}
    }
    ##########
    call QC1.CountTheBeans as NoMissingBeans { input:
        bam=uBAM,
        bam_descriptor=bam_descriptor,
        gcs_out_root_dir=bc_specific_metrics_out,
        use_local_ssd=disk_type=='LOCAL'
    }
    Map[String, String] methyl_out = {"missing_methyl_tag_reads":
                                      NoMissingBeans.methyl_tag_simple_stats['files_holding_reads_without_tags']}
    ##########
    call QC2.DystPeaker { input:
        input_file=uBAM, input_is_bam=true,
        id=readgroup_id, short_reads_threshold=short_reads_threshold,
        gcs_out_root_dir=bc_specific_metrics_out
    }
    Map[String, String] read_len_out = {"read_len_hist": DystPeaker.read_len_hist,
                                        "raw_rl_file": DystPeaker.read_len_summaries['raw_rl_file']}

    ##########
    # seqkit stats, detail depends on if FASTQ is available/desired or not
    if (i_donot_want_fq) {
        call QC3.FASTQstats as SeqKitOnBAM { input: reads = uBAM, file_type = "BAM" }
    }
    if (!i_donot_want_fq) {
        File use_this_hifi_fq = select_first([input_hifi_fq, BamToFastq.reads_fq])
        call QC3.FASTQstats as SeqKitOnFASTQ { input: reads=use_this_hifi_fq, file_type='FASTQ' }
    }

    ###################################################################################
    if (run_nanoplot) {
        call SAVE.SaveFilestoDestination as SaveRest { input:
            instructions = select_all([a]),
            already_finalized = select_all([methyl_out,
                                            read_len_out])
        }
    }
    if (!run_nanoplot) {
        call GU.MergeMaps { input:
            one = methyl_out,
            two = read_len_out,
        }
    }
    Map[String, String] files_out = select_first([SaveRest.result, MergeMaps.merged])

    call GU.GetTodayDate as today {}
}
