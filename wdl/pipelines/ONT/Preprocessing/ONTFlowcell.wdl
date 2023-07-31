version 1.0

import "../../../tasks/Utility/ONTUtils.wdl" as ONT
import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Alignment/AlignReads.wdl" as AR
import "../../../tasks/QC/AlignedMetrics.wdl" as AM
import "../../../tasks/Visualization/NanoPlot.wdl" as NP
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ONTFlowcell {

    meta {
        description: "Align ONT reads to a reference genome"
    }
    parameter_meta {
        final_summary:      "GCS path to '*final_summary*.txt*' file for basecalled fastq files"
        sequencing_summary: "GCS path to '*sequencing_summary*.txt*' file for basecalled fastq files"
        fastq_dir:          "GCS path to fastq directory"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"

        SM:                 "the value to place in the BAM read group's SM field"
        ID:                 "the value to place in the BAM read group's ID field"

        num_shards:         "[default-valued] number of shards into which fastq files should be batched"
        experiment_type:    "[default-valued] type of experiment run (DNA, RNA, R2C2)"
        dir_prefix:         "directory prefix for output files"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        File? final_summary
        File? sequencing_summary
        String? fastq_dir

        File ref_map_file

        String SM
        String ID

        Int num_shards = 300
        String experiment_type
        String dir_prefix

        String gcs_out_root_dir
    }

    Map[String, String] ref_map = read_map(ref_map_file)
    Map[String, String] map_presets = {
        'DNA':  'map-ont',
        'RNA':  'splice',
        'R2C2': 'splice'
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTFlowcell/~{dir_prefix}"

    if (defined(final_summary)) {
        call ONT.GetRunInfo { input: final_summary = select_first([final_summary]) }
    }
    Map[String, String] runinfo = select_first([GetRunInfo.run_info, { "instrument": "unknown", "started": "2021-01-01T12:00:00.000000-05:00" }])
    String PU = runinfo['instrument']
    String DT = runinfo['started']

    if (defined(sequencing_summary)) {
        call ONT.ListFiles as ListFastqs {
            input:
                sequencing_summary = select_first([sequencing_summary]),
                suffix = "fastq"
        }

        call NP.NanoPlotFromSummary { input: summary_files = [ select_first([sequencing_summary]) ] }
    }

    if (defined(fastq_dir)) {
        call Utils.ListFilesOfType {
            input:
                gcs_dir = select_first([fastq_dir]),
                suffixes = [ '.fastq', '.fastq.gz', '.fq', '.fq.gz' ],
                recurse = true
        }

        call NP.NanoPlotFromRichFastqs { input: fastqs = ListFilesOfType.files }
    }

    Map[String, Float] nanoplot_map = select_first([NanoPlotFromRichFastqs.stats_map, NanoPlotFromSummary.stats_map ])

    File manifest = select_first([ListFilesOfType.manifest, ListFastqs.manifest])

    String PL  = "ONT"
    String RG = "@RG\\tID:~{ID}\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

    call ONT.PartitionManifest as PartitionFastqManifest { input: manifest = manifest, N = num_shards }

    scatter (manifest_chunk in PartitionFastqManifest.manifest_chunks) {
        Array[File] reads_arr = read_lines(manifest_chunk)
        scatter (dummy in reads_arr) {
            String read_basename = basename(dummy)
        }
        Array[String] reads_file_basenames = read_basename
        call AR.Minimap2 as AlignReads {
            input:
                reads      = reads_arr,
                reads_file_basenames = reads_file_basenames,
                ref_fasta  = ref_map['fasta'],
                RG         = RG,
                map_preset = map_presets[experiment_type]
        }
    }

    call Utils.MergeBams as MergeAlignedReads { input: bams = AlignReads.aligned_bam, prefix = ID }

    call AM.AlignedMetrics as PerFlowcellMetrics {
        input:
            aligned_bam    = MergeAlignedReads.merged_bam,
            aligned_bai    = MergeAlignedReads.merged_bai,
            ref_fasta      = ref_map['fasta'],
            ref_dict       = ref_map['dict'],
            gcs_output_dir = outdir + "/metrics"
    }

    call NP.NanoPlotFromBam { input: bam = MergeAlignedReads.merged_bam, bai = MergeAlignedReads.merged_bai }
    call Utils.ComputeGenomeLength { input: fasta = ref_map['fasta'] }

    # Finalize data
    String dir = outdir + "/alignments"

    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = dir, file = MergeAlignedReads.merged_bam }
    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = dir, file = MergeAlignedReads.merged_bai }

    output {
        # Flowcell stats
        Float active_channels = nanoplot_map['active_channels']

        # Aligned BAM file
        File aligned_bam = FinalizeAlignedBam.gcs_path
        File aligned_bai = FinalizeAlignedBai.gcs_path

        # Unaligned read stats
        Float num_reads = nanoplot_map['number_of_reads']
        Float num_bases = nanoplot_map['number_of_bases']
        Float raw_est_fold_cov = nanoplot_map['number_of_bases']/ComputeGenomeLength.length

        Float read_length_mean = nanoplot_map['mean_read_length']
        Float read_length_median = nanoplot_map['median_read_length']
        Float read_length_stdev = nanoplot_map['read_length_stdev']
        Float read_length_N50 = nanoplot_map['n50']

        Float read_qual_mean = nanoplot_map['mean_qual']
        Float read_qual_median = nanoplot_map['median_qual']

        Float num_reads_Q5 = nanoplot_map['Reads_Q5']
        Float num_reads_Q7 = nanoplot_map['Reads_Q7']
        Float num_reads_Q10 = nanoplot_map['Reads_Q10']
        Float num_reads_Q12 = nanoplot_map['Reads_Q12']
        Float num_reads_Q15 = nanoplot_map['Reads_Q15']

        # Aligned read stats
        Float aligned_num_reads = NanoPlotFromBam.stats_map['number_of_reads']
        Float aligned_num_bases = NanoPlotFromBam.stats_map['number_of_bases_aligned']
        Float aligned_frac_bases = NanoPlotFromBam.stats_map['fraction_bases_aligned']
        Float aligned_est_fold_cov = NanoPlotFromBam.stats_map['number_of_bases_aligned']/ComputeGenomeLength.length

        Float aligned_read_length_mean = NanoPlotFromBam.stats_map['mean_read_length']
        Float aligned_read_length_median = NanoPlotFromBam.stats_map['median_read_length']
        Float aligned_read_length_stdev = NanoPlotFromBam.stats_map['read_length_stdev']
        Float aligned_read_length_N50 = NanoPlotFromBam.stats_map['n50']

        Float average_identity = NanoPlotFromBam.stats_map['average_identity']
        Float median_identity = NanoPlotFromBam.stats_map['median_identity']
    }
}