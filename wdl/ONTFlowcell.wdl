version 1.0

import "tasks/ONTUtils.wdl" as ONT
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/NanoPlot.wdl" as NP
import "tasks/Finalize.wdl" as FF

workflow ONTFlowcell {
    input {
        File final_summary
        File sequencing_summary
        File ref_map_file

        String SM
        String ID

        Int num_shards = 50
        String experiment_type

        String gcs_out_root_dir
    }

    parameter_meta {
        final_summary:      "GCS path to '*final_summary*.txt*' file for basecalled fastq files"
        sequencing_summary: "GCS path to '*sequencing_summary*.txt*' file for basecalled fastq files"
        ref_map_file:       "table indicating reference sequence and auxillary file locations"

        SM:                 "the value to place in the BAM read group's SM field"
        ID:                 "the value to place in the BAM read group's ID field"

        num_shards:         "[default-valued] number of shards into which fastq files should be batched"
        experiment_type:    "[default-valued] type of experiment run (DNA, RNA, R2C2)"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)
    Map[String, String] map_presets = {
        'DNA':  'map-ont',
        'RNA':  'splice',
        'R2C2': 'splice'
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTFlowcell/~{ID}"

    call ONT.GetRunInfo { input: final_summary = final_summary }

    call ONT.ListFiles as ListFastqs { input: sequencing_summary = sequencing_summary, suffix = "fastq" }

    String PL  = "ONT"
    String PU  = GetRunInfo.run_info["instrument"]
    String DT  = GetRunInfo.run_info["started"]
    String RG = "@RG\\tID:~{ID}\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

    call ONT.PartitionManifest as PartitionFastqManifest { input: manifest = ListFastqs.manifest, N = num_shards }

    scatter (manifest_chunk in PartitionFastqManifest.manifest_chunks) {
        call AR.Minimap2 as AlignReads {
            input:
                reads      = read_lines(manifest_chunk),
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

    call NP.NanoPlotFromSummary { input: summary_files = [sequencing_summary] }
    call NP.NanoPlotFromBam { input: bam = MergeAlignedReads.merged_bam, bai = MergeAlignedReads.merged_bai }
    call Utils.ComputeGenomeLength { input: fasta = ref_map['fasta'] }

    # Finalize data
    String dir = outdir + "/alignments"

    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = dir, file = MergeAlignedReads.merged_bam }
    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = dir, file = MergeAlignedReads.merged_bai }

    output {
        # Flowcell stats
        Float active_channels = NanoPlotFromSummary.stats_map['active_channels']

        # Aligned BAM file
        File aligned_bam = FinalizeAlignedBam.gcs_path
        File aligned_bai = FinalizeAlignedBai.gcs_path

        # Unaligned read stats
        Float num_reads = NanoPlotFromSummary.stats_map['number_of_reads']
        Float num_bases = NanoPlotFromSummary.stats_map['number_of_bases']
        Float raw_est_fold_cov = NanoPlotFromSummary.stats_map['number_of_bases']/ComputeGenomeLength.length

        Float read_length_mean = NanoPlotFromSummary.stats_map['mean_read_length']
        Float read_length_median = NanoPlotFromSummary.stats_map['median_read_length']
        Float read_length_stdev = NanoPlotFromSummary.stats_map['read_length_stdev']
        Float read_length_N50 = NanoPlotFromSummary.stats_map['n50']

        Float read_qual_mean = NanoPlotFromSummary.stats_map['mean_qual']
        Float read_qual_median = NanoPlotFromSummary.stats_map['median_qual']

        Float num_reads_Q5 = NanoPlotFromSummary.stats_map['Reads_Q5']
        Float num_reads_Q7 = NanoPlotFromSummary.stats_map['Reads_Q7']
        Float num_reads_Q10 = NanoPlotFromSummary.stats_map['Reads_Q10']
        Float num_reads_Q12 = NanoPlotFromSummary.stats_map['Reads_Q12']
        Float num_reads_Q15 = NanoPlotFromSummary.stats_map['Reads_Q15']

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