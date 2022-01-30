version 1.0

##########################################################################################
## A workflow that performs CCS correction on PacBio HiFi reads from a single flow cell.
## The workflow shards the subreads into clusters and performs CCS in parallel on each cluster.
## Ultimately, all the corrected reads (and uncorrected) are gathered into a single BAM.
## Various metrics are produced along the way.
##########################################################################################

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Longbow.wdl" as Longbow
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/NanoPlot.wdl" as NP
import "tasks/Finalize.wdl" as FF

workflow PBFlowcell {
    input {
        File bam
        File pbi
        File ref_map_file

        String SM
        String LB

        Boolean drop_per_base_N_pulse_tags = true

        Int? num_shards
        String experiment_type
        String dir_prefix

        String gcs_out_root_dir
    }

    parameter_meta {
        bam:                "GCS path to raw subread bam"
        pbi:                "GCS path to pbi index for raw subread bam"
        ref_map_file:       "table indicating reference sequence and auxillary file locations"

        SM:                 "the value to place in the BAM read group's SM field"
        LB:                 "the value to place in the BAM read group's LB (library) field"

        num_shards:         "number of shards into which fastq files should be batched [optional]"
        experiment_type:    "type of experiment run (CLR, CCS, ISOSEQ, MASSEQ)"
        dir_prefix:         "directory prefix for output files"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    DataTypeParameters clr    = { 'num_shards': select_first([num_shards, 100]), 'map_preset': 'SUBREAD' }
    DataTypeParameters ccs    = { 'num_shards': select_first([num_shards, 100]), 'map_preset': 'CCS'     }
    DataTypeParameters isoseq = { 'num_shards': select_first([num_shards,  50]), 'map_preset': 'ISOSEQ'  }
    DataTypeParameters masseq = { 'num_shards': select_first([num_shards, 300]), 'map_preset': 'ISOSEQ'  }
    Map[String, DataTypeParameters] data_presets = { 'CLR': clr, 'CCS': ccs, 'ISOSEQ': isoseq, 'MASSEQ': masseq }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBFlowcell/~{dir_prefix}"

    call PB.GetRunInfo { input: bam = bam, SM = SM }
    String PU = GetRunInfo.run_info['PU']

    # break one raw BAM into fixed number of shards
    call PB.ShardLongReads {
        input:
            unaligned_bam = bam,
            unaligned_pbi = pbi,
            num_shards = data_presets[experiment_type].num_shards,
            drop_per_base_N_pulse_tags = drop_per_base_N_pulse_tags
    }

    # for MAS-seq data, automatically detect the array model to use
    if (experiment_type == "MASSEQ") {
        call Longbow.Peek { input: bam = ShardLongReads.unmapped_shards[0], n = 100 }
    }

    # then perform correction and alignment on each of the shard
    scatter (unmapped_shard in ShardLongReads.unmapped_shards) {
        if (experiment_type != "CLR") {
            if (!GetRunInfo.is_corrected) { call PB.CCS { input: subreads = unmapped_shard } }

            if (experiment_type != "MASSEQ") {
                call PB.ExtractHifiReads {
                    input:
                        bam = select_first([CCS.consensus, unmapped_shard]),
                        sample_name = SM,
                        library     = LB
                }
            }

            if (experiment_type == 'MASSEQ') {
                call Longbow.Process { input: bam = select_first([unmapped_shard, CCS.consensus]), model = Peek.model }
            }
        }

        File unaligned_bam = select_first([Process.extracted_bam, ExtractHifiReads.hifi_bam, CCS.consensus, unmapped_shard])

        call PB.Align as AlignReads {
            input:
                bam         = unaligned_bam,
                ref_fasta   = ref_map['fasta'],
                sample_name = SM,
                library     = LB,
                map_preset  = data_presets[experiment_type].map_preset,
                drop_per_base_N_pulse_tags = drop_per_base_N_pulse_tags
        }

        call Utils.BamToFastq { input: bam = unaligned_bam, prefix = basename(unaligned_bam, ".bam") }
    }

    call Utils.MergeFastqs as MergeAllFastqs { input: fastqs = BamToFastq.reads_fq }

    # merge corrected, unaligned reads
    String cdir = outdir + "/reads/ccs/unaligned"
    if (experiment_type != "CLR" && experiment_type != "MASSEQ") {
        call Utils.MergeBams as MergeCCSUnalignedReads { input: bams = select_all(ExtractHifiReads.hifi_bam), prefix = "~{PU}.reads" }
        call PB.PBIndex as IndexCCSUnalignedReads { input: bam = MergeCCSUnalignedReads.merged_bam }

        call FF.FinalizeToFile as FinalizeCCSUnalignedBam { input: outdir = cdir, file = MergeCCSUnalignedReads.merged_bam }
        call FF.FinalizeToFile as FinalizeCCSUnalignedPbi {
            input:
                outdir = cdir,
                file = IndexCCSUnalignedReads.pbi,
                name = basename(MergeCCSUnalignedReads.merged_bam) + ".pbi"
        }
    }

    if (experiment_type != "CLR" && !GetRunInfo.is_corrected) {
        call PB.MergeCCSReports as MergeCCSReports { input: reports = select_all(CCS.report), prefix = PU }
        call PB.SummarizeCCSReport { input: report = MergeCCSReports.report }

        call FF.FinalizeToFile as FinalizeCCSReport { input: outdir = cdir, file = MergeCCSReports.report }
    }

    # merge the corrected per-shard BAM/report into one, corresponding to one raw input BAM
    call Utils.MergeBams as MergeAlignedReads { input: bams = AlignReads.aligned_bam, prefix = PU }
    call PB.PBIndex as IndexAlignedReads { input: bam = MergeAlignedReads.merged_bam }

    call AM.AlignedMetrics as PerFlowcellMetrics {
        input:
            aligned_bam    = MergeAlignedReads.merged_bam,
            aligned_bai    = MergeAlignedReads.merged_bai,
            ref_fasta      = ref_map['fasta'],
            ref_dict       = ref_map['dict'],
            gcs_output_dir = outdir + "/metrics"
    }

    call PB.SummarizePBI as SummarizeSubreadsPBI   { input: pbi = pbi, runtime_attr_override = { 'mem_gb': 72 } }
    call PB.SummarizePBI as SummarizeAlignedPBI    { input: pbi = IndexAlignedReads.pbi }
    call PB.SummarizePBI as SummarizeAlignedQ5PBI  { input: pbi = IndexAlignedReads.pbi, qual_threshold = 5 }
    call PB.SummarizePBI as SummarizeAlignedQ7PBI  { input: pbi = IndexAlignedReads.pbi, qual_threshold = 7 }
    call PB.SummarizePBI as SummarizeAlignedQ10PBI { input: pbi = IndexAlignedReads.pbi, qual_threshold = 10 }
    call PB.SummarizePBI as SummarizeAlignedQ12PBI { input: pbi = IndexAlignedReads.pbi, qual_threshold = 12 }
    call PB.SummarizePBI as SummarizeAlignedQ15PBI { input: pbi = IndexAlignedReads.pbi, qual_threshold = 15 }

    call NP.NanoPlotFromBam { input: bam = MergeAlignedReads.merged_bam, bai = MergeAlignedReads.merged_bai }
    call Utils.ComputeGenomeLength { input: fasta = ref_map['fasta'] }

    # Finalize data
    String dir = outdir + "/" + if (experiment_type != "CLR") then "reads/ccs/aligned" else "reads/subreads/aligned"

    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = dir, file = MergeAlignedReads.merged_bam }
    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = dir, file = MergeAlignedReads.merged_bai }
    call FF.FinalizeToFile as FinalizeAlignedPbi {
        input:
            outdir = dir,
            file = IndexAlignedReads.pbi,
            name = basename(MergeAlignedReads.merged_bam) + ".pbi"
    }

    String fqdir = outdir + "/" + if (experiment_type != "CLR") then "reads/ccs/unaligned" else "reads/subreads/unaligned"
    call FF.FinalizeToFile as FinalizeFastq {
        input:
            outdir = fqdir,
            file = MergeAllFastqs.merged_fastq,
            name = basename(MergeAlignedReads.merged_bam, ".bam") + ".fq.gz"
    }

    output {
        # Flowcell stats
        File? ccs_report = FinalizeCCSReport.gcs_path
        Float? ccs_zmws_input = SummarizeCCSReport.zmws_input
        Float? ccs_zmws_pass_filters = SummarizeCCSReport.zmws_pass_filters
        Float? ccs_zmws_fail_filters = SummarizeCCSReport.zmws_fail_filters
        Float? ccs_zmws_shortcut_filters = SummarizeCCSReport.zmws_shortcut_filters
        Float? ccs_zmws_pass_filters_pct = SummarizeCCSReport.zmws_pass_filters_pct
        Float? ccs_zmws_fail_filters_pct = SummarizeCCSReport.zmws_fail_filters_pct
        Float? ccs_zmws_shortcut_filters_pct = SummarizeCCSReport.zmws_shortcut_filters_pct

        Float polymerase_read_length_mean = SummarizeSubreadsPBI.results['polymerase_mean']
        Float polymerase_read_length_N50 = SummarizeSubreadsPBI.results['polymerase_n50']

        Float subread_read_length_mean = SummarizeSubreadsPBI.results['subread_mean']
        Float subread_read_length_N50 = SummarizeSubreadsPBI.results['subread_n50']

        # Unaligned reads
        File? fq = FinalizeFastq.gcs_path

        # Unaligned BAM file
        File? ccs_bam = FinalizeCCSUnalignedBam.gcs_path
        File? ccs_pbi = FinalizeCCSUnalignedPbi.gcs_path

        # Aligned BAM file
        File aligned_bam = FinalizeAlignedBam.gcs_path
        File aligned_bai = FinalizeAlignedBai.gcs_path
        File aligned_pbi = FinalizeAlignedPbi.gcs_path

        # Unaligned read stats
        Float num_reads = SummarizeSubreadsPBI.results['reads']
        Float num_bases = SummarizeSubreadsPBI.results['bases']
        Float raw_est_fold_cov = SummarizeSubreadsPBI.results['bases']/ComputeGenomeLength.length

        Float read_length_mean = SummarizeSubreadsPBI.results['subread_mean']
        Float read_length_median = SummarizeSubreadsPBI.results['subread_median']
        Float read_length_stdev = SummarizeSubreadsPBI.results['subread_stdev']
        Float read_length_N50 = SummarizeSubreadsPBI.results['subread_n50']

        Float read_qual_mean = SummarizeSubreadsPBI.results['mean_qual']
        Float read_qual_median = SummarizeSubreadsPBI.results['median_qual']

        Float num_reads_Q5 = SummarizeAlignedQ5PBI.results['reads']
        Float num_reads_Q7 = SummarizeAlignedQ7PBI.results['reads']
        Float num_reads_Q10 = SummarizeAlignedQ10PBI.results['reads']
        Float num_reads_Q12 = SummarizeAlignedQ12PBI.results['reads']
        Float num_reads_Q15 = SummarizeAlignedQ15PBI.results['reads']

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
