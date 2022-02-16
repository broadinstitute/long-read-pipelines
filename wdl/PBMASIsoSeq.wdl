version 1.0

######################################################################################
## A workflow that performs processing of MAS-ISO-seq data on a single sample from
## one or more flow cells. The workflow merges multiple samples into a single BAM
## prior to processing.
######################################################################################

import "tasks/Utils.wdl" as Utils
import "tasks/PBUtils.wdl" as PB
import "tasks/AlignReads.wdl" as AR
#import "tasks/NanoPlot.wdl" as NP
#import "tasks/SQANTI.wdl"
import "tasks/StringTie2.wdl" as ST2
import "tasks/Ten_X_Tool.wdl" as TENX
import "tasks/TranscriptAnalysis/UMI_Tools.wdl" as UMI_TOOLS
import "tasks/TranscriptAnalysis/Postprocessing_Tasks.wdl" as TX_POST
import "tasks/Finalize.wdl" as FF

workflow PBMASIsoSeq {
    input {
        Array[File] aligned_bams
        Array[File] aligned_bais
        String participant_name

        File ref_map_file
        File ref_gtf

        String gcs_out_root_dir
    }

    parameter_meta {
        aligned_bams:     "GCS path to aligned BAM files"
        aligned_bais:     "GCS path to aligned BAM file indices"
        participant_name: "name of the participant from whom these samples were obtained"

        ref_map_file:     "table indicating reference sequence and auxillary file locations"
        ref_gtf:          "GTF file to use for quantification"

        gcs_out_root_dir: "GCS bucket to store the corrected/uncorrected reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBMASIsoSeq/~{participant_name}"

    # gather across (potential multiple) input CCS BAMs
    if (length(aligned_bams) > 1) {
        call Utils.MergeBams as MergeAllReads { input: bams = aligned_bams, prefix = participant_name }
    }

    call Utils.SubsetBam {
        input:
            bam = select_first([MergeAllReads.merged_bam, aligned_bams[0]]),
            bai = select_first([MergeAllReads.merged_bai, aligned_bais[0]]),
            locus = ["chr21"]
    }

    File bam = SubsetBam.subset_bam
    File bai = SubsetBam.subset_bai

#    File bam = select_first([MergeAllReads.merged_bam, aligned_bams[0]])
#    File bai = select_first([MergeAllReads.merged_bai, aligned_bais[0]])

    call PB.PBIndex as IndexCCSAlignedReads { input: bam = bam }
    File pbi = IndexCCSAlignedReads.pbi

    call ST2.Quantify {
        input:
            aligned_bam = bam,
            aligned_bai = bai,
            gtf = ref_gtf,
            keep_retained_introns = false,
            prefix = participant_name + "_StringTie2_Quantify",
    }

    call ST2.ExtractTranscriptSequences {
        input:
            ref_fasta = ref_map['fasta'],
            ref_fasta_fai = ref_map['fai'],
            gtf = Quantify.st_gtf,
            prefix = participant_name + "_StringTie2_ExtractTranscriptSequences",
    }

    call ST2.CompareTranscriptomes {
        input:
            guide_gtf = ref_gtf,
            new_gtf = Quantify.st_gtf,
            prefix = participant_name + "_StringTie2_CompareTranscriptome",
    }

    # Set our transcriptome files:
    File transcriptome_reference_for_quant = ExtractTranscriptSequences.transcripts_fa
    File transcriptome_reference_index_for_quant = ExtractTranscriptSequences.transcripts_fai
    File transcriptome_reference_dict_for_quant = ExtractTranscriptSequences.transcripts_dict

    # break one raw BAM into fixed number of shards
#    call PB.ShardLongReads { input: unaligned_bam = bam, unaligned_pbi = pbi, num_shards = 50 }
    Array[String] default_filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']
    call Utils.MakeChrIntervalList { input: ref_dict = ref_map['dict'], filter = default_filter }

    # Now we have to align the array elements to the new transcriptome.
    scatter (c in MakeChrIntervalList.chrs) {
        String contig = c[0]

        call Utils.SubsetBam { input: bam = bam, bai = bai, locus = contig }

        # Align our array elements:
        call PB.Align as AlignToTranscriptome {
            input:
                bam = SubsetBam.subset_bam,
                ref_fasta = transcriptome_reference_for_quant,
                sample_name = participant_name,
                drop_per_base_N_pulse_tags = true,
                map_preset = 'SUBREAD',
        }

        # To properly count our transcripts we must throw away the non-primary and unaligned reads:
        call Utils.FilterReadsBySamFlags as RemoveUnmappedAndNonPrimaryReads {
            input:
                bam = AlignToTranscriptome.aligned_bam,
                sam_flags = "2308",
                prefix = participant_name + "_ArrayElements_Annotated_Aligned_PrimaryOnly",
                runtime_attr_override = object { cpu_cores: 4, preemptible_tries: 0 }
        }

        # Filter reads with no UMI tag:
        call Utils.FilterReadsWithTagValues as FilterReadsWithNoUMI {
            input:
                bam = RemoveUnmappedAndNonPrimaryReads.output_bam,
                tag = "ZU",
                value_to_remove = ".",
                prefix = participant_name + "_ArrayElements_Annotated_Aligned_PrimaryOnly_WithUMIs",
                runtime_attr_override = object { cpu_cores: 4, preemptible_tries: 0 }
        }

        # Copy the contig to a tag.
        # By this point in the pipeline, array elements are aligned to a transcriptome, so this tag will
        # actually indicate the transcript to which each array element aligns.
        call TENX.CopyContigNameToReadTag as CopyContigNameToReadTag {
            input:
                aligned_bam_file = FilterReadsWithNoUMI.output_bam,
                prefix = participant_name + "_ArrayElements_Annotated_Aligned_PrimaryOnly_WithUMIs"
        }
    }

    # Now we merge together our TX-ome aligned stuff:
    call Utils.MergeBams as MergeTranscriptomeAlignedExtractedArrayElements {
        input:
            bams = AlignToTranscriptome.aligned_bam,
            prefix = participant_name + "_array_elements_longbow_extracted_tx_aligned",
            runtime_attr_override = object { cpu_cores: 4 }
    }

    call Utils.MergeBams as MergePrimaryTranscriptomeAlignedArrayElements {
        input:
            bams = CopyContigNameToReadTag.output_bam,
            prefix = participant_name + "_array_elements_longbow_extracted_tx_aligned_primary_alignments",
            runtime_attr_override = object { cpu_cores: 4 }
    }

    ##########
    # Quantify Transcripts:
    ##########

    call UMI_TOOLS.Run_Group as UMIToolsGroup {
        input:
            aligned_transcriptome_reads = MergePrimaryTranscriptomeAlignedArrayElements.merged_bam,
            aligned_transcriptome_reads_index = MergePrimaryTranscriptomeAlignedArrayElements.merged_bai,
            do_per_cell = true,
            prefix = "~{participant_name}_umi_tools_group"
    }

    call TX_POST.CreateCountMatrixFromAnnotatedBam as t_57_CreateCountMatrixFromAnnotatedBam {
        input:
            annotated_transcriptome_bam = UMIToolsGroup.output_bam,
            prefix = "~{participant_name}_gene_tx_expression_count_matrix"
    }

    # Only create the anndata objects if we're looking at real genomic data:
    call TX_POST.CreateCountMatrixAnndataFromTsv as t_58_CreateCountMatrixAnndataFromTsv {
        input:
            count_matrix_tsv = t_57_CreateCountMatrixFromAnnotatedBam.count_matrix,
            genome_annotation_gtf_file = select_first([Quantify.st_gtf]),
            gencode_reference_gtf_file = ref_gtf,
            prefix = "~{participant_name}_gene_tx_expression_count_matrix"
    }

#    call SQANTI.QC {
#        input:
#            bam = MergeAligned.merged_bam,
#            ref_fasta = ref_map['fasta'],
#            ref_gtf = ref_gtf,
#            prefix = prefix
#    }

    # Finalize data
    String adir = outdir + "/alignments"

    call FF.FinalizeToFile as FinalizeBam { input: outdir = adir, file = bam, name = "~{participant_name}.bam" }
    call FF.FinalizeToFile as FinalizeBai { input: outdir = adir, file = bai, name = "~{participant_name}.bam.bai" }

#    String updir = outdir + "/stats/unfiltered/png"
#    String usdir = outdir + "/stats/unfiltered/svg"
#    call FF.FinalizeToDir as FinalizeStatsUnfilteredPng { input: outdir = updir, files = StatsUnfiltered.pngs }
#    call FF.FinalizeToDir as FinalizeStatsUnfilteredSvg { input: outdir = usdir, files = StatsUnfiltered.svgs }
#
#    String fpdir = outdir + "/stats/filtered/png"
#    String fsdir = outdir + "/stats/filtered/svg"
#    call FF.FinalizeToDir as FinalizeStatsFilteredPng { input: outdir = fpdir, files = StatsFiltered.pngs }
#    call FF.FinalizeToDir as FinalizeStatsFilteredSvg { input: outdir = fsdir, files = StatsFiltered.svgs }
#
#    call FF.FinalizeToFile as FinalizeNPRichFqStats { input: outdir = outdir + "/stats/nanoplot/fastq", file = NanoPlotFromRichFastqs.stats }
#    call FF.FinalizeToDir as FinalizeNPRichFqPlots { input: outdir = outdir + "/stats/nanoplot/fastq", files = NanoPlotFromRichFastqs.plots }
#
#    call FF.FinalizeToFile as FinalizeNPBamStats { input: outdir = outdir + "/stats/nanoplot/bam", file = NanoPlotFromBam.stats }
#    call FF.FinalizeToDir as FinalizeNPBamPlots { input: outdir = outdir + "/stats/nanoplot/bam", files = NanoPlotFromBam.plots }
#
#    String qdir = outdir + "/stats/transcripts/"
#    call FF.FinalizeToFile as FinalizeClassifications { input: outdir = qdir, file = QC.classification }
#    call FF.FinalizeToFile as FinalizeJunctions { input: outdir = qdir, file = QC.junctions }
#    call FF.FinalizeToFile as FinalizeReport { input: outdir = qdir, file = QC.report }

    output {
        File merged_bam = FinalizeBam.gcs_path
        File merged_bai = FinalizeBai.gcs_path

        # String stats_unfiltered_pngs = FinalizeStatsUnfilteredPng.gcs_dir
        # String stats_unfiltered_svgs = FinalizeStatsUnfilteredSvg.gcs_dir
        # String stats_filtered_pngs = FinalizeStatsFilteredPng.gcs_dir
        # String stats_filtered_svgs = FinalizeStatsFilteredSvg.gcs_dir

        # File nanoplot_fq_stats = FinalizeNPRichFqStats.gcs_path
        # File nanoplot_fq_dir = FinalizeNPRichFqPlots.gcs_dir

        # File nanoplot_bam_stats = FinalizeNPBamStats.gcs_path
        # File nanoplot_bam_dir = FinalizeNPBamPlots.gcs_dir

        # File classifications = FinalizeClassifications.gcs_path
        # File junctions = FinalizeJunctions.gcs_path
        # File report = FinalizeReport.gcs_path
    }
}

task Process {
    input {
        Array[File] fastq_files
        String model_name
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 * ceil(size(fastq_files, "GB"))

    command <<<
        set -euxo pipefail

        DIR=$(dirname ~{fastq_files[0]})

        longbow convert $DIR | \
            longbow annotate -m ~{model_name} | \
                tee ~{prefix}.annotated_unfiltered.bam | \
            longbow filter | \
                tee ~{prefix}.annotated_filtered.bam | \
            longbow segment | \
            longbow extract -o ~{prefix}.extracted.bam
    >>>

    output {
        File annotated_unfiltered = "~{prefix}.annotated_unfiltered.bam"
        File annotated_filtered = "~{prefix}.annotated_filtered.bam"
        File extracted = "~{prefix}.extracted.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             4,
        disk_gb:            1,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.4.7-kvg6"
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

task Stats {
    input {
        File bam
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2 * ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        longbow stats -o ~{prefix}.stats ~{bam}
    >>>

    output {
        Array[File] pngs = glob("*.png")
        Array[File] svgs = glob("*.svg")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longbow:0.4.7-kvg6"
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
