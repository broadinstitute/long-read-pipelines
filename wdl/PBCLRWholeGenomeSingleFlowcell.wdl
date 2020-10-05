version 1.0

##########################################################################################
## A workflow that performs CCS correction and variant calling on PacBio HiFi reads from a
## single flow cell. The workflow shards the subreads into clusters and performs CCS in
## parallel on each cluster.  Error-corrected reads are then variant-called.  A number of
## metrics and figures are produced along the way.
##########################################################################################

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.36/wdl/tasks/PBUtils.wdl" as PB
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.36/wdl/tasks/ShardUtils.wdl" as SU
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.36/wdl/tasks/Utils.wdl" as Utils
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.36/wdl/tasks/AlignReads.wdl" as AR
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.36/wdl/tasks/AlignedMetrics.wdl" as AM
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.36/wdl/tasks/CallSVs.wdl" as SV
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.36/wdl/tasks/Figures.wdl" as FIG
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.36/wdl/tasks/Finalize.wdl" as FF
import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.36/wdl/tasks/CallSmallVariants.wdl" as SMV

workflow PBCLRWholeGenomeSingleFlowcell {
    input {
        String raw_reads_gcs_bucket
        String? sample_name

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File tandem_repeat_bed
        File ref_flat
        File dbsnp_vcf
        File dbsnp_tbi

        String mt_chr_name
        File metrics_locus

        Int? num_reads_per_split = 2000000

        String gcs_out_root_dir
    }

    parameter_meta {
        raw_reads_gcs_bucket: "GCS bucket holding subreads BAMs (and other related files) holding the sequences to be CCS-ed"
        sample_name:          "[optional] name of sample this FC is sequencing"

        ref_fasta:            "Reference fasta file"
        ref_fasta_fai:        "Index (.fai) for the reference fasta file"
        ref_dict:             "Sequence dictionary (.dict) for the reference fasta file"

        tandem_repeat_bed:    "BED file specifying the location of tandem repeats in the reference"
        ref_flat:             "Gene predictions in refFlat format (https://genome.ucsc.edu/goldenpath/gbdDescriptions.html)"
        dbsnp_vcf:            "dbSNP vcf"
        dbsnp_tbi:            "Index (.tbi) for dbSNP vcf"

        mt_chr_name:          "Contig name for the mitochondrial sequence in the reference"
        metrics_locus:        "Loci over which some summary metrics should be computed"

        num_reads_per_split:  "[default-valued] number of subreads each sharded BAM contains (tune for performance)"

        gcs_out_root_dir :    "GCS bucket to store the corrected/uncorrected reads and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = raw_reads_gcs_bucket }

    # double scatter: one FC may generate multiple raw BAMs, we perform another layer scatter on each of these BAMs
    scatter (subread_bam in FindBams.subread_bams) {
        call PB.GetRunInfo { input: subread_bam = subread_bam }

        String SM  = select_first([sample_name, GetRunInfo.run_info["SM"]])
        String PL  = "PACBIO"
        String PU  = GetRunInfo.run_info["PU"]
        String DT  = GetRunInfo.run_info["DT"]
        String ID  = PU
        String DS  = GetRunInfo.run_info["DS"]
        String DIR = SM + "." + ID
        String RG = "@RG\\tID:~{ID}\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        # shard one raw BAM into fixed chunk size (num_reads_per_split)
        call Utils.ShardLongReads { input: unmapped_files = [ subread_bam ], num_reads_per_split = num_reads_per_split }

        # then perform correction and alignment on each of the shard
        scatter (subreads in ShardLongReads.unmapped_shards) {
            call AR.Minimap2 as AlignChunk {
                input:
                    reads      = [ subreads ],
                    ref_fasta  = ref_fasta,
                    RG         = RG,
                    map_preset = "map-pb"
            }
        }

        # merge the aligned per-shard BAM into one, corresponding to one raw input BAM
        call Utils.MergeBams as MergeChunks { input: bams = AlignChunk.aligned_bam }

        # compute alignment metrics
        call AM.AlignedMetrics as PerFlowcellSubRunMetrics {
            input:
                aligned_bam    = MergeChunks.merged_bam,
                aligned_bai    = MergeChunks.merged_bai,
                ref_fasta      = ref_fasta,
                ref_dict       = ref_dict,
                ref_flat       = ref_flat,
                dbsnp_vcf      = dbsnp_vcf,
                dbsnp_tbi      = dbsnp_tbi,
                metrics_locus  = metrics_locus,
                per            = "flowcell",
                type           = "subrun",
                label          = ID,
                gcs_output_dir = outdir + "/" + DIR
        }

#        call FIG.Figures as PerFlowcellSubRunFigures {
#            input:
#                summary_files  = [ summary_file ],
#
#                per            = "flowcell",
#                type           = "subrun",
#                label          = SID,
#
#                gcs_output_dir = outdir + "/" + DIR
#        }
    }

    # gather across (potential multiple) input raw BAMs
    if (length(FindBams.subread_bams) > 1) {
        call Utils.MergeBams as MergeRuns { input: bams = MergeChunks.merged_bam, prefix = "~{SM[0]}.~{ID[0]}" }
    }

    File merged_bam = select_first([ MergeRuns.merged_bam, MergeChunks.merged_bam[0] ])
    File merged_bai = select_first([ MergeRuns.merged_bai, MergeChunks.merged_bai[0] ])

    # compute alignment metrics
    call AM.AlignedMetrics as PerFlowcellRunMetrics {
        input:
            aligned_bam    = merged_bam,
            aligned_bai    = merged_bai,
            ref_fasta      = ref_fasta,
            ref_dict       = ref_dict,
            ref_flat       = ref_flat,
            dbsnp_vcf      = dbsnp_vcf,
            dbsnp_tbi      = dbsnp_tbi,
            metrics_locus  = metrics_locus,
            per            = "flowcell",
            type           = "run",
            label          = ID[0],
            gcs_output_dir = outdir + "/" + DIR[0]
    }

#    call FIG.Figures as PerFlowcellRunFigures {
#        input:
#            summary_files  = FindSequencingSummaryFiles.summary_files,
#
#            per            = "flowcell",
#            type           = "run",
#            label          = ID[0],
#
#            gcs_output_dir = outdir + "/" + DIR[0]
#    }

    # call SVs
    call SV.CallSVs as CallSVs {
        input:
            bam               = merged_bam,
            bai               = merged_bai,

            ref_fasta         = ref_fasta,
            ref_fasta_fai     = ref_fasta_fai,
            tandem_repeat_bed = tandem_repeat_bed
    }

    # call SNVs and small indels
    call SMV.CallSmallVariants as CallSmallVariants {
        input:
            bam               = merged_bam,
            bai               = merged_bai,

            ref_fasta         = ref_fasta,
            ref_fasta_fai     = ref_fasta_fai,
            ref_dict          = ref_dict
    }

    ##########
    # store the results into designated bucket
    ##########

    call FF.FinalizeToDir as FinalizeSVs {
        input:
            files = [ CallSVs.pbsv_vcf, CallSVs.sniffles_vcf, CallSVs.svim_vcf ],
            outdir = outdir + "/" + DIR[0] + "/variants"
    }

    call FF.FinalizeToDir as FinalizeSmallVariants {
        input:
            files = [ CallSmallVariants.longshot_vcf, CallSmallVariants.longshot_tbi ],
            outdir = outdir + "/" + DIR[0] + "/variants"
    }

    call FF.FinalizeToDir as FinalizeMergedRuns {
        input:
            files = [ merged_bam, merged_bai ],
            outdir = outdir + "/" + DIR[0] + "/alignments"
    }
}
