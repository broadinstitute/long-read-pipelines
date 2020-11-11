version 1.0

##########################################################################################
## A workflow that performs CCS correction and variant calling on PacBio HiFi reads from a
## single flow cell. The workflow shards the subreads into clusters and performs CCS in
## parallel on each cluster.  Error-corrected reads are then variant-called.  A number of
## metrics and figures are produced along the way.
##########################################################################################

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/CallSVs.wdl" as SV
import "tasks/Figures.wdl" as FIG
import "tasks/Finalize.wdl" as FF
import "tasks/CallSmallVariants.wdl" as SMV

workflow PBCLRWholeGenome {
    input {
        Array[File] bams
        File ref_map_file

        String participant_name
        Int num_shards = 300

        String? gcs_out_root_dir
    }

    parameter_meta {
        bams:                      "GCS path to raw subreads or CCS data"
        ref_map_file:              "table indicating reference sequence and auxillary file locations"

        participant_name:          "name of the participant from whom these samples were obtained"
        num_shards:                "[default-valued] number of sharded BAMs to create (tune for performance)"

        gcs_out_root_dir:          "[optional] GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    call Utils.GetDefaultDir { input: workflow_name = "PBCLRWholeGenome" }
    String outdir = sub(select_first([gcs_out_root_dir, GetDefaultDir.path]), "/$", "") + "/" + participant_name

    # scatter over all sample BAMs
    scatter (bam in bams) {
        File pbi = sub(bam, ".bam$", ".bam.pbi")

        call PB.GetRunInfo { input: bam = bam }
        String ID = GetRunInfo.run_info["PU"]

        # break one raw BAM into fixed number of shards
        call PB.ShardLongReads { input: unaligned_bam = bam, unaligned_pbi = pbi, num_shards = num_shards }

        # then perform correction and alignment on each of the shard
        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.Align as AlignUncorrected {
                input:
                    bam         = subreads,
                    ref_fasta   = ref_map['fasta'],
                    sample_name = participant_name,
                    map_preset  = "SUBREAD",
                    runtime_attr_override = { 'mem_gb': 64 }
            }
        }

        # merge the aligned per-shard BAM into one, corresponding to one raw input BAM
        call Utils.MergeBams as MergeUncorrected { input: bams = AlignUncorrected.aligned_bam }

        # compute alignment metrics
        call AM.AlignedMetrics as PerFlowcellMetrics {
            input:
                aligned_bam    = MergeUncorrected.merged_bam,
                aligned_bai    = MergeUncorrected.merged_bai,
                ref_fasta      = ref_map['fasta'],
                ref_dict       = ref_map['dict'],
                ref_flat       = ref_map['flat'],
                dbsnp_vcf      = ref_map['dbsnp_vcf'],
                dbsnp_tbi      = ref_map['dbsnp_tbi'],
                metrics_locus  = ref_map['metrics_locus'],
                gcs_output_dir = outdir + "/metrics/per_flowcell/" + ID
        }
    }

    # gather across (potential multiple) input raw BAMs
    if (length(bams) > 1) {
        call Utils.MergeBams as MergeAllUncorrected { input: bams = MergeUncorrected.merged_bam, prefix = "~{participant_name}.uncorrected" }
    }

    File uncorrected_bam = select_first([ MergeAllUncorrected.merged_bam, MergeUncorrected.merged_bam[0] ])
    File uncorrected_bai = select_first([ MergeAllUncorrected.merged_bai, MergeUncorrected.merged_bai[0] ])

    # compute alignment metrics
    call AM.AlignedMetrics as PerSampleMetrics {
        input:
            aligned_bam    = uncorrected_bam,
            aligned_bai    = uncorrected_bai,
            ref_fasta      = ref_map['fasta'],
            ref_dict       = ref_map['dict'],
            ref_flat       = ref_map['flat'],
            dbsnp_vcf      = ref_map['dbsnp_vcf'],
            dbsnp_tbi      = ref_map['dbsnp_tbi'],
            metrics_locus  = ref_map['metrics_locus'],
            gcs_output_dir = outdir + "/metrics/combined/" + participant_name
    }

    # call SVs
    call SV.CallSVs as CallSVs {
        input:
            bam               = uncorrected_bam,
            bai               = uncorrected_bai,

            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
            tandem_repeat_bed = ref_map['tandem_repeat_bed'],

            preset            = "hifi"
    }

    # call SNVs and small indels
    call SMV.CallSmallVariants as CallSmallVariants {
        input:
            bam               = uncorrected_bam,
            bai               = uncorrected_bai,

            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
            ref_dict          = ref_map['dict'],
    }

    ##########
    # store the results into designated bucket
    ##########

    call FF.FinalizeToDir as FinalizeSVs {
        input:
            files = [ CallSVs.pbsv_vcf, CallSVs.sniffles_vcf, CallSVs.svim_vcf, CallSVs.cutesv_vcf ],
            outdir = outdir + "/variants"
    }

    call FF.FinalizeToDir as FinalizeSmallVariants {
        input:
            files = [ CallSmallVariants.longshot_vcf, CallSmallVariants.longshot_tbi ],
            outdir = outdir + "/variants"
    }

    call FF.FinalizeToDir as FinalizeMergedRuns {
        input:
            files = [ uncorrected_bam, uncorrected_bai ],
            outdir = outdir + "/alignments"
    }
}
