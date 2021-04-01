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
import "tasks/Hifiasm.wdl" as HA
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/Figures.wdl" as FIG
import "tasks/CallSVs.wdl" as SV
import "tasks/CallSmallVariants.wdl" as SMV
import "tasks/CallAssemblyVariants.wdl" as AV
import "tasks/Finalize.wdl" as FF

workflow PBCCSWholeGenome {
    input {
        Array[File] bams
        Array[File] pbis

        File ref_map_file

        String participant_name
        Int num_shards = 25
        Boolean extract_uncorrected_reads = false

        String gcs_out_root_dir
    }

    parameter_meta {
        bams:                      "GCS path to CCS-corrected data"
        pbis:                      "GCS path to .pbi index for CCS-corrected data"

        ref_map_file:              "table indicating reference sequence and auxillary file locations"

        participant_name:          "name of the participant from whom these samples were obtained"
        num_shards:                "[default-valued] number of sharded BAMs to create (tune for performance)"
        extract_uncorrected_reads: "[default-valued] extract reads that were not CCS-corrected to a separate file"

        gcs_out_root_dir:          "GCS bucket to store the corrected/uncorrected reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCCSWholeGenome/" + participant_name

    # scatter over all sample BAMs
    scatter (p in zip(bams, pbis)) {
        File bam = p.left
        File pbi = p.right

        call PB.GetRunInfo { input: bam = bam, SM = participant_name }
        String ID = GetRunInfo.run_info["PU"]

        # break one raw BAM into fixed number of shards
        call PB.ShardLongReads { input: unaligned_bam = bam, unaligned_pbi = pbi, num_shards = num_shards }

        # then perform correction and alignment on each of the shard
        scatter (reads in ShardLongReads.unmapped_shards) {
            call Utils.BamToFastq { input: bam = reads, prefix = basename(reads, ".bam") }

            call PB.Align as AlignCorrected {
                input:
                    bam         = reads,
                    ref_fasta   = ref_map['fasta'],
                    sample_name = participant_name,
                    map_preset  = "CCS"
            }
        }

        # merge the corrected per-shard BAM/fastq into one, corresponding to one raw input BAM
        call Utils.MergeBams as MergeCorrected { input: bams = AlignCorrected.aligned_bam, prefix = "~{participant_name}.~{ID}.corrected" }
        call Utils.MergeFastqs { input: fqs = BamToFastq.reads_fq, prefix = "~{participant_name}.~{ID}" }

        # compute alignment metrics
        call AM.AlignedMetrics as PerFlowcellMetrics {
            input:
                aligned_bam    = MergeCorrected.merged_bam,
                aligned_bai    = MergeCorrected.merged_bai,
                ref_fasta      = ref_map['fasta'],
                ref_dict       = ref_map['dict'],
                gcs_output_dir = outdir + "/metrics/per_flowcell/" + ID
        }
    }

    # gather across (potential multiple) input raw BAMs
    if (length(bams) > 1) {
        call Utils.MergeBams as MergeAllCorrected { input: bams = MergeCorrected.merged_bam, prefix = "~{participant_name}.corrected" }
        call Utils.MergeFastqs as MergeAllFastqs { input: fqs = MergeFastqs.merged_fq, prefix = "~{participant_name}.corrected" }
    }

    File ccs_bam = select_first([ MergeAllCorrected.merged_bam, MergeCorrected.merged_bam[0] ])
    File ccs_bai = select_first([ MergeAllCorrected.merged_bai, MergeCorrected.merged_bai[0] ])
    File ccs_fq  = select_first([ MergeAllFastqs.merged_fq, MergeFastqs.merged_fq[0] ])

#    # assemble genome
#    call HA.Hifiasm {
#        input:
#            reads = ccs_fq,
#            prefix = participant_name
#    }

    # compute alignment metrics
    call AM.AlignedMetrics as PerSampleMetrics {
        input:
            aligned_bam    = ccs_bam,
            aligned_bai    = ccs_bai,
            ref_fasta      = ref_map['fasta'],
            ref_dict       = ref_map['dict'],
            gcs_output_dir = outdir + "/metrics/combined/" + participant_name
    }

    # call SVs
    call SV.CallSVs as CallSVs {
        input:
            bam               = ccs_bam,
            bai               = ccs_bai,

            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
            tandem_repeat_bed = ref_map['tandem_repeat_bed'],

            preset            = "hifi"
    }

    # call SNVs and small indels
    call SMV.CallSmallVariants as CallSmallVariants {
        input:
            bam               = ccs_bam,
            bai               = ccs_bai,

            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
            ref_dict          = ref_map['dict'],

            preset            = "hifi"
    }

#    # call variants in assemblies
#    call AV.CallAssemblyVariants {
#        input:
#            asm_fasta = Hifiasm.fa,
#            ref_fasta = ref_map['fasta'],
#            participant_name = participant_name,
#            prefix = basename(ccs_bam, ".bam") + ".hifiasm"
#    }

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
            files = [ CallSmallVariants.longshot_vcf, CallSmallVariants.longshot_tbi,
                      CallSmallVariants.deepvariant_vcf, CallSmallVariants.deepvariant_tbi,
                      CallSmallVariants.deepvariant_gvcf, CallSmallVariants.deepvariant_gtbi ],
            outdir = outdir + "/variants"
    }

#    call FF.FinalizeToDir as FinalizeAssemblyVariants {
#        input:
#            files = [ CallAssemblyVariants.paftools_vcf ],
#            outdir = outdir + "/variants"
#    }

    call FF.FinalizeToDir as FinalizeMergedRuns {
        input:
            files = [ ccs_bam, ccs_bai ],
            outdir = outdir + "/alignments"
    }

#    call FF.FinalizeToDir as FinalizeAssembly {
#        input:
#            files = [ Hifiasm.gfa, Hifiasm.fa, CallAssemblyVariants.paf ],
#            outdir = outdir + "/assembly"
#    }

    output {
#        # Assembly
#        File hifiasm_gfa = Hifiasm.gfa
#        File hifiasm_fa = Hifiasm.fa
#        File hifiasm_paf = CallAssemblyVariants.paf
#
#        # Assembly-based variants
#        File paftools_vcf = CallAssemblyVariants.paftools_vcf

        ### Ref
        ### ===
        # BAMs
        File corrected_bam = ccs_bam
        File corrected_bai = ccs_bai

        # SVs
        File pbsv_vcf = CallSVs.pbsv_vcf
        File sniffles_vcf = CallSVs.sniffles_vcf
        File svim_vcf = CallSVs.svim_vcf
        File cutesv_vcf = CallSVs.cutesv_vcf

        # SNVs/indels
        File longshot_vcf = CallSmallVariants.longshot_vcf
        File longshot_tbi = CallSmallVariants.longshot_tbi
        File deepvariant_vcf = CallSmallVariants.deepvariant_vcf
        File deepvariant_tbi = CallSmallVariants.deepvariant_tbi
        File deepvariant_gvcf = CallSmallVariants.deepvariant_gvcf
        File deepvariant_gtbi = CallSmallVariants.deepvariant_gtbi
    }
}
