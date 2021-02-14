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
        File ref_map_file
        File? alt_map_file

        String participant_name
        Int num_shards = 300
        Boolean extract_uncorrected_reads = false

        String gcs_out_root_dir
    }

    parameter_meta {
        bams:                      "GCS path to raw subreads or CCS data"
        ref_map_file:              "table indicating reference sequence and auxillary file locations"
        alt_map_file:              "table indicating alternate reference sequence and auxillary file locations"

        participant_name:          "name of the participant from whom these samples were obtained"
        num_shards:                "[default-valued] number of sharded BAMs to create (tune for performance)"
        extract_uncorrected_reads: "[default-valued] extract reads that were not CCS-corrected to a separate file"

        gcs_out_root_dir:          "GCS bucket to store the corrected/uncorrected reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)
    Map[String, String] alt_map = if defined(alt_map_file) then read_map(select_first([alt_map_file])) else {'fasta': 'none'}

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCCSWholeGenome/" + participant_name

    # scatter over all sample BAMs
    scatter (bam in bams) {
        File pbi = sub(bam, ".bam$", ".bam.pbi")

        call PB.GetRunInfo { input: bam = bam }
        String ID = GetRunInfo.run_info["PU"]

        # break one raw BAM into fixed number of shards
        call PB.ShardLongReads { input: unaligned_bam = bam, unaligned_pbi = pbi, num_shards = num_shards }

        # then perform correction and alignment on each of the shard
        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads }

            call Utils.BamToFastq { input: bam = CCS.consensus, prefix = basename(CCS.consensus, ".bam") }

            if (extract_uncorrected_reads) {
                call PB.ExtractUncorrectedReads { input: subreads = subreads, consensus = CCS.consensus }

                call PB.Align as AlignUncorrected {
                    input:
                        bam         = ExtractUncorrectedReads.uncorrected,
                        ref_fasta   = ref_map['fasta'],
                        sample_name = participant_name,
                        map_preset  = "SUBREAD",
                        runtime_attr_override = { 'mem_gb': 64 }
                }
            }

            call PB.Align as AlignCorrected {
                input:
                    bam         = CCS.consensus,
                    ref_fasta   = ref_map['fasta'],
                    sample_name = participant_name,
                    map_preset  = "CCS"
            }

            if (defined(alt_map_file)) {
                call PB.Align as AlignAltCorrected {
                    input:
                        bam         = CCS.consensus,
                        ref_fasta   = alt_map['fasta'],
                        sample_name = participant_name,
                        map_preset  = "CCS"
                }
            }
        }

        # merge the corrected per-shard BAM/report/fastq into one, corresponding to one raw input BAM
        call Utils.MergeBams as MergeCorrected { input: bams = AlignCorrected.aligned_bam, prefix = "~{participant_name}.~{ID}.corrected" }

        if (length(select_all(AlignUncorrected.aligned_bam)) > 0) {
            call Utils.MergeBams as MergeUncorrected {
                input:
                    bams = select_all(AlignUncorrected.aligned_bam),
                    prefix = "~{participant_name}.~{ID}.uncorrected"
            }
        }

        if (length(select_all(AlignAltCorrected.aligned_bam)) > 0) {
            call Utils.MergeBams as MergeAltCorrected {
                input:
                    bams = select_all(AlignAltCorrected.aligned_bam),
                    prefix = "~{participant_name}.~{ID}.alt.corrected"
            }
        }

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

        call PB.MergeCCSReports as MergeCCSReports { input: reports = CCS.report, prefix = "~{participant_name}.~{ID}" }

        call FF.FinalizeToDir as FinalizeCCSReport {
            input:
                files = [ MergeCCSReports.report ],
                outdir = outdir + "/metrics/per_flowcell/" + ID + "/ccs_metrics"
        }
    }

    # gather across (potential multiple) input raw BAMs
    if (length(bams) > 1) {
        call Utils.MergeBams as MergeAllCorrected { input: bams = MergeCorrected.merged_bam, prefix = "~{participant_name}.corrected" }
        call PB.MergeCCSReports as MergeAllCCSReports { input: reports = MergeCCSReports.report }

        if (length(select_all(MergeUncorrected.merged_bam)) > 0) {
            call Utils.MergeBams as MergeAllUncorrected {
                input:
                    bams = select_all(MergeUncorrected.merged_bam),
                    prefix = "~{participant_name}.uncorrected"
            }
        }

        if (length(select_all(MergeAltCorrected.merged_bam)) > 0) {
            call Utils.MergeBams as MergeAllAltCorrected {
                input:
                    bams = select_all(MergeAltCorrected.merged_bam),
                    prefix = "~{participant_name}.alt.corrected"
            }
        }

        call Utils.MergeFastqs as MergeAllFastqs { input: fqs = MergeFastqs.merged_fq, prefix = "~{participant_name}.corrected" }
    }

    File ccs_bam = select_first([ MergeAllCorrected.merged_bam, MergeCorrected.merged_bam[0] ])
    File ccs_bai = select_first([ MergeAllCorrected.merged_bai, MergeCorrected.merged_bai[0] ])
    File ccs_report = select_first([ MergeAllCCSReports.report, MergeCCSReports.report[0] ])

    if (extract_uncorrected_reads) {
        File? uncorrected_bam = select_first([ MergeAllUncorrected.merged_bam, MergeUncorrected.merged_bam[0] ])
        File? uncorrected_bai = select_first([ MergeAllUncorrected.merged_bai, MergeUncorrected.merged_bai[0] ])
    }

    if (defined(alt_map_file)) {
        File? alt_bam = select_first([ MergeAllAltCorrected.merged_bam, MergeAltCorrected.merged_bam[0] ])
        File? alt_bai = select_first([ MergeAllAltCorrected.merged_bai, MergeAltCorrected.merged_bai[0] ])
    }

    File ccs_fq = select_first([ MergeAllFastqs.merged_fq, MergeFastqs.merged_fq[0] ])

    # assemble genome
    call HA.Hifiasm {
        input:
            reads = ccs_fq,
            prefix = participant_name
    }

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

    # call variants in assemblies
    call AV.CallAssemblyVariants {
        input:
            asm_fasta = Hifiasm.fa,
            ref_fasta = ref_map['fasta'],
            participant_name = participant_name,
            prefix = basename(ccs_bam, ".bam") + ".hifiasm"
    }

    if (defined(alt_bam) && defined(alt_bai)) {
        # call SVs
        call SV.CallSVs as CallAltSVs {
            input:
                bam               = select_first([alt_bam]),
                bai               = select_first([alt_bai]),

                ref_fasta         = alt_map['fasta'],
                ref_fasta_fai     = alt_map['fai'],

                preset            = "hifi"
        }

        # call SNVs and small indels
        call SMV.CallSmallVariants as CallAltSmallVariants {
            input:
                bam               = select_first([alt_bam]),
                bai               = select_first([alt_bai]),

                ref_fasta         = alt_map['fasta'],
                ref_fasta_fai     = alt_map['fai'],
                ref_dict          = alt_map['dict'],

                preset            = "hifi"
        }

        call FF.FinalizeToDir as FinalizeAltSVs {
            input:
                files = [ CallAltSVs.pbsv_vcf, CallAltSVs.sniffles_vcf, CallAltSVs.svim_vcf, CallAltSVs.cutesv_vcf ],
                outdir = outdir + "/alt_variants"
        }

        call FF.FinalizeToDir as FinalizeAltSmallVariants {
            input:
                files = [ CallAltSmallVariants.longshot_vcf, CallAltSmallVariants.longshot_tbi,
                          CallAltSmallVariants.deepvariant_vcf, CallAltSmallVariants.deepvariant_tbi,
                          CallAltSmallVariants.deepvariant_gvcf, CallAltSmallVariants.deepvariant_gtbi ],
                outdir = outdir + "/alt_variants"
        }
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
            files = [ CallSmallVariants.longshot_vcf, CallSmallVariants.longshot_tbi,
                      CallSmallVariants.deepvariant_vcf, CallSmallVariants.deepvariant_tbi,
                      CallSmallVariants.deepvariant_gvcf, CallSmallVariants.deepvariant_gtbi ],
            outdir = outdir + "/variants"
    }

    call FF.FinalizeToDir as FinalizeAssemblyVariants {
        input:
            files = [ CallAssemblyVariants.paftools_vcf ],
            outdir = outdir + "/variants"
    }

    call FF.FinalizeToDir as FinalizeMergedRuns {
        input:
            files = [ ccs_bam, ccs_bai ],
            outdir = outdir + "/alignments"
    }

    call FF.FinalizeToDir as FinalizeMergedCCSReport {
        input:
            files = [ ccs_report ],
            outdir = outdir + "/metrics/combined/" + participant_name + "/ccs_metrics"
    }

    if (extract_uncorrected_reads) {
        call FF.FinalizeToDir as FinalizeMergedUncorrectedRuns {
            input:
                files = select_all([ uncorrected_bam, uncorrected_bai ]),
                outdir = outdir + "/alignments"
        }
    }

    call FF.FinalizeToDir as FinalizeAssembly {
        input:
            files = [ Hifiasm.gfa, Hifiasm.fa, CallAssemblyVariants.paf ],
            outdir = outdir + "/assembly"
    }

    output {
        File corrected_bam = ccs_bam
        File corrected_bai = ccs_bai
    }
}
