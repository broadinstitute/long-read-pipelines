version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/PBUtils.wdl" as PB
import "tasks/HiFi.wdl" as HIFI
import "tasks/AlignReads.wdl" as AR
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/AnnotateAdapters.wdl" as AA
import "tasks/Figures.wdl" as FIG
import "tasks/Finalize.wdl" as FF

workflow PB10xSingleFlowcell {
    input {
        String gcs_input_dir
        String? sample_name
        Int fastq_shards = 50

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File ref_flat
        File dbsnp_vcf
        File dbsnp_tbi

        File metrics_locus

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PB.FindBams { input: gcs_input_dir = gcs_input_dir }

    scatter (subread_bam in FindBams.subread_bams) {
        call PB.GetRunInfo { input: subread_bam = subread_bam }

        String SM  = select_first([sample_name, GetRunInfo.run_info["SM"]])
        String PL  = "PACBIO"
        String PU  = GetRunInfo.run_info["PU"]
        String DT  = GetRunInfo.run_info["DT"]
        String ID  = PU
        String DS  = GetRunInfo.run_info["DS"]
        String DIR = SM + "." + ID

        String rg_subreads  = "@RG\\tID:~{ID}.subreads\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"
        String rg_consensus = "@RG\\tID:~{ID}.consensus\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}\\tDT:~{DT}"

        call PB.ShardLongReads { input: unmapped_files = [ subread_bam ], num_reads_per_split = 2000000 }

        scatter (subreads in ShardLongReads.unmapped_shards) {
            call PB.CCS { input: subreads = subreads }

#            call Utils.FastaToSam as FastaToSam { input: fasta = Correct.consensus }
#
#            call AA.AnnotateAdapters { input: bam = FastaToSam.output_bam }
#
#            call AR.Minimap2 as AlignSubreads {
#                input:
#                    reads = [ subreads ],
#                    ref_fasta = ref_fasta,
#                    RG = rg_subreads,
#                    map_preset = "splice"
#            }
#
#            call AR.Minimap2 as AlignConsensus {
#                input:
#                    reads = [ AnnotateAdapters.annotated_fq ],
#                    ref_fasta = ref_fasta,
#                    RG = rg_consensus,
#                    map_preset = "splice"
#            }
        }
    }

#    call PD.PrepareData as PrepareData { input: gcs_input_dir = gcs_input_dir, sample_name = sample_name }
#
#    scatter (manifest_chunk in PrepareData.manifest_chunks) {
#        call PR.ProcessReads as ProcessReads {
#            input:
#                manifest_chunk = manifest_chunk,
#                ref_fasta = ref_fasta,
#                run_info = PrepareData.run_info
#        }
#
#        call AA.AnnotateAdapters { input: bam = FastaToSam.output_bam }
#    }
#
#    call PR.FastqsToUnmappedBam as FastqsChunkToUnmappedBam {
#        input:
#            manifest_chunks = PrepareData.manifest_chunks,
#            SM=PrepareData.run_info['SM'],
#            ID=PrepareData.run_info['ID'],
#            PL=PrepareData.run_info['PL']
#    }
#
#    call PR.MergeBams as MergeRemaining {
#        input:
#            shards = ProcessReads.remaining_shard,
#            merged_name="~{sample_name}.~{PrepareData.run_info['ID']}.unaligned.bam",
#    }
#
#    call PR.MergeBams as MergeAligned {
#        input:
#            shards = ProcessReads.aligned_shard,
#            merged_name="~{sample_name}.~{PrepareData.run_info['ID']}.aligned.bam",
#    }
#
#    call UM.UnalignedMetrics as PerFlowcellUnalignedMetrics {
#        input:
#            unaligned_bam = FastqsChunkToUnmappedBam.unmapped,
#
#            per = "flowcell",
#            type = "raw",
#            label = PrepareData.run_info['ID'],
#
#            gcs_output_dir = gcs_out_root_dir
#    }
#
#    call AM.AlignedMetrics as PerFlowcellAlignedMetrics {
#        input:
#            aligned_bam = MergeAligned.merged,
#            aligned_bai = MergeAligned.merged_bai,
#            ref_fasta = ref_fasta,
#            ref_dict = ref_dict,
#            ref_flat = ref_flat,
#            dbsnp = dbsnp_vcf,
#            dbsnp_tbi = dbsnp_tbi,
#            metrics_locus = metrics_locus,
#            per = "flowcell",
#            type = "raw",
#            label = PrepareData.run_info['ID'],
#            gcs_output_dir = gcs_out_root_dir
#    }
#
#    call AR.MergeBams as MergeAllSubreads  { input: bams = MergeSubreads.merged_bam,  prefix = "~{SM[0]}.~{ID[0]}.subreads"  }
#    call AR.MergeBams as MergeAllConsensus { input: bams = MergeConsensus.merged_bam, prefix = "~{SM[0]}.~{ID[0]}.consensus" }
#    call AR.MergeBams as MergeAllAnnotated { input: bams = MergeAnnotated.merged_bam, prefix = "~{SM[0]}.~{ID[0]}.annotated" }
#
#    call AM.AlignedMetrics as PerFlowcellRunSubreadMetrics {
#        input:
#            aligned_bam    = MergeAllSubreads.merged_bam,
#            aligned_bai    = MergeAllSubreads.merged_bai,
#            ref_fasta      = ref_fasta,
#            ref_dict       = ref_dict,
#            ref_flat       = ref_flat,
#            dbsnp_vcf      = dbsnp_vcf,
#            dbsnp_tbi      = dbsnp_tbi,
#            metrics_locus  = metrics_locus,
#            per            = "flowcell",
#            type           = "run",
#            label          = ID[0] + ".subreads",
#            gcs_output_dir = outdir + "/" + DIR[0]
#    }
#
#    call AM.AlignedMetrics as PerFlowcellRunConsensusMetrics {
#        input:
#            aligned_bam    = MergeAllConsensus.merged_bam,
#            aligned_bai    = MergeAllConsensus.merged_bai,
#            ref_fasta      = ref_fasta,
#            ref_dict       = ref_dict,
#            ref_flat       = ref_flat,
#            dbsnp_vcf      = dbsnp_vcf,
#            dbsnp_tbi      = dbsnp_tbi,
#            metrics_locus  = metrics_locus,
#            per            = "flowcell",
#            type           = "run",
#            label          = ID[0] + ".consensus",
#            gcs_output_dir = outdir + "/" + DIR[0]
#    }
#
#    call AM.AlignedMetrics as PerFlowcellRunAnnotatedMetrics {
#        input:
#            aligned_bam    = MergeAllAnnotated.merged_bam,
#            aligned_bai    = MergeAllAnnotated.merged_bai,
#            ref_fasta      = ref_fasta,
#            ref_dict       = ref_dict,
#            ref_flat       = ref_flat,
#            dbsnp_vcf      = dbsnp_vcf,
#            dbsnp_tbi      = dbsnp_tbi,
#            metrics_locus  = metrics_locus,
#            per            = "flowcell",
#            type           = "run",
#            label          = ID[0] + ".annotated",
#            gcs_output_dir = outdir + "/" + DIR[0]
#    }
#
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
#
#    ##########
#    # Finalize
#    ##########
#
#    call FF.FinalizeToDir as FinalizeMergedRuns {
#        input:
#            files = [ MergeAllConsensus.merged_bam, MergeAllConsensus.merged_bai,
#                      MergeAllAnnotated.merged_bam, MergeAllAnnotated.merged_bai ],
#            outdir = outdir + "/" + DIR[0] + "/alignments"
#    }
#
#    call FF.FinalizeToDir as FinalizeMergedSubreads {
#        input:
#            files = [ MergeAllSubreads.merged_bam,  MergeAllSubreads.merged_bai ],
#            outdir = outdir + "/" + DIR[0] + "/alignments"
#    }
}
