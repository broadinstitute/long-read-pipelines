version 1.0

import "tasks/TranscriptAnalysis/Preprocessing_Tasks.wdl" as TX_PRE
import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as UTILS
import "tasks/Finalize.wdl" as FF

workflow MasSeqDownsampleArrayElementBam {

    meta {
        description : "This workflow is designed to downsample reads from a given bam file containing MAS-seq extracted array elements.  This workflow will downsample the input by several factors, all of which will be produced as output."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File array_element_bam
        File ccs_array_element_bam
        String sample_name

        String gcs_out_root_dir = "gs://broad-dsde-methods-long-reads-outgoing/MasSeqDownsampleArrayElementBam"
    }

    parameter_meta {
        array_element_bam : "Bam file containing MAS-seq array elements."
        ccs_array_element_bam : "Bam file containing CCS corrected MAS-seq array elements."
        sample_name : "Name of the sample being processed."
        gcs_out_root_dir  : "GCS bucket to store the output."
    }

    call UTILS.GetCurrentTimestampString as t_01_WdlExecutionStartTimestamp { input: }
    String outdir = sub(gcs_out_root_dir, "/$", "")

    # Get the number of reads in the given bam
    call UTILS.CountBamRecords as t_02_CountBamRecords {
        input:
            bam = array_element_bam
    }
    Float num_reads = t_02_CountBamRecords.num_records
    # Get the number of ZMWs in the given bam
    call PB.CountZMWs as t_03_CountZMWs {
        input:
            bam = array_element_bam
    }

    RuntimeAttr fat_mem_attrs = object {
        mem_gb: 32
    }

    # Downsample by various factors:
    call TX_PRE.DownsampleToIsoSeqEquivalent as t_04_DownsampleCcsReadsToIsoSeqEquivalent {
        input:
            array_element_bam = ccs_array_element_bam,
            prefix = sample_name + "_downsampled_to_isoseq"
    }
    if (num_reads > 1000000 ) {
        call UTILS.DownsampleSam as t_05_DownsampleSamTo1mReads {
            input:
                bam = array_element_bam,
                probability = 1000000.0 / num_reads,
                strategy = "Chained",
                prefix = sample_name + "_downsampled_to_1m"
        }
    }
    if (num_reads > 5000000 ) {
        call UTILS.DownsampleSam as t_06_DownsampleSamTo5mReads {
            input:
                bam = array_element_bam,
                probability = 5000000.0 / num_reads,
                strategy = "Chained",
                prefix = sample_name + "_downsampled_to_5m"
        }
    }
    if (num_reads > 10000000 ) {
        call UTILS.DownsampleSam as t_07_DownsampleSamTo10mReads {
            input:
                bam = array_element_bam,
                probability = 10000000.0 / num_reads,
                strategy = "Chained",
                runtime_attr_override = fat_mem_attrs,
                prefix = sample_name + "_downsampled_to_10m"
        }
    }
    if (num_reads > 20000000 ) {
        call UTILS.DownsampleSam as t_08_DownsampleSamTo20mReads {
            input:
                bam = array_element_bam,
                probability = 20000000.0 / num_reads,
                strategy = "Chained",
                runtime_attr_override = fat_mem_attrs,
                prefix = sample_name + "_downsampled_to_20m"
        }
    }
    if (num_reads > 30000000 ) {
        call UTILS.DownsampleSam as t_09_DownsampleSamTo30mReads {
            input:
                bam = array_element_bam,
                probability = 30000000.0 / num_reads,
                strategy = "Chained",
                runtime_attr_override = fat_mem_attrs,
                prefix = sample_name + "_downsampled_to_30m"
        }
    }

    ##########
    # store the results into designated bucket
    ##########

    String base_out_dir = outdir + "/" + sample_name + "/" + t_01_WdlExecutionStartTimestamp.timestamp_string + "/"

    call FF.FinalizeToDir as t_10_FinalizeDownsampledToIsoseqReads {
        input:
            files = [t_04_DownsampleCcsReadsToIsoSeqEquivalent.downsampled_bam],
            outdir = base_out_dir + "downsampled_to_isoseq",
    }

    if (num_reads > 1000000 ) {
        call FF.FinalizeToDir as t_11_FinalizeDownsampledTo1mReads {
            input:
                files = select_all([t_05_DownsampleSamTo1mReads.output_bam,
                                   t_05_DownsampleSamTo1mReads.output_bam_index]),
                outdir = base_out_dir + "downsampled_to_01m"
        }
    }

    if (num_reads > 5000000 ) {
        call FF.FinalizeToDir as t_12_FinalizeDownsampledTo5mReads {
            input:
                files = select_all([t_06_DownsampleSamTo5mReads.output_bam,
                                   t_06_DownsampleSamTo5mReads.output_bam_index]),
                outdir = base_out_dir + "downsampled_to_05m"
        }
    }

    if (num_reads > 10000000 ) {
        call FF.FinalizeToDir as t_13_FinalizeDownsampledTo10mReads {
            input:
                files = select_all([t_07_DownsampleSamTo10mReads.output_bam,
                                   t_07_DownsampleSamTo10mReads.output_bam_index]),
                outdir = base_out_dir + "downsampled_to_10m"
        }
    }

    if (num_reads > 20000000 ) {
        call FF.FinalizeToDir as t_14_FinalizeDownsampledTo20mReads {
            input:
                files = select_all([t_08_DownsampleSamTo20mReads.output_bam,
                                   t_08_DownsampleSamTo20mReads.output_bam_index]),
                outdir = base_out_dir + "downsampled_to_20m"
        }
    }

    if (num_reads > 30000000 ) {
        call FF.FinalizeToDir as t_15_FinalizeDownsampledTo30mReads {
            input:
                files = select_all([t_09_DownsampleSamTo30mReads.output_bam,
                                   t_09_DownsampleSamTo30mReads.output_bam_index]),
                outdir = base_out_dir + "downsampled_to_30m"
        }
    }

    # Write out some metadata here:
    call FF.WriteNamedFile as t_16_WriteNumReadsFile {
        input:
            name = "Num_reads_" + num_reads,
            outdir = base_out_dir
    }
    call FF.WriteNamedFile as t_17_WriteNumZMWsFile {
        input:
            name = "Num_ZMWs_" + t_03_CountZMWs.num_zmws,
            outdir = base_out_dir
    }
}
