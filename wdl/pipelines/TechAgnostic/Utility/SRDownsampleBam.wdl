version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow SRDownsampleBam {

    meta {
        author: "Raphael Brosula"
        description: "Downsample a bam file with a given probability."
    }

    parameter_meta {
        bam: "Bam file for which to create an index."
        probability: "Probability that a read will be emitted (default = 0.01)."

        dir_prefix: "Directory prefix to use for finalized location."
        gcs_out_root_dir: "GCS Bucket into which to finalize outputs."
    }

    input {
        File bam
        Float probability

        String dir_prefix
        String gcs_out_root_dir
    }
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/SRDownsampleBam/~{dir_prefix}"
    String ds_bam_prefix = "~{dir_prefix}.downsampled_reads"

    call Utils.DownsampleSam { input: bam = bam, probability = probability, prefix = ds_bam_prefix }
    Float ds_bam_size = size(DownsampleSam.output_bam, "GB")

    call FF.FinalizeToFile as FinalizeDownsampledBam { input: outdir = outdir, file = DownsampleSam.output_bam }
    call FF.FinalizeToFile as FinalizeDownsampledBamIndex { input: outdir = outdir, file = DownsampleSam.output_bam_index }

    output {
        File downsampled_bam = FinalizeDownsampledBam.gcs_path
        File downsampled_bai = FinalizeDownsampledBamIndex.gcs_path

        Float downsampled_bam_size = ds_bam_size
        Float downsampled_factor = ds_bam_size / size(bam, "GB")
    }
}
