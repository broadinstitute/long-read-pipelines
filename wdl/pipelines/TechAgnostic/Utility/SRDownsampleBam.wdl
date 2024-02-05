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

        outdir: "GCS Bucket into which to finalize outputs."
    }

    input {
        File bam

        Float probability = 0.01
        String outdir
    }

    call Utils.DownsampleSam { input: bam = bam, probability = probability }
    Float downsampled_bam_size = size(DownsampleSam.output_bam, "GB")

    call FF.FinalizeToFile as FinalizeDownsampledBam { input: outdir = outdir, file = DownsampleSam.output_bam }
    call FF.FinalizeToFile as FinalizeDownsampledBamIndex { input: outdir = outdir, file = DownsampleSam.output_bam_index }

    output {
        File downsampled_bam = FinalizeDownsampledBam.gcs_path
        File downsampled_bai = FinalizeDownsampledBamIndex.gcs_path
        Float downsampled_bam_size = downsampled_bam_size
    }
}
