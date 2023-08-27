version 1.0

import "../../../tasks/Utility/PBUtils.wdl" as PB
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow PBIndex {

    meta {
        description: "The workflow generates a .pbi index for a PacBio .bam file."
    }
    parameter_meta {
        bam:                "GCS path to bam"
        pbi:                "GCS path to pbi index for bam"
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        File bam
        String gcs_out_root_dir
    }

    String dir_prefix = basename(bam, ".bam")
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBIndex/~{dir_prefix}"

    call PB.PBIndex as IndexReads { input: bam = bam }

    # Finalize data
    call FF.FinalizeToFile as FinalizePbi {
        input:
            outdir = outdir,
            file = IndexReads.pbi,
            name = basename(bam) + ".pbi"
    }

    output {
        File pbi = FinalizePbi.gcs_path
    }
}
