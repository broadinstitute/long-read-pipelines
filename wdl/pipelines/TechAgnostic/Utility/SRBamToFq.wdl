version 1.0

import "../../../tasks/Utility/SRUtils.wdl" as SRUtils
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow SRBamToFq {
    input {
        File bam
        String participant_name

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/SRBamToFq/~{participant_name}"

    call SRUtils.BamToFq { input: bam = bam, prefix = participant_name }

    call FF.FinalizeToFile as FinalizeFqEnd1 { input: outdir = outdir, file = BamToFq.fq_end1 }
    call FF.FinalizeToFile as FinalizeFqEnd2 { input: outdir = outdir, file = BamToFq.fq_end2 }
    call FF.FinalizeToFile as FinalizeFqUnpaired { input: outdir = outdir, file = BamToFq.fq_unpaired }

    output {
        File fq_end1 = FinalizeFqEnd1.gcs_path
        File fq_end2 = FinalizeFqEnd2.gcs_path
        File fq_unpaired = FinalizeFqUnpaired.gcs_path
    }
}
