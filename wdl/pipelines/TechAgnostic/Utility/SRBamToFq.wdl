version 1.0

import "../../../tasks/Utility/SRUtils.wdl" as SRUtils
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow SRBamToFq {
    input {
        File bam
        File? bam_index

        File? reference_fasta
        File? reference_fasta_index
        File? reference_dict

        String participant_name

        String? gcs_out_root_dir
    }

    call SRUtils.BamToFq { 
        input: 
            bam = bam, 
            bam_index = bam_index, 
            reference_fasta = reference_fasta, 
            reference_fasta_index = reference_fasta_index, 
            reference_dict = reference_dict, 
            prefix = participant_name 
    }

    if (defined(gcs_out_root_dir)) {
        String outdir = sub(select_first([gcs_out_root_dir]), "/$", "") + "/SRBamToFq/~{participant_name}"

        call FF.FinalizeToFile as FinalizeFqEnd1 { input: outdir = outdir, file = BamToFq.fq_end1 }
        call FF.FinalizeToFile as FinalizeFqEnd2 { input: outdir = outdir, file = BamToFq.fq_end2 }
        call FF.FinalizeToFile as FinalizeFqUnpaired { input: outdir = outdir, file = BamToFq.fq_unpaired }
    }

    output {
        File fq_end1 = select_first([FinalizeFqEnd1.gcs_path, BamToFq.fq_end1])
        File fq_end2 = select_first([FinalizeFqEnd2.gcs_path, BamToFq.fq_end2])
        File fq_unpaired = select_first([FinalizeFqUnpaired.gcs_path, BamToFq.fq_unpaired])
    }
}
