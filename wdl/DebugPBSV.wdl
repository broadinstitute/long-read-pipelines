version 1.0

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/CallSVs.wdl" as SV
import "tasks/Finalize.wdl" as FF

workflow DebugPBSV {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai
        File tandem_repeat_bed
    }

    call SV.PBSVDiscover {
        input:
            bam = bam,
            bai = bai,

            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            tandem_repeat_bed = tandem_repeat_bed,

            prefix = "test"
    }

    call SV.PBSVCall {
        input:
            bam = bam,
            bai = bai,

            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            svsig = PBSVDiscover.svsig,

            prefix = "test"
    }
}
