version 1.0

import "tasks/ONTPepper.wdl"
import "tasks/Finalize.wdl" as FF

workflow PhaseBams {
    input {
        File bam
        File bai

        String sample_name

        File ref_map_file

        Int threads
        Int memory

        String out_root = "gs://broad-dsde-methods-long-reads-outgoing/aou-long-reads-ont/ONTWholeGenome"
    }
    Map[String, String] ref_map = read_map(ref_map_file)

    String dir = out_root + '/~{sample_name}/alignments/'

    call ONTPepper.Pepper {
        input:
            bam = bam, bai = bai,
            ref_fasta = ref_map['fasta'], ref_fasta_fai = ref_map['fai'],
            threads = threads, memory = memory
    }

    call FF.FinalizeToFile as FinalizeBam { input: outdir = dir, file = Pepper.hap_tagged_bam, name = "~{sample_name}.MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam" }
    call FF.FinalizeToFile as FinalizeBai { input: outdir = dir, file = Pepper.hap_tagged_bai, name = "~{sample_name}.MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai" }

    output {
        File haplotagged_bam = FinalizeBam.gcs_path
        File haplotagged_bai = FinalizeBai.gcs_path
    }
}