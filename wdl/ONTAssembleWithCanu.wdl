version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/AsmUtils.wdl" as AU
import "tasks/Canu.wdl" as Canu
import "tasks/Medaka.wdl" as Medaka
import "tasks/Quast.wdl" as Quast
import "tasks/Finalize.wdl" as FF

workflow ONTAssembleWithCanu {
    input {
        String gcs_fastq_dir

        File ref_map_file

        Float correct_error_rate = 0.15
        Float trim_error_rate = 0.15
        Float assemble_error_rate = 0.15
        String medaka_model = "r941_prom_high_g360"

        String prefix

        String gcs_out_root_dir
    }

    Map[String, String] ref_map = read_map(ref_map_file)
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTAssembleWithCanu/~{prefix}"

    call Utils.ComputeGenomeLength { input: fasta = ref_map['fasta'] }

    call Utils.ListFilesOfType { input: gcs_dir = gcs_fastq_dir, suffixes = [".fastq", ".fq", ".fastq.gz", ".fq.gz"] }
    call AU.MergeFastqs { input: fastqs = ListFilesOfType.files }

    call Canu.Canu {
        input:
            reads = MergeFastqs.merged_fastq,
            prefix = prefix,
            genome_size = ComputeGenomeLength.length,
            correct_error_rate = correct_error_rate,
            trim_error_rate = trim_error_rate,
            assemble_error_rate = assemble_error_rate
    }

    call Medaka.MedakaPolish {
        input:
            basecalled_reads = MergeFastqs.merged_fastq,
            draft_assembly = Canu.fa,
            model = medaka_model,
            prefix = basename(Canu.fa, ".fasta") + ".consensus.fasta",
            n_rounds = 3
    }

    call FF.FinalizeToFile as FinalizeAsmUnpolished {
        input:
            file    = Canu.fa,
            outfile = outdir + "/assembly/" + basename(Canu.fa)
    }

    call FF.FinalizeToFile as FinalizeAsmPolished {
        input:
            file    = MedakaPolish.polished_assembly,
            outfile = outdir + "/assembly/" + basename(MedakaPolish.polished_assembly)
    }

    output {
        File asm_unpolished = FinalizeAsmUnpolished.gcs_path
        File asm_polished = FinalizeAsmPolished.gcs_path
    }
}