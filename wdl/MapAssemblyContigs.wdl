version 1.0

import "tasks/AlignReads.wdl"

workflow MapAssemblyContigs {
    input {
        File contigs_fa
        String read_group
        String sample
        File ref_map_file
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    call AlignReads.Minimap2 { input:
        reads = [contigs_fa],
        ref_fasta = ref_map['fasta'],
        RG = "@RG\\tID:~{read_group}\\tSM:~{sample}",
        map_preset = "asm5"
    }
    output {
        File aligned_bam = Minimap2.aligned_bam
        File aligned_bai = Minimap2.aligned_bai
    }
}