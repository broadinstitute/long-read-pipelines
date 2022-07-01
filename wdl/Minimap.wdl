version 1.0

import "tasks/AlignReads.wdl" as AR

workflow CallMinimap {
    input {
        Array[File] reads
        String map_preset
        String prefix
        File ref_map_file
    }

    parameter_meta {
        map_preset: "Options include: map-ont, map-hifi, map-pb, asm5, asm10, asm20"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    call AR.Minimap2 as Align {
        input:
            reads = reads,
            ref_fasta = ref_map['fasta'],
            RG = prefix,
            map_preset = map_preset
    }

    output {
        File aligned_bam = Align.aligned_bam
        File aligned_bai = Align.aligned_bai
    }
}
