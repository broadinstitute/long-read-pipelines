version 1.0

import "../../../tasks/Alignment/AlignReads.wdl" as AR
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow AlignAssemblies {

    input {
        File fa_gz
        File ref_map_file

        String prefix

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/AlignAssemblies/~{prefix}"

    Map[String, String] ref_map = read_map(ref_map_file)

    call AR.Minimap2 {
        input:
            reads      = [ fa_gz ],
            ref_fasta  = ref_map['fasta'],
            RG         = "@RG\\tID:~{prefix}\\tSM:~{prefix}",
            map_preset = "asm10"
    }

    # Finalize
    call FF.FinalizeToFile as FinalizeBAM { input: outdir = outdir, file = Minimap2.aligned_bam }
    call FF.FinalizeToFile as FinalizeBAI { input: outdir = outdir, file = Minimap2.aligned_bai }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File aligned_bam = FinalizeBAM.gcs_path
        File aligned_bai = FinalizeBAI.gcs_path
    }
}
