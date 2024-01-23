version 1.0

import "../../../tasks/Utility/Hail.wdl" as Hail
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ConvertToHailMTT2T {

    meta {
        description: "Convert a VCF to a Hail MatrixTable"
    }
    parameter_meta {
    }

    input {
        File whole_genome_vcf
        File whole_genome_vcf_tbi
        String prefix
        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/Hail/~{prefix}"

    call Hail.ConvertToHailMT as RunConvertToHailMT {
        input:
            gvcf = whole_genome_vcf,
            tbi = whole_genome_vcf_tbi,
            prefix = prefix,
            outdir = outdir

    }
    

    output {
        String joint_mt = RunConvertToHailMT.gcs_path
    }
}

