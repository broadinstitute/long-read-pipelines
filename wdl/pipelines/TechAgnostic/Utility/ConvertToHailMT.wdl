version 1.0

import "../../../tasks/Utility/Hail.wdl" as Hail

workflow ConvertToHailMT {

    meta {
        description: "Convert a (g)VCF to a Hail MatrixTable"
    }
    parameter_meta {
        joint_gvcf:       "joint-called gVCF file"
        joint_gvcf_tbi:   ".tbi index for joint-called gVCF file"
        prefix:           "prefix for output Hail MatrixTable"
        gcs_out_root_dir: "GCS bucket in which to store the Hail MatrixTable"
    }

    input {
        File joint_gvcf
        File joint_gvcf_tbi
        String prefix

        String gcs_out_root_dir

        String reference
        String? ref_fasta
        String? ref_fai
    }

    output {
        String joint_mt = ConvertToHailMT.gcs_path
    }

    call Hail.ConvertToHailMT { input:
        gvcf = joint_gvcf,
        tbi  = joint_gvcf_tbi,
        prefix = prefix,
        outdir = sub(gcs_out_root_dir, "/$", ""),

        reference = reference,
        ref_fasta = ref_fasta,
        ref_fai = ref_fai
    }
}
