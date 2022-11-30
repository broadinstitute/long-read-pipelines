version 1.0

import "../../../tasks/QC/FPCheckAoU.wdl" as FP

workflow VerifyBamFingerprint {
    meta {
        desciption:
        "Verify fingerprint of a single BAM where it's assumed the BAM holds data from a single entity"
        warn:
        "So far, we've verified that this works for CCS/Hifi and ONT data, for legacy CLR data, it's not supported."
    }

    parameter_meta {
        fp_store: "GCS storage bucket and folder where the fingperprint VCF files are stored."
        sample_id_at_store: "sample id of the data at the storage; CRITICAL: it's assumsed that the fingerprint VCF file follow the naming convention that start with this sample id."
        ref_specific_haplotype_map: "Reference-specific haplotype map file to be passed on to Picard's `CheckFingerprint`"

        lod_pass_threshold: "Threshold for LOD above which the BAM will be declared as PASSing this QC check"
        lod_fail_threshold: "Threshold for LOD below which the BAM will be declared as FAILing this QC check; LOD between the two thresholds will lead to the BAM's QC status as BORDERLINE"
    }

    input {
        File aligned_bam
        File aligned_bai

        String fp_store
        String sample_id_at_store

        File ref_specific_haplotype_map

        Float lod_pass_threshold =  6.0
        Float lod_fail_threshold = -3.0
    }

    output {
        Map[String, String] fingerprint_check = {"status": core.FP_status,
                                                 "LOD": core.lod_expected_sample}
    }

    call FP.FPCheckAoU as core {
        input:
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai,
            fp_store = fp_store,
            sample_id_at_store = sample_id_at_store,
            ref_specific_haplotype_map = ref_specific_haplotype_map,
            lod_pass_threshold = lod_pass_threshold,
            lod_fail_threshold = lod_fail_threshold
    }
}
