version 1.0

############################################################################################
## A workflow that merges SV VCFs for downstream evaluation.
############################################################################################

import "tasks/VariantUtils.wdl"
import "tasks/Truvari.wdl"
import "tasks/SVTK.wdl"
import "tasks/Finalize.wdl" as FF

workflow LRMergeSVVCFs {
    input {
        Array[File] vcfs
        Array[File] tbis
        File ref_map_file

        String prefix
        String caller

        String gcs_out_root_dir
    }

    parameter_meta {
        vcfs:             "GCS paths to VCF files"
        tbis:             "GCS paths to VCF tbi files"
        ref_map_file:     "table indicating reference sequence and auxillary file locations"
        prefix:           "prefix for output joint-called gVCF and tabix index"
        caller:           "SV caller whose output we're standardizing"
        gcs_out_root_dir: "GCS bucket to store the reads, variants, and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/LRMergeSVVCFs/~{prefix}"

    Map[String, String] ref_map = read_map(ref_map_file)

    call VariantUtils.MergeVCFs {
        input:
            vcfs = vcfs,
            tbis = tbis,
            prefix = prefix
    }

    call Truvari.Collapse {
        input:
            vcf = MergeVCFs.merged_vcf,
            tbi = MergeVCFs.merged_tbi,
            ref_fasta = ref_map['fasta'],
            prefix = prefix
    }

    call SVTK.Standardize {
        input:
            vcf = Collapse.collapsed_vcf,
            tbi = Collapse.collapsed_tbi,
            ref_fai = ref_map['fai'],
            prefix = prefix,
            caller = caller
    }

    # Finalize
    call FF.FinalizeToFile as FinalizeStandardizedVCF { input: outdir = outdir, file = Standardize.standardized_vcf }
    call FF.FinalizeToFile as FinalizeStandardizedTBI { input: outdir = outdir, file = Standardize.standardized_tbi }

    ##########
    # store the results into designated bucket
    ##########

    output {
        File standardized_vcf = FinalizeStandardizedVCF.gcs_path
        File standardized_tbi = FinalizeStandardizedTBI.gcs_path
    }
}
