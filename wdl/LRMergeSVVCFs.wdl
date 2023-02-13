version 1.0

################################################################################
## A workflow that merges SV VCFs for downstream evaluation.
################################################################################

import "tasks/VariantUtils.wdl"
import "tasks/Truvari.wdl"
import "tasks/SVTK.wdl"
import "tasks/Finalize.wdl" as FF
import "LRRegenotype.wdl"


workflow LRMergeSVVCFs {
    input {
        Array[File] vcfs
        Array[File] tbis
        File bam_addresses
        File ref_map_file

        String prefix
        String caller
        
        Int n_nodes
        Int n_cpus
        Int bam_size_gb

        String gcs_out_root_dir
    }

    parameter_meta {
        vcfs:             "GCS paths to VCF files"
        tbis:             "GCS paths to VCF tbi files"
        ref_map_file:     "table indicating reference sequence and auxillary file locations"
        prefix:           "prefix for output joint-called gVCF and tabix index"
        caller:           "SV caller whose output we're standardizing"
        n_nodes:          "Use this number of nodes to regenotype in parallel."
        n_cpus:           "Lower bound on the number of CPUs per regenotype node."
        bam_size_gb:      "Upper bound on the size of a single BAM."
        gcs_out_root_dir: "GCS bucket to store the reads, variants, and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/LRMergeSVVCFs/~{prefix}"

    Map[String, String] ref_map = read_map(ref_map_file)

    call VariantUtils.MergeVCFs { input: vcfs = vcfs, tbis = tbis, reference_fa = ref_map['fasta'], prefix = prefix }

    #call VariantUtils.GetContigNames { input: vcf = MergeVCFs.merged_vcf }

    #scatter (contig_name in GetContigNames.contig_names) {
        
    #call VariantUtils.SubsetVCF { input: vcf_gz = MergeVCFs.merged_vcf, vcf_tbi = MergeVCFs.merged_tbi, locus = contig_name }
    
    call Truvari.Collapse {
        input:
            vcf = MergeVCFs.merged_vcf,
            tbi = MergeVCFs.merged_tbi,
            ref_fasta = ref_map['fasta'],
            prefix = prefix
    }
    #}
    
    #call VariantUtils.ConcatVCFs { input: vcfs = Collapse.collapsed_vcf, tbis = Collapse.collapsed_tbi, prefix = prefix }

    call LRRegenotype.LRRegenotype {
        input:
            merged_vcf_gz = Collapse.collapsed_vcf,
            bam_addresses = bam_addresses,
            use_lrcaller = 1,
            use_cutesv = 0,
            reference_fa = ref_map['fasta'],
            reference_fai = ref_map['fai'],
            n_nodes = n_nodes,
            n_cpus = n_cpus,
            bam_size_gb = bam_size_gb
    }

    call SVTK.Standardize {
        input:
            vcf = LRRegenotype.vcf_gz,
            tbi = LRRegenotype.tbi,
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
        File counts_MergeVCFs = MergeVCFs.counts
        File counts_TruvariCollapse = Collapse.counts
    }
}
