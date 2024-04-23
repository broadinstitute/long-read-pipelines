version 1.0

# import "../../../structs/Structs.wdl"
import "../../../tasks/Utility/Utils.wdl"
# import "../../../tasks/Utility/SRUtils.wdl" as SRUTIL
import "../../../tasks/VariantCalling/SRJointGenotyping.wdl" as SRJOINT
import "../../../tasks/VariantCalling/HaplotypeCaller.wdl" as HC



workflow ReblockGVCF {
    meta {
        author: "Raphael Brosula"
        description: "A workflow for reblocking GVCFs, such as those produced by HaplotypeCaller."
    }

    input {
        File gvcf
        File gvcf_index

        File ref_map_file
        File? dbsnp_vcf
    
        String prefix

        String mito_contig = "chrM"
        Boolean call_vars_on_mitochondria = true
        Array[String] contigs_names_to_ignore = ["RANDOM_PLACEHOLDER_VALUE"]  ## Required for ignoring any filtering - this is kind of a hack - TODO: fix the task!
    }
    Map[String, String] ref_map = read_map(ref_map_file)
    
    # Make sure interval list is concordant from HaplotypeCaller
    Array[String] use_filter = if (call_vars_on_mitochondria) then contigs_names_to_ignore else flatten([[mito_contig], contigs_names_to_ignore])
    
    call Utils.MakeChrIntervalList as SmallVariantsScatterPrep {
        input:
            ref_dict = ref_map['dict'],
            filter = use_filter
    }

    call HC.ReblockGVCF as ReblockRawGVCF {
        input:
            gvcf = gvcf,
            gvcf_index = gvcf_index,
            ref_fasta = ref_map['fasta'],
            ref_fasta_fai = ref_map['fai'],
            ref_dict = ref_map['dict'],
            prefix = prefix
    }

    # Collapse the reblocked GVCF into a regular VCF
    call SRJOINT.GenotypeGVCFs as CollapseGVCFtoVCF {
        input:
            input_gvcf_data = ReblockRawGVCF.output_gvcf,
            input_gvcf_index = ReblockRawGVCF.output_gvcf_index,
            interval_list = SmallVariantsScatterPrep.interval_list,
            ref_fasta = ref_map['fasta'],
            ref_fasta_fai = ref_map['fai'],
            ref_dict = ref_map['dict'],
            dbsnp_vcf = dbsnp_vcf,
            prefix = prefix,
    }

    output {
        File output_gvcf = ReblockRawGVCF.output_gvcf
        File output_gvcf_index = ReblockRawGVCF.output_gvcf_index
        File output_vcf = CollapseGVCFtoVCF.output_vcf
        File output_vcf_index = CollapseGVCFtoVCF.output_vcf_index
        Float reblock_frac_mem = size(ReblockRawGVCF.output_gvcf) / size(gvcf)
    }

}