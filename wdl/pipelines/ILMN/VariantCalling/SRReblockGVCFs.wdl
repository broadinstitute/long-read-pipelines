version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/VariantCalling/SRJointGenotyping.wdl" as SRJOINT
import "../../../tasks/VariantCalling/HaplotypeCaller.wdl" as HC
import "../../../tasks/Utility/Finalize.wdl" as FF
import "../../../tasks/Utility/VariantUtils.wdl" as VARUTIL



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
    
        String participant_name
        String gcs_out_root_dir

        String mito_contig = "chrM"
        Boolean call_vars_on_mitochondria = true
        Array[String] contigs_names_to_ignore = ["RANDOM_PLACEHOLDER_VALUE"]  ## Required for ignoring any filtering - this is kind of a hack - TODO: fix the task!
    }
    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/SRReblockGVCFs/~{prefix}"
    
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
            prefix = participant_name + ".reblocked"
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
            prefix = participant_name,
    }

    #Make sure our sample name is correct
    call VARUTIL.RenameSingleSampleVcf as RenameReblockGvcf {
        input:
            vcf = ReblockRawGVCF.output_gvcf,
            vcf_index = ReblockRawGVCF.output_gvcf_index,
            prefix = participant_name + ".reblocked.renamed",
            is_gvcf=true,
            new_sample_name = participant_name
    }

    call VARUTIL.RenameSingleSampleVcf as RenameReblockVcf {
        input:
            vcf = CollapseGVCFtoVCF.output_vcf,
            vcf_index = CollapseGVCFtoVCF.output_vcf_index,
            prefix = participant_name + ".reblocked.renamed",
            new_sample_name = participant_name
    }

    call FF.FinalizeToFile as FinalizeReblockGVCF { input: outdir = outdir, file = RenameReblockGvcf.new_sample_name_vcf }
    call FF.FinalizeToFile as FinalizeReblockGVCFIndex { input: outdir = outdir, file = RenameReblockGvcf.new_sample_name_vcf_index }
    call FF.FinalizeToFile as FinalizeReblockVCF { input: outdir = outdir, file = RenameReblockVcf.new_sample_name_vcf }
    call FF.FinalizeToFile as FInalizeReblockVCFIndex { input: outdir = outdir, file = RenameReblockVcf.new_sample_name_vcf_index }

    output {
        File output_gvcf = FinalizeReblockGVCF.gcs_path
        File output_gvcf_index = FinalizeReblockGVCFIndex.gcs_path
        File output_vcf = FinalizeReblockVCF.gcs_path
        File output_vcf_index = FInalizeReblockVCFIndex.gcs_path
        Float reblock_frac_mem = size(ReblockRawGVCF.output_gvcf) / size(gvcf)
    }

}