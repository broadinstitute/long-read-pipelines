version 1.0

import "../../../structs/Structs.wdl"
import "../../../tasks/TertiaryAnalysis/HmmIBDTasks.wdl" as HMMIBD

workflow FilterVcfForHmmIBD {

    meta {
        author: "Jonn Smith"
        description: "Standalone filtration of a VCF/BCF for IBD analysis, runnable on its own or called by the HmmIBD workflow. Masks low-depth genotypes, (optionally) preserves original allele frequencies, recomputes AN/AC/AF, selects variant types / biallelic sites, and optionally thins to a variant cap. The filtered callset is emitted in the requested output_format: a compressed BCF (default), a bgzipped VCF, or the hmmibd-rs text genotype table (+ matching allele-frequency file)."

        outputs: {
            filtered_bcf:       "Filtered, AN/AC/AF-recomputed callset as a compressed BCF (always produced)",
            filtered_bcf_index: "CSI index for filtered_bcf",
            filtered_vcf:       "Filtered callset as a bgzipped VCF; only when output_format='vcf'",
            filtered_vcf_index: "Tabix index for filtered_vcf; only when output_format='vcf'",
            hmmibd_gt_table:    "hmmIBD text genotype table; only when output_format='hmmibd_table'",
            hmmibd_freq_table:  "Matching bi-allelic allele-frequency file; only when output_format='hmmibd_table'"
        }
    }

    parameter_meta {
        input_vcf:           "VCF/BCF to filter prior to IBD inference. (required)"
        prefix:              "Basename / sample-set name for all outputs. (required)"
        output_format:       "Format to emit: 'bcf' (default), 'vcf', or 'hmmibd_table' (hmmibd-rs text genotype + freq tables). The BCF is always produced; vcf/table are additional conversions of it."

        min_depth:           "Genotypes with FORMAT/DP below this value are set to missing. (default: 5)"
        keep_original_af:    "Preserve the original allele-frequency annotations (via rename_annots_tsv) before recomputing AN/AC/AF. (default: false)"
        variant_types:       "Which variant types to keep: 'snps', 'indels', or 'both'. (default: snps)"
        biallelic_only:      "Restrict to biallelic sites. (default: true)"
        split_multiallelics: "Split multiallelics into biallelic records to recover SNP alleles. (default: true)"
        max_variants:        "Optional cap; when >0 the callset is thinned evenly to at most this many variants. (default: 0 = no limit)"
        populations_file:    "Optional sample-to-population file for per-population AN/AC/AF. (default: none)"
        rename_annots_tsv:   "Required when keep_original_af=true: old-name<TAB>new-name TSV. (default: none)"
        filter_extra_args:   "Additional args appended verbatim to the final bcftools view invocation. (default: empty)"

        filter_runtime_attr_override:  "Override runtime for the filtration task. (default: none)"
        convert_runtime_attr_override: "Override runtime for the vcf/table conversion task. (default: none)"
    }

    input {
        File input_vcf
        String prefix
        String output_format = "bcf"

        Int min_depth = 5
        Boolean keep_original_af = false
        String variant_types = "snps"
        Boolean biallelic_only = true
        Boolean split_multiallelics = true
        Int max_variants = 0
        File? populations_file
        File? rename_annots_tsv
        String filter_extra_args = ""

        RuntimeAttr? filter_runtime_attr_override
        RuntimeAttr? convert_runtime_attr_override
    }

    # Filtration always produces a BCF.
    call HMMIBD.FilterVcfForHmmIBD as t_01_Filter {
        input:
            input_vcf           = input_vcf,
            prefix              = prefix,
            min_depth           = min_depth,
            keep_original_af    = keep_original_af,
            variant_types       = variant_types,
            biallelic_only      = biallelic_only,
            split_multiallelics = split_multiallelics,
            max_variants        = max_variants,
            populations_file    = populations_file,
            rename_annots_tsv   = rename_annots_tsv,
            extra_args          = filter_extra_args,
            runtime_attr_override = filter_runtime_attr_override
    }

    # Optional conversions of the filtered BCF.
    if (output_format == "vcf") {
        call HMMIBD.BcfToVcf as t_02_ToVcf {
            input:
                input_bcf = t_01_Filter.filtered_bcf,
                prefix    = prefix,
                runtime_attr_override = convert_runtime_attr_override
        }
    }

    if (output_format == "hmmibd_table") {
        call HMMIBD.BcfToHmmIBDTable as t_02_ToTable {
            input:
                input_bcf = t_01_Filter.filtered_bcf,
                prefix    = prefix,
                runtime_attr_override = convert_runtime_attr_override
        }
    }

    output {
        File filtered_bcf        = t_01_Filter.filtered_bcf
        File filtered_bcf_index  = t_01_Filter.filtered_bcf_index
        File? filtered_vcf       = t_02_ToVcf.filtered_vcf
        File? filtered_vcf_index = t_02_ToVcf.filtered_vcf_index
        File? hmmibd_gt_table    = t_02_ToTable.hmmibd_gt_table
        File? hmmibd_freq_table  = t_02_ToTable.hmmibd_freq_table
    }
}
