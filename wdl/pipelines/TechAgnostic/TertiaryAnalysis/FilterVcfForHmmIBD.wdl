version 1.0

import "../../../structs/Structs.wdl"
import "../../../tasks/TertiaryAnalysis/HmmIBDTasks.wdl" as HMMIBD

workflow FilterVcfForHmmIBD {

    meta {
        author: "Jonn Smith"
        description: "Filter a LIST of VCFs/BCFs for IBD analysis and (for the hmmIBD-table format) combine them into a SINGLE input for hmmibd-rs. Each input file is filtered independently (scattered): mask low-depth genotypes, optionally preserve original allele frequencies, recompute AN/AC/AF, select variant types / biallelic sites. The list is meant to be disjoint regions of ONE joint call, split across files so that no multi-hundred-GB VCF is ever handled whole. Output depends on output_format: 'hmmibd_table' (recommended at scale) reduces each filtered file to a compact integer genotype table and STITCHES them into one coordinate-sorted table + allele-frequency file (a merged BCF is never materialized); 'bcf' (default) and 'vcf' return the per-file filtered callsets as arrays (not combined). Pass a single file as a one-element array."

        outputs: {
            filtered_bcfs:        "Per-input filtered, AN/AC/AF-recomputed callsets as compressed BCFs (always produced, one per input VCF)",
            filtered_bcf_indices: "CSI indices for filtered_bcfs",
            filtered_vcfs:        "Per-input filtered callsets as bgzipped VCFs; only when output_format='vcf'",
            filtered_vcf_indices: "Tabix indices for filtered_vcfs; only when output_format='vcf'",
            sample_gt_table:      "SINGLE combined hmmIBD text genotype table across all inputs; only when output_format='hmmibd_table'",
            sample_freq_table:    "Matching bi-allelic allele-frequency file for the combined table; only when output_format='hmmibd_table'"
        }
    }

    parameter_meta {
        input_vcfs:          "VCFs/BCFs to filter prior to IBD inference. Meant to be disjoint regions of one joint call (same samples across files). Pass a single VCF as a one-element array. (required)"
        prefix:              "Basename / sample-set name for all outputs. (required)"
        output_format:       "Format to emit: 'bcf' (default) or 'vcf' return the per-input filtered callsets as arrays (not combined); 'hmmibd_table' reduces each filtered file to an hmmIBD genotype table and STITCHES all of them into one combined table + freq file (the scalable path for large cohorts)."

        min_depth:           "Genotypes with FORMAT/DP below this value are set to missing. (default: 5)"
        keep_original_af:    "Preserve the original allele-frequency annotations (via rename_annots_tsv) before recomputing AN/AC/AF. (default: false)"
        variant_types:       "Which variant types to keep: 'snps', 'indels', or 'both'. Default keeps SNPs only; set 'both' to also run IBD on indels (hmmibd-rs codes by allele index, so bi-allelic indels are fine). (default: snps)"
        biallelic_only:      "Restrict to biallelic sites. (default: true)"
        split_multiallelics: "Split multiallelics into biallelic records to recover SNP alleles. (default: true)"
        max_variants:        "Optional cap on the number of variants; when >0 the callset is thinned evenly to at most this many. For output_format='hmmibd_table' the cap is applied to the COMBINED, genome-wide callset (per-file filtration runs uncapped); for 'bcf'/'vcf' it is applied per input file. (default: 0 = no limit)"
        max_alt:             "Optional pre-norm site width cap: when >0, sites with more than this many ALT alleles (after trimming unobserved ones) are dropped before `bcftools norm` splits them, bounding norm's memory on Pf hyper-multiallelic sites. ~6-8 is safe for Pf SNP IBD. (default: 0 = no cap)"
        gt_mode:             "output_format='hmmibd_table' only: how to derive each per-sample call. 'dominant-allele' (default) = max-depth allele from FORMAT/AD with hmmibd-rs read_dom gating (right for polyclonal Pf; requires AD to survive filtration); 'first-ploidy' = first allele of GT. (default: dominant-allele)"
        dom_min_depth:       "dominant-allele gate: minimum total AD depth (total > dom_min_depth), else the call is missing. Matches hmmibd-rs min_depth. (default: 5)"
        dom_min_ratio:       "dominant-allele gate: minimum major-allele fraction (major/total >= dom_min_ratio). Matches hmmibd-rs min_ratio. (default: 0.7)"
        dom_min_r1_r2:       "dominant-allele gate: accept only if minor/major < 1/dom_min_r1_r2. Matches hmmibd-rs min_r1_r2. (default: 3.0)"
        mask_gq0_genotypes:  "Also set GQ0 genotypes to missing (beyond the DP < min_depth mask). Reproduces the GATK GQ0-hom-ref fix (issue #7792 / PR #8741) for VCFs from pre-4.6.0.0 GenotypeGVCFs or GnarlyGenotyper, where no-/low-confidence hom-refs are 0/0:GQ=0 instead of ./. — the DP mask alone misses GQ0 calls with DP >= min_depth. Conservative: drops all GQ0 hom-refs (use a GVCF cross-reference to keep well-covered ones). Masks GT only: affects the first-ploidy table (and AC/AF/site selection); the dominant-allele table ignores GT (reads AD), so it is a no-op there — the dom_* AD gates handle GQ0 hom-refs instead. (default: false)"
        populations_file:    "Optional sample-to-population file for per-population AN/AC/AF. (default: none)"
        rename_annots_tsv:   "Required when keep_original_af=true: old-name<TAB>new-name TSV. (default: none)"
        filter_extra_args:   "Additional args appended verbatim to the final bcftools view invocation of each per-file filtration. (default: empty)"

        filter_runtime_attr_override:  "Override runtime for the (scattered) filtration task. (default: none)"
        convert_runtime_attr_override: "Override runtime for the vcf/table conversion task. (default: none)"
        combine_runtime_attr_override: "Override runtime for the table-stitching task. (default: none)"
    }

    input {
        Array[File] input_vcfs
        String prefix
        String output_format = "bcf"

        Int min_depth = 5
        Boolean keep_original_af = false
        String variant_types = "snps"
        Boolean biallelic_only = true
        Boolean split_multiallelics = true
        Int max_variants = 0
        Int max_alt = 0
        Boolean mask_gq0_genotypes = false
        File? populations_file
        File? rename_annots_tsv
        String filter_extra_args = ""

        String gt_mode = "dominant-allele"
        Int dom_min_depth = 5
        Float dom_min_ratio = 0.7
        Float dom_min_r1_r2 = 3.0

        RuntimeAttr? filter_runtime_attr_override
        RuntimeAttr? convert_runtime_attr_override
        RuntimeAttr? combine_runtime_attr_override
    }

    # For the combined-table path the per-file filtration runs UNCAPPED (max_variants=0) and the cap
    # is applied once to the stitched, genome-wide callset. For the array (bcf/vcf) paths there is no
    # combine step, so the cap is applied per input file.
    Int per_file_max_variants = if output_format == "hmmibd_table" then 0 else max_variants

    # Filter each input VCF independently. Filtration always produces a BCF.
    scatter (idx in range(length(input_vcfs))) {
        String shard_prefix = "~{prefix}.shard-~{idx}"

        call HMMIBD.FilterVcfForHmmIBD as t_01_Filter {
            input:
                input_vcf           = input_vcfs[idx],
                prefix              = shard_prefix,
                min_depth           = min_depth,
                keep_original_af    = keep_original_af,
                variant_types       = variant_types,
                biallelic_only      = biallelic_only,
                split_multiallelics = split_multiallelics,
                max_variants        = per_file_max_variants,
                max_alt             = max_alt,
                mask_gq0_genotypes  = mask_gq0_genotypes,
                populations_file    = populations_file,
                rename_annots_tsv   = rename_annots_tsv,
                extra_args          = filter_extra_args,
                runtime_attr_override = filter_runtime_attr_override
        }

        # For the array outputs, convert each filtered BCF per input file.
        if (output_format == "vcf") {
            call HMMIBD.BcfToVcf as t_02_ToVcf {
                input:
                    input_bcf = t_01_Filter.filtered_bcf,
                    prefix    = shard_prefix,
                    runtime_attr_override = convert_runtime_attr_override
            }
        }

        # For the combined-table path, reduce each filtered BCF to its compact hmmIBD table first;
        # only these small tables are stitched below (a merged BCF is never built).
        if (output_format == "hmmibd_table") {
            call HMMIBD.BcfToSampleTable as t_02_ToTable {
                input:
                    input_bcf     = t_01_Filter.filtered_bcf,
                    prefix        = shard_prefix,
                    gt_mode       = gt_mode,
                    dom_min_depth = dom_min_depth,
                    dom_min_ratio = dom_min_ratio,
                    dom_min_r1_r2 = dom_min_r1_r2,
                    runtime_attr_override = convert_runtime_attr_override
            }
        }
    }

    # Stitch the per-file genotype tables into ONE combined table + freq file (applies the cap
    # to the combined callset). select_all drops the None entries from the scattered optionals.
    if (output_format == "hmmibd_table") {
        call HMMIBD.StitchHmmIBDTables as t_03_Stitch {
            input:
                gt_tables    = select_all(t_02_ToTable.sample_gt_table),
                prefix       = prefix,
                max_variants = max_variants,
                runtime_attr_override = combine_runtime_attr_override
        }
    }

    output {
        Array[File]  filtered_bcfs        = t_01_Filter.filtered_bcf
        Array[File]  filtered_bcf_indices = t_01_Filter.filtered_bcf_index
        Array[File?] filtered_vcfs        = t_02_ToVcf.filtered_vcf
        Array[File?] filtered_vcf_indices = t_02_ToVcf.filtered_vcf_index
        File? sample_gt_table   = t_03_Stitch.sample_gt_table
        File? sample_freq_table = t_03_Stitch.sample_freq_table
    }
}
