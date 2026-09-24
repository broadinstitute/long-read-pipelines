version 1.0

import "../../../tasks/TertiaryAnalysis/HmmIBDTasks.wdl" as HMMIBD
import "FilterVcfForHmmIBD.wdl" as FILTER

workflow HmmIBD {

    meta {
        author: "Jonn Smith"
        description: "Infer identity-by-descent (IBD) from a LIST of VCFs with hmmibd-rs (the Rust reimplementation of hmmIBD), designed to scale to cohorts split across many large files (hundreds of GB). With run_filtration=true, each input VCF is filtered independently and reduced to a compact hmmIBD genotype table, and all tables are STITCHED into one combined, coordinate-sorted table + allele-frequency file — a merged multi-hundred-GB BCF is never materialized. hmmibd-rs then runs on that single combined table. The inputs are meant to be disjoint regions of ONE joint call (same samples across files), so the stitched frequencies are exact. Every hmmibd-rs tuning parameter is exposed with its upstream default."

        outputs: {
            ibd_segments:       "Per sample-pair IBD/non-IBD segments (<prefix>.hmm.txt); empty when the HMM is skipped via a bcf-to-bin mode",
            ibd_fraction:       "Per sample-pair IBD-fraction summary (<prefix>.hmm_fract.txt); empty when suppress_frac=true or a bcf-to-bin mode is used",
            binary_genotypes:   "Binary genotype file(s) (<prefix>*.bin); populated only when bcf_to_bin_file or bcf_to_bin_file_by_chromosome is set",
            filtered_bcfs:      "The per-input filtered BCFs (one per input VCF); populated when run_filtration=true",
            combined_gt_table:  "The single stitched hmmIBD genotype table hmmibd-rs was run on (when the combined-table path is used)",
            combined_freq_table:"The matching combined allele-frequency file (when the combined-table path is used)"
        }
    }

    parameter_meta {
        input_vcfs:       "Genotype inputs, one or more. Meant to be disjoint regions of ONE joint call (same samples across files); pass a single file as a one-element array. With run_filtration=true each is filtered + reduced to a table, then all are stitched into one combined table for hmmibd-rs. With run_filtration=false: if hmmibd_input_format='hmmibd_table' the elements are treated as pre-built per-region hmmIBD tables and STITCHED (re-run IBD without re-filtering hundreds of GB); if 'bcf'/'vcf' the FIRST element is fed straight to hmmibd-rs via --from-bcf (a list is NOT combined in this mode). (required)"
        prefix:           "Basename / sample-set name for all outputs. (required)"

        run_filtration:      "Filter + reduce + stitch input_vcfs via the FilterVcfForHmmIBD workflow before IBD. Set false to feed a prepared input straight to hmmibd-rs (see input_vcfs / hmmibd_input_format). (default: true)"
        hmmibd_input_format: "Only used when run_filtration=false: format of input_vcfs — 'hmmibd_table' (per-region hmmIBD text tables, stitched into one) or 'bcf'/'vcf' (first element read via --from-bcf). (default: bcf)"

        # ---- Step 1: filtration (bcftools) ----
        min_depth:         "Genotypes with FORMAT/DP below this value are set to missing. (default: 5)"
        keep_original_af:  "Preserve the original allele-frequency annotations (via rename_annots_tsv) before recomputing AN/AC/AF. (default: false)"
        variant_types:     "Which variant types to keep: 'snps', 'indels', or 'both'. Default keeps SNPs only; set 'both' to also run IBD on indels (hmmibd-rs codes by allele index, so bi-allelic indels are fine). (default: snps)"
        biallelic_only:    "Restrict to biallelic sites. Recommended true: multiallelic sites (especially indels) can exceed hmmibd-rs --max-all and crash it until patched. (default: true)"
        split_multiallelics: "Split multiallelics into biallelic records (bcftools norm -m-any) to recover SNP alleles from multiallelic/spanning-deletion sites instead of dropping them. (default: true)"
        max_variants:      "Optional cap on the number of variants; when >0, sites are thinned evenly across the genome to at most this many. Applied to the COMBINED, genome-wide callset (per-file filtration runs uncapped). (default: 0 = no limit)"
        max_alt:           "Optional pre-norm site width cap: when >0, sites with more than this many ALT alleles (after trimming unobserved ones) are dropped before `bcftools norm` splits them — bounds the filtration task's memory on Pf hyper-multiallelic (var-gene/indel) sites that OOM-kill norm on whole-cohort joint calls. ~6-8 is safe for Pf SNP IBD. (default: 0 = no cap)"
        mask_gq0_genotypes: "Also set GQ0 genotypes to missing (beyond the DP < min_depth mask). Reproduces the GATK GQ0-hom-ref fix (issue #7792 / PR #8741) for VCFs from pre-4.6.0.0 GenotypeGVCFs or GnarlyGenotyper, where no-/low-confidence hom-refs are 0/0:GQ=0 instead of ./. Conservative (drops all GQ0 hom-refs). Only applied when run_filtration=true. Masks GT only, so with bcf_read_mode='dominant-allele' (the default) it is a no-op for the table's call values (that path reads AD, ignores GT) — the dom_min_depth/dom_min_ratio/dom_min_r1_r2 AD gates handle GQ0 hom-refs there; it bites in first-ploidy mode. (default: false)"
        populations_file:  "Optional sample-to-population file for per-population AN/AC/AF (bcftools +fill-tags -S). (default: none)"
        rename_annots_tsv: "Required when keep_original_af=true: old-name<TAB>new-name TSV for bcftools annotate --rename-annots. (default: none)"
        filter_extra_args: "Extra `bcftools view` filters added to the final per-file filtration view (e.g. \"-e MAF<0.01\"); whitespace-separated, no internal spaces. Only applied when run_filtration=true. (default: empty)"

        # ---- Step 2: hmmibd-rs ----
        data_file2:           "Optional second-population genotypes (-I). (default: none)"
        freq_file1:           "Optional allele-frequency file for population 1 (-f); overrides the stitched combined frequencies when supplied. (default: none; computed from the combined table)"
        freq_file2:           "Optional allele-frequency file for population 2 (-F). (default: none)"
        bad_samples_file:     "Optional list of sample IDs to exclude (-b). (default: none)"
        good_pairs_file:      "Optional list of sample pairs to analyze (-g). (default: none)"
        bcf_filter_config:    "Optional TOML BCF filter config (--bcf-filter-config); only used in --from-bcf mode. (default: none)"
        genome:               "Optional genome/recombination-map spec (--genome); mutually exclusive with rec_rate. (default: none)"
        bcf_read_mode:        "How genotypes are read. In --from-bcf mode: dominant-allele | first-ploidy | each-ploidy. When run_filtration=true this ALSO governs how the combined table is built (dominant-allele | first-ploidy; each-ploidy is --from-bcf only) — so dominant-allele gives the same majority-clone calls whether via table or --from-bcf. (default: dominant-allele)"
        dom_min_depth:        "dominant-allele table gate: minimum total AD depth (total > dom_min_depth), else missing. Matches hmmibd-rs min_depth. (default: 5)"
        dom_min_ratio:        "dominant-allele table gate: minimum major-allele fraction (major/total >= dom_min_ratio). Matches hmmibd-rs min_ratio. (default: 0.7)"
        dom_min_r1_r2:        "dominant-allele table gate: accept only if minor/major < 1/dom_min_r1_r2. Matches hmmibd-rs min_r1_r2. (default: 3.0)"
        max_iter:             "Max EM iterations (-m). (default: 5)"
        k_rec_max:            "Cap on inferred generations (-n). (default: none; hmmibd-rs uses Inf = no cap)"
        eps:                  "Genotyping error rate (--eps). (default: 0.001)"
        min_inform:           "Minimum informative sites per pair (--min-inform). (default: 10)"
        min_discord:          "Minimum discordance fraction (--min-discord). (default: 0.0)"
        max_discord:          "Maximum discordance fraction (--max-discord). (default: 1.0)"
        min_snp_sep:          "Minimum bp between SNPs (--min-snp-sep). (default: 5)"
        fit_thresh_dpi:       "Convergence threshold on pi (--fit-thresh-dpi). (default: 0.001)"
        fit_thresh_dk:        "Convergence threshold on k (--fit-thresh-dk). (default: 0.01)"
        fit_thresh_drelk:     "Convergence threshold on relative k (--fit-thresh-drelk). (default: 0.001)"
        rec_rate:             "Constant recombination rate per generation per bp (-r); ignored when genome is set. (default: 0.00000074)"
        max_all:              "Maximum unique alleles per site (--max-all). (default: 8)"
        buffer_size_segments: "Segments-output write buffer in bytes (--buffer-size-segments). (default: none; hmmibd-rs uses 8Kb)"
        buffer_size_frac:     "Fraction-output write buffer in bytes (--buffer-size-frac). (default: none; hmmibd-rs uses 8Kb)"
        filt_min_seg_cm:      "Drop output segments shorter than this many cM (--filt-min-seg-cm). (default: none)"
        filt_max_tmrca:       "Drop pairs with k_rec above this threshold (--filt-max-tmrca). (default: none)"
        filt_ibd_only:        "Omit non-IBD segments (--filt-ibd-only). (default: false)"
        num_threads:          "Threads; 0 = all CPUs (--num-threads). (default: 0)"
        par_mode:             "Parallelization mode: 0 small / 1 large sample sets (--par-mode). (default: 0)"
        par_chunk_size:       "Pairs-per-chunk (mode 0) or samples-per-chunk (mode 1) (--par-chunk-size). (default: 120)"

        suppress_frac:                 "Suppress the per-pair fraction output (--suppress-frac); leaves ibd_fraction empty. (default: false)"
        bcf_to_bin_file:               "Convert to a single binary genotype file and skip the HMM (--bcf-to-bin-file). (default: false)"
        bcf_to_bin_file_by_chromosome: "Convert to one binary genotype file per chromosome and skip the HMM (--bcf-to-bin-file-by-chromosome). (default: false)"
    }

    input {
        Array[File] input_vcfs
        String prefix

        Boolean run_filtration = true
        String hmmibd_input_format = "bcf"

        # ---- Step 1: filtration ----
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

        # ---- Step 2: hmmibd-rs ----
        File? data_file2
        File? freq_file1
        File? freq_file2
        File? bad_samples_file
        File? good_pairs_file
        File? bcf_filter_config
        File? genome

        String bcf_read_mode = "dominant-allele"
        Int dom_min_depth = 5
        Float dom_min_ratio = 0.7
        Float dom_min_r1_r2 = 3.0

        Int max_iter = 5
        Float? k_rec_max
        Float eps = 0.001
        Int min_inform = 10
        Float min_discord = 0.0
        Float max_discord = 1.0
        Int min_snp_sep = 5
        Float fit_thresh_dpi = 0.001
        Float fit_thresh_dk = 0.01
        Float fit_thresh_drelk = 0.001
        String rec_rate = "0.00000074"

        Int max_all = 8
        Int? buffer_size_segments
        Int? buffer_size_frac

        Float? filt_min_seg_cm
        Float? filt_max_tmrca
        Boolean filt_ibd_only = false

        Int num_threads = 0
        Int par_mode = 0
        Int par_chunk_size = 120

        Boolean suppress_frac = false
        Boolean bcf_to_bin_file = false
        Boolean bcf_to_bin_file_by_chromosome = false
    }

    # Whether the genotype input for hmmibd-rs is a text table (-i) rather than a BCF/VCF (--from-bcf):
    # true when we filter (always emits + stitches a table) or when given pre-built per-region tables.
    Boolean use_table = run_filtration || (hmmibd_input_format == "hmmibd_table")

    # Step 1a (run_filtration): filter each input VCF, reduce to a table, and stitch into ONE
    # combined table + freq file. The FilterVcfForHmmIBD workflow does the scatter + stitch.
    if (run_filtration) {
        call FILTER.FilterVcfForHmmIBD as t_01_FilterAndCombine {
            input:
                input_vcfs          = input_vcfs,
                prefix              = prefix,
                output_format       = "hmmibd_table",
                min_depth           = min_depth,
                keep_original_af    = keep_original_af,
                variant_types       = variant_types,
                biallelic_only      = biallelic_only,
                split_multiallelics = split_multiallelics,
                max_variants        = max_variants,
                max_alt             = max_alt,
                mask_gq0_genotypes  = mask_gq0_genotypes,
                populations_file    = populations_file,
                rename_annots_tsv   = rename_annots_tsv,
                filter_extra_args   = filter_extra_args,
                gt_mode             = bcf_read_mode,
                dom_min_depth       = dom_min_depth,
                dom_min_ratio       = dom_min_ratio,
                dom_min_r1_r2       = dom_min_r1_r2
        }
    }

    # Step 1b (run_filtration=false + pre-built tables): stitch the provided per-region hmmIBD
    # tables into one combined table, so IBD can be re-run without re-filtering the raw VCFs.
    if (!run_filtration && hmmibd_input_format == "hmmibd_table") {
        call HMMIBD.StitchHmmIBDTables as t_01b_Stitch {
            input:
                gt_tables    = input_vcfs,
                prefix       = prefix,
                max_variants = max_variants
        }
    }

    # Resolve the hmmibd-rs genotype input:
    #  - combined table (filtered, or stitched from pre-built tables) -> text mode (-i)
    #  - otherwise (run_filtration=false + bcf/vcf) -> the first prepared file via --from-bcf
    File hmmibd_input    = select_first([t_01_FilterAndCombine.sample_gt_table, t_01b_Stitch.sample_gt_table, input_vcfs[0]])
    Boolean use_from_bcf = !use_table

    # Use the caller-supplied freq file if given, else the stitched combined frequencies (text
    # mode), else none (bcf/vcf mode computes frequencies internally).
    File? resolved_freq_file1 = if defined(freq_file1) then freq_file1
                                else if defined(t_01_FilterAndCombine.sample_freq_table) then t_01_FilterAndCombine.sample_freq_table
                                else t_01b_Stitch.sample_freq_table

    # Step 2: infer IBD on the resolved input.
    call HMMIBD.HmmIBDrs as t_02_HmmIBDrs {
        input:
            input_bcf            = hmmibd_input,
            from_bcf             = use_from_bcf,
            prefix               = prefix,
            data_file2           = data_file2,
            freq_file1           = resolved_freq_file1,
            freq_file2           = freq_file2,
            bad_samples_file     = bad_samples_file,
            good_pairs_file      = good_pairs_file,
            bcf_filter_config    = bcf_filter_config,
            genome               = genome,
            bcf_read_mode        = bcf_read_mode,
            max_iter             = max_iter,
            k_rec_max            = k_rec_max,
            eps                  = eps,
            min_inform           = min_inform,
            min_discord          = min_discord,
            max_discord          = max_discord,
            min_snp_sep          = min_snp_sep,
            fit_thresh_dpi       = fit_thresh_dpi,
            fit_thresh_dk        = fit_thresh_dk,
            fit_thresh_drelk     = fit_thresh_drelk,
            rec_rate             = rec_rate,
            max_all              = max_all,
            buffer_size_segments = buffer_size_segments,
            buffer_size_frac     = buffer_size_frac,
            filt_min_seg_cm      = filt_min_seg_cm,
            filt_max_tmrca       = filt_max_tmrca,
            filt_ibd_only        = filt_ibd_only,
            num_threads          = num_threads,
            par_mode             = par_mode,
            par_chunk_size       = par_chunk_size,
            suppress_frac                 = suppress_frac,
            bcf_to_bin_file               = bcf_to_bin_file,
            bcf_to_bin_file_by_chromosome = bcf_to_bin_file_by_chromosome
    }

    output {
        # Which of these are populated depends on the hmmibd-rs output-mode flags:
        # normal HMM run -> ibd_segments + ibd_fraction; --suppress-frac -> ibd_segments only;
        # --bcf-to-bin-file[-by-chromosome] -> binary_genotypes only (HMM skipped).
        Array[File] ibd_segments     = t_02_HmmIBDrs.ibd_segments
        Array[File] ibd_fraction     = t_02_HmmIBDrs.ibd_fraction
        Array[File] binary_genotypes = t_02_HmmIBDrs.binary_genotypes

        Array[File]? filtered_bcfs        = t_01_FilterAndCombine.filtered_bcfs
        File?        combined_gt_table    = if defined(t_01_FilterAndCombine.sample_gt_table) then t_01_FilterAndCombine.sample_gt_table else t_01b_Stitch.sample_gt_table
        File?        combined_freq_table  = resolved_freq_file1
    }
}
