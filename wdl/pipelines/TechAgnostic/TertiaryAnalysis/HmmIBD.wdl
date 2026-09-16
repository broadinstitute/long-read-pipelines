version 1.0

import "../../../structs/Structs.wdl"
import "../../../tasks/TertiaryAnalysis/HmmIBDTasks.wdl" as HMMIBD

workflow HmmIBD {

    meta {
        author: "Jonn Smith"
        description: "Infer identity-by-descent (IBD) from a VCF. Step one pre-filters the VCF with bcftools (mask low-depth genotypes, recompute AN/AC/AF, trim unrepresented alleles/sites) into a BCF; step two runs hmmibd-rs (the Rust reimplementation of hmmIBD) on that BCF. The hmmibd-rs IBD segments and per-pair IBD fractions are returned as workflow outputs. Every hmmibd-rs tuning parameter is exposed with its upstream default so runs are fully configurable from Terra."

        outputs: {
            ibd_segments:     "Per sample-pair IBD/non-IBD segments (<prefix>.hmm.txt); empty when the HMM is skipped via a bcf-to-bin mode",
            ibd_fraction:     "Per sample-pair IBD-fraction summary (<prefix>.hmm_fract.txt); empty when suppress_frac=true or a bcf-to-bin mode is used",
            binary_genotypes: "Binary genotype file(s) (<prefix>*.bin); populated only when bcf_to_bin_file or bcf_to_bin_file_by_chromosome is set",
            filtered_bcf:     "The pre-filtered BCF that hmmibd-rs was run on"
        }
    }

    parameter_meta {
        input_vcf:        "VCF/BCF to infer IBD from. (required)"
        prefix:           "Basename / sample-set name for all outputs. (required)"

        # ---- Step 1: filtration (bcftools) ----
        min_depth:         "Genotypes with FORMAT/DP below this value are set to missing. (default: 5)"
        keep_original_af:  "Preserve the original allele-frequency annotations (via rename_annots_tsv) before recomputing AN/AC/AF. (default: false)"
        variant_types:     "Which variant types to keep: 'snps', 'indels', or 'both'. SNPs are the standard hmmIBD marker set. (default: snps)"
        biallelic_only:    "Restrict to biallelic sites. Recommended true: multiallelic sites (especially indels) can exceed hmmibd-rs --max-all and crash it until patched. (default: true)"
        max_variants:      "Optional cap on the number of variants; when >0, the filtered callset is thinned evenly across the genome to at most this many. (default: 0 = no limit)"
        populations_file:  "Optional sample-to-population file for per-population AN/AC/AF (bcftools +fill-tags -S). (default: none)"
        rename_annots_tsv: "Required when keep_original_af=true: old-name<TAB>new-name TSV for bcftools annotate --rename-annots. (default: none)"

        # ---- Step 2: hmmibd-rs ----
        data_file2:           "Optional second-population genotypes (-I). (default: none)"
        freq_file1:           "Optional allele-frequency file for population 1 (-f). (default: none)"
        freq_file2:           "Optional allele-frequency file for population 2 (-F). (default: none)"
        bad_samples_file:     "Optional list of sample IDs to exclude (-b). (default: none)"
        good_pairs_file:      "Optional list of sample pairs to analyze (-g). (default: none)"
        bcf_filter_config:    "Optional TOML BCF filter config (--bcf-filter-config). (default: none)"
        genome:               "Optional genome/recombination-map spec (--genome); mutually exclusive with rec_rate. (default: none)"
        bcf_read_mode:        "BCF genotype read mode: dominant-allele | first-ploidy | each-ploidy. (default: dominant-allele)"
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
        bcf_to_bin_file:               "Convert the BCF to a single binary genotype file and skip the HMM (--bcf-to-bin-file). (default: false)"
        bcf_to_bin_file_by_chromosome: "Convert the BCF to one binary genotype file per chromosome and skip the HMM (--bcf-to-bin-file-by-chromosome). (default: false)"
    }

    input {
        File input_vcf
        String prefix

        # ---- Step 1: filtration ----
        Int min_depth = 5
        Boolean keep_original_af = false
        String variant_types = "snps"
        Boolean biallelic_only = true
        Int max_variants = 0
        File? populations_file
        File? rename_annots_tsv

        # ---- Step 2: hmmibd-rs ----
        File? data_file2
        File? freq_file1
        File? freq_file2
        File? bad_samples_file
        File? good_pairs_file
        File? bcf_filter_config
        File? genome

        String bcf_read_mode = "dominant-allele"

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

    # Step 1: filter the VCF down to the features/samples of interest and recompute allele stats.
    call HMMIBD.FilterVcfForHmmIBD as t_01_FilterVcfForHmmIBD {
        input:
            input_vcf         = input_vcf,
            prefix            = prefix,
            min_depth           = min_depth,
            keep_original_af    = keep_original_af,
            variant_types       = variant_types,
            biallelic_only      = biallelic_only,
            max_variants        = max_variants,
            populations_file    = populations_file,
            rename_annots_tsv   = rename_annots_tsv
    }

    # Step 2: infer IBD on the filtered BCF.
    call HMMIBD.HmmIBDrs as t_02_HmmIBDrs {
        input:
            input_bcf            = t_01_FilterVcfForHmmIBD.filtered_bcf,
            prefix               = prefix,
            data_file2           = data_file2,
            freq_file1           = freq_file1,
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
        File filtered_bcf            = t_01_FilterVcfForHmmIBD.filtered_bcf
    }
}
