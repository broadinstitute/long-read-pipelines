version 1.0

struct HumanReferenceBundle {
    File fasta
    File fai
    File dict

    File tandem_repeat_bed
    File PAR_bed

    File? size_balanced_scatter_interval_ids
    File? size_balanced_scatter_intervallists_locators

    File intervallists_autosomes
    File intervallists_allosomes

    File chromosome_ploidy_priors

    String mt_chr_name

    File? haplotype_map
}
