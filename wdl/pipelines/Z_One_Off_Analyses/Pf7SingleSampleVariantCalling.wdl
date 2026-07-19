version 1.0

## Pf7SingleSampleVariantCalling.wdl
## Per-sample variant calling following the MalariaGEN Pf7 methods
## (Wellcome Open Res 2023, doi:10.12688/wellcomeopenres.18681.1, pp.1-2).
## STARTS FROM AN ANALYSIS-READY BAM — alignment (bwa mem 0.7.15 -M),
## samtools fixmate, Picard MarkDuplicates, and GATK BQSR are all assumed done
## upstream. This workflow runs, for one sample:
##   HaplotypeCaller 4.1.4.0  -contamination 0 -ERC GVCF   -> per-sample gVCF
##   GenotypeGVCFs 4.1.4.0    (Pf7 flags, whole-genome)    -> genotyped VCF
##   VQSR (SNP & INDEL separate; train = Pf crosses 1.0 PASS, prior 15;
##         -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum, --maxGaussians 8)
##   VQSLOD<2.0 filter in core, keep non-core, merge              -> filtered VCF
##
## NOTE ON SINGLE-SAMPLE VQSR: VariantRecalibrator is built for cohort-scale data
## and may fail ("no data / zero variance in annotations") on a single sample.
## The Pf7 headline callset is a JOINT-genotyped multisample VCF — for that, run
## Pf7JointGenotyping.wdl over the full cohort. This workflow produces the
## single-sample VCF for the sample-level QC / comparison case.
##
## Each task uses only tools in its image: GATK tasks in the gatk image,
## VCF-wrangling (concat/filter/merge) in the bcftools image.

workflow Pf7SingleSampleVariantCalling {
  input {
    String sample_name
    # Analysis-ready reads: aligned, dup-marked, BQSR-recalibrated BAM (+ index).
    File analysis_ready_bam
    File analysis_ready_bai

    # Reference: PlasmoDB-61 3D7 v3 (defaults = broad-malaria-public / sr-malaria workspace)
    File ref_fasta     = "gs://broad-malaria-public/short_read_workspace_data/reference/PlasmoDB-61_Pfalciparum3D7_Genome.fasta"
    File ref_fasta_fai = "gs://broad-malaria-public/short_read_workspace_data/reference/PlasmoDB-61_Pfalciparum3D7_Genome.fasta.fai"
    File ref_dict      = "gs://broad-malaria-public/short_read_workspace_data/reference/PlasmoDB-61_Pfalciparum3D7_Genome.dict"

    # Pf crosses 1.0 PASS (VQSR training) + core-region bed (defaults = workspace constants)
    File pfcrosses_pass_vcf       = "gs://broad-malaria-public/short_read_workspace_data/ALL_CROSSES.sites_only.PASS.vcf.gz"
    File pfcrosses_pass_vcf_index = "gs://broad-malaria-public/short_read_workspace_data/ALL_CROSSES.sites_only.PASS.vcf.gz.tbi"
    File core_bed                 = "gs://broad-malaria-public/short_read_workspace_data/regions/regions-20130225.Core.bed"

    Float vqsr_prior         = 15.0
    Int   vqsr_max_gaussians = 8
    Float vqslod_threshold   = 2.0

    String gatk_docker     = "broadinstitute/gatk:4.1.4.0"
    String bcftools_docker = "quay.io/biocontainers/bcftools:1.13--h3a49de5_0"
  }

  # 1) analysis-ready BAM -> per-sample gVCF
  call HaplotypeCallerGVCF {
    input: sample_name = sample_name,
           input_bam = analysis_ready_bam, input_bai = analysis_ready_bai,
           ref_fasta = ref_fasta, ref_fasta_fai = ref_fasta_fai, ref_dict = ref_dict,
           gatk_docker = gatk_docker
  }

  # 2) gVCF -> genotyped VCF (single-sample GenotypeGVCFs, whole genome, Pf7 flags)
  call GenotypeSample {
    input: sample_name = sample_name,
           gvcf = HaplotypeCallerGVCF.gvcf, gvcf_index = HaplotypeCallerGVCF.gvcf_index,
           ref_fasta = ref_fasta, ref_fasta_fai = ref_fasta_fai, ref_dict = ref_dict,
           gatk_docker = gatk_docker
  }

  # 3) VQSR annotate (SNP & indel separate) -> VQSLOD in INFO
  call VQSRAnnotate {
    input: cohort_vcf = GenotypeSample.vcf, cohort_vcf_index = GenotypeSample.vcf_index,
           ref_fasta = ref_fasta, ref_fasta_fai = ref_fasta_fai, ref_dict = ref_dict,
           pfcrosses_pass_vcf = pfcrosses_pass_vcf, pfcrosses_pass_vcf_index = pfcrosses_pass_vcf_index,
           vqsr_prior = vqsr_prior, vqsr_max_gaussians = vqsr_max_gaussians, gatk_docker = gatk_docker
  }

  # 4) VQSLOD<threshold in core, keep non-core, merge -> final single-sample VCF
  call FilterMerge {
    input: cohort_name = sample_name,
           snp_vqslod_vcf = VQSRAnnotate.snp_vqslod_vcf, snp_vqslod_vcf_index = VQSRAnnotate.snp_vqslod_vcf_index,
           indel_vqslod_vcf = VQSRAnnotate.indel_vqslod_vcf, indel_vqslod_vcf_index = VQSRAnnotate.indel_vqslod_vcf_index,
           core_bed = core_bed, vqslod_threshold = vqslod_threshold, bcftools_docker = bcftools_docker
  }

  output {
    File gvcf                 = HaplotypeCallerGVCF.gvcf
    File gvcf_index           = HaplotypeCallerGVCF.gvcf_index
    File genotyped_vcf        = GenotypeSample.vcf
    File genotyped_vcf_index  = GenotypeSample.vcf_index
    File filtered_vcf         = FilterMerge.filtered_vcf
    File filtered_vcf_index   = FilterMerge.filtered_vcf_index
  }
}

# GATK 4.1.4.0 HaplotypeCaller -ERC GVCF (Pf7 defaults: ploidy 2, het 1e-3, indel-het 1.25e-4)
task HaplotypeCallerGVCF {
  input {
    String sample_name
    File input_bam
    File input_bai
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    String gatk_docker
    Int memory_gb = 12
    Int disk_gb = 60
  }
  command <<<
    set -euo pipefail
    gatk --java-options "-Xmx~{memory_gb - 2}g" HaplotypeCaller \
      -R ~{ref_fasta} -I ~{input_bam} \
      -contamination 0 -ERC GVCF \
      -O ~{sample_name}.g.vcf.gz
  >>>
  runtime {
    docker: gatk_docker
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
  }
  output {
    File gvcf = "~{sample_name}.g.vcf.gz"
    File gvcf_index = "~{sample_name}.g.vcf.gz.tbi"
  }
}

# Single-sample GenotypeGVCFs (whole genome). No interval scatter / padding needed
# for one sample — Pf7's 10kb tiling is a joint-calling parallelization; for a lone
# gVCF the whole-genome call is identical. Same Pf7 genotyping annotation flags.
task GenotypeSample {
  input {
    String sample_name
    File gvcf
    File gvcf_index
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    String gatk_docker
    Int memory_gb = 12
    Int disk_gb = 50
  }
  command <<<
    set -euo pipefail
    gatk --java-options "-Xmx~{memory_gb - 2}g" GenotypeGVCFs \
      -R ~{ref_fasta} -V ~{gvcf} \
      --use-new-qual-calculator \
      --annotation-group StandardAnnotation --annotation-group AS_StandardAnnotation \
      -O ~{sample_name}.genotyped.vcf.gz
  >>>
  runtime {
    docker: gatk_docker
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
  }
  output {
    File vcf = "~{sample_name}.genotyped.vcf.gz"
    File vcf_index = "~{sample_name}.genotyped.vcf.gz.tbi"
  }
}

# GATK VQSR: VariantRecalibrator + ApplyVQSR for SNP and INDEL separately.
# ApplyVQSR at ts 100 => annotate VQSLOD without tranche-filtering (thresholded later).
task VQSRAnnotate {
  input {
    File cohort_vcf
    File cohort_vcf_index
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    File pfcrosses_pass_vcf
    File pfcrosses_pass_vcf_index
    Float vqsr_prior
    Int   vqsr_max_gaussians
    String gatk_docker
    Int memory_gb = 16
    Int disk_gb = 50
  }
  command <<<
    set -euo pipefail
    for mode in SNP INDEL; do
      tag=$(echo "$mode" | tr '[:upper:]' '[:lower:]')
      gatk --java-options "-Xmx~{memory_gb - 4}g" VariantRecalibrator \
        -R ~{ref_fasta} -V ~{cohort_vcf} -mode "$mode" \
        --resource:pfcrosses,known=false,training=true,truth=true,prior=~{vqsr_prior} ~{pfcrosses_pass_vcf} \
        -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
        --trust-all-polymorphic --max-gaussians ~{vqsr_max_gaussians} \
        -O "$tag.recal" --tranches-file "$tag.tranches"
      gatk --java-options "-Xmx~{memory_gb - 4}g" ApplyVQSR \
        -R ~{ref_fasta} -V ~{cohort_vcf} -mode "$mode" \
        --recal-file "$tag.recal" --tranches-file "$tag.tranches" \
        --truth-sensitivity-filter-level 100.0 \
        -O "$tag.vqslod.vcf.gz"
    done
  >>>
  runtime {
    docker: gatk_docker
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
  }
  output {
    File snp_vqslod_vcf = "snp.vqslod.vcf.gz"
    File snp_vqslod_vcf_index = "snp.vqslod.vcf.gz.tbi"
    File indel_vqslod_vcf = "indel.vqslod.vcf.gz"
    File indel_vqslod_vcf_index = "indel.vqslod.vcf.gz.tbi"
  }
}

# Take SNPs from SNP-recal + indels from INDEL-recal; VQSLOD<threshold filtered in
# core; non-core kept unconditionally; merge -> final VCF.
task FilterMerge {
  input {
    String cohort_name
    File snp_vqslod_vcf
    File snp_vqslod_vcf_index
    File indel_vqslod_vcf
    File indel_vqslod_vcf_index
    File core_bed
    Float vqslod_threshold
    String bcftools_docker
    Int memory_gb = 8
    Int disk_gb = 50
  }
  command <<<
    set -euo pipefail
    bcftools view -v snps   -Oz -o snps.vcf.gz   ~{snp_vqslod_vcf};   bcftools index -t snps.vcf.gz
    bcftools view -v indels -Oz -o indels.vcf.gz ~{indel_vqslod_vcf}; bcftools index -t indels.vcf.gz
    bcftools concat -a snps.vcf.gz indels.vcf.gz | bcftools sort -Oz -o vqslod.vcf.gz
    bcftools index -t vqslod.vcf.gz

    # core: keep VQSLOD >= threshold ; non-core: keep unconditionally
    bcftools view -R ~{core_bed} vqslod.vcf.gz \
      | bcftools view -e "INFO/VQSLOD < ~{vqslod_threshold}" -Oz -o core.pass.vcf.gz
    bcftools index -t core.pass.vcf.gz
    bcftools view -T "^~{core_bed}" vqslod.vcf.gz -Oz -o noncore.keep.vcf.gz
    bcftools index -t noncore.keep.vcf.gz

    # Pf7 used GATK3.8 CombineVariants; bcftools concat here (documented divergence).
    bcftools concat -a core.pass.vcf.gz noncore.keep.vcf.gz \
      | bcftools sort -Oz -o ~{cohort_name}.filtered.vcf.gz
    bcftools index -t ~{cohort_name}.filtered.vcf.gz
  >>>
  runtime {
    docker: bcftools_docker
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
  }
  output {
    File filtered_vcf = "~{cohort_name}.filtered.vcf.gz"
    File filtered_vcf_index = "~{cohort_name}.filtered.vcf.gz.tbi"
  }
}
