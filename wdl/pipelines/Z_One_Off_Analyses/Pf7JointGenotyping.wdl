version 1.0

## Pf7JointGenotyping.wdl
## Joint genotyping following the MalariaGEN Pf7 methods (pp.1-2):
##   genome tiled in 10 kbp intervals (2,342 across 3D7 v3)
##   per interval: GenomicsDBImport (-ip 500, --sample-name-map)
##                 -> GenotypeGVCFs 4.1.4.0 (--only-output-calls-starting-in-intervals
##                    --use-new-qual-calculator -L <int + 500bp before>
##                    --annotation-group StandardAnnotation --annotation-group AS_StandardAnnotation)
##                 -> excise 500bp padding (SelectVariants -L <pure window>)
##   gather -> VQSR (SNP & indel separate; train = Pf crosses 1.0 PASS, prior 15;
##             -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum, no DP/MQ;
##             --maxGaussians 8) -> VQSLOD<2.0 filter in core -> merge core/non-core.
## Ingests the gVCFs produced by Pf7SingleSampleVariantCalling.wdl.
##
## SCATTER WIDTH: the 2,342 intervals are BINNED into at most `max_shards`
## contiguous shards (default 10). Each shard task processes its intervals
## sequentially, so the fan-out is <=10 VMs, not 2,342. Contiguous binning keeps
## the intervals genome-ordered, so the final concat needs no re-sort.
##
## Each task uses only tools present in its image: GATK tasks in the gatk image,
## VCF-wrangling (concat/filter/merge) in the bcftools image. Pf7 excised padding
## with `bcftools view -r`; we use `gatk SelectVariants -L` (equivalent) to keep
## the shard task in the gatk image. Pf7's final core/non-core merge used GATK3.8
## CombineVariants; we use bcftools concat (documented divergence).

workflow Pf7JointGenotyping {
  input {
    String cohort_name = "pf7_cohort"
    # Sample names are DERIVED from the gVCF headers (bcftools query -l) — no
    # hand-maintained name array (mirrors long-read-pipelines CreateSampleNameMap).
    Array[File]   gvcfs
    Array[File]   gvcf_indices

    # Reference: PlasmoDB-61 3D7 v3 (defaults = broad-malaria-public / sr-malaria workspace)
    File ref_fasta     = "gs://broad-malaria-public/short_read_workspace_data/reference/PlasmoDB-61_Pfalciparum3D7_Genome.fasta"
    File ref_fasta_fai = "gs://broad-malaria-public/short_read_workspace_data/reference/PlasmoDB-61_Pfalciparum3D7_Genome.fasta.fai"
    File ref_dict      = "gs://broad-malaria-public/short_read_workspace_data/reference/PlasmoDB-61_Pfalciparum3D7_Genome.dict"

    # Pf crosses 1.0 PASS (VQSR training) + core-region bed (defaults = workspace constants)
    File pfcrosses_pass_vcf       = "gs://broad-malaria-public/short_read_workspace_data/ALL_CROSSES.sites_only.PASS.vcf.gz"
    File pfcrosses_pass_vcf_index = "gs://broad-malaria-public/short_read_workspace_data/ALL_CROSSES.sites_only.PASS.vcf.gz.tbi"
    File core_bed                 = "gs://broad-malaria-public/short_read_workspace_data/regions/regions-20130225.Core.bed"

    Int    interval_size   = 10000
    Int    interval_pad    = 500
    Int    max_shards      = 10
    Float  vqsr_prior      = 15.0
    Int    vqsr_max_gaussians = 8
    Float  vqslod_threshold = 2.0

    String gatk_docker = "broadinstitute/gatk:4.1.4.0"
    String bcftools_docker = "quay.io/biocontainers/bcftools:1.13--h3a49de5_0"
  }

  # 1) tile genome into interval_size windows, then bin into <=max_shards shard TSVs
  call MakeIntervals {
    input: ref_fasta_fai = ref_fasta_fai, interval_size = interval_size,
           max_shards = max_shards, bcftools_docker = bcftools_docker
  }

  # 2) derive the sample-name-map ONCE from the gVCF headers (streams the gVCFs,
  #    no localization). name<TAB>gs://path — consumed by every shard's import.
  call CreateSampleNameMap {
    input: gvcfs = gvcfs, gvcf_indices = gvcf_indices, prefix = cohort_name
  }

  # 3) per shard (<=10): loop its intervals -> GenomicsDBImport -> GenotypeGVCFs
  #    -> excise padding -> merge the shard's interval VCFs into one shard VCF
  scatter (shard_tsv in MakeIntervals.shard_tsvs) {
    call GenotypeShard {
      input:
        shard_tsv = shard_tsv, pad = interval_pad,
        sample_name_map = CreateSampleNameMap.sample_name_map,
        ref_fasta = ref_fasta, ref_fasta_fai = ref_fasta_fai, ref_dict = ref_dict,
        gatk_docker = gatk_docker
    }
  }

  # 3) gather shard VCFs (genome order — shards are contiguous & lexically ordered)
  call GatherVcfs {
    input: cohort_name = cohort_name,
           interval_vcfs = GenotypeShard.shard_vcf, interval_vcf_indices = GenotypeShard.shard_vcf_index,
           bcftools_docker = bcftools_docker
  }

  # 4a) VQSR annotate (SNP & indel separate) -> VQSLOD in INFO
  call VQSRAnnotate {
    input: cohort_vcf = GatherVcfs.cohort_vcf, cohort_vcf_index = GatherVcfs.cohort_vcf_index,
           ref_fasta = ref_fasta, ref_fasta_fai = ref_fasta_fai, ref_dict = ref_dict,
           pfcrosses_pass_vcf = pfcrosses_pass_vcf, pfcrosses_pass_vcf_index = pfcrosses_pass_vcf_index,
           vqsr_prior = vqsr_prior, vqsr_max_gaussians = vqsr_max_gaussians, gatk_docker = gatk_docker
  }

  # 4b) VQSLOD<threshold in core, keep non-core, merge -> final multisample VCF
  call FilterMerge {
    input: cohort_name = cohort_name,
           snp_vqslod_vcf = VQSRAnnotate.snp_vqslod_vcf, snp_vqslod_vcf_index = VQSRAnnotate.snp_vqslod_vcf_index,
           indel_vqslod_vcf = VQSRAnnotate.indel_vqslod_vcf, indel_vqslod_vcf_index = VQSRAnnotate.indel_vqslod_vcf_index,
           core_bed = core_bed, vqslod_threshold = vqslod_threshold, bcftools_docker = bcftools_docker
  }

  output {
    File cohort_genotyped_vcf = GatherVcfs.cohort_vcf
    File cohort_genotyped_vcf_index = GatherVcfs.cohort_vcf_index
    File cohort_filtered_vcf  = FilterMerge.filtered_vcf
    File cohort_filtered_vcf_index = FilterMerge.filtered_vcf_index
  }
}

# tile each contig into fixed windows (contig<TAB>start<TAB>end, 1-based inclusive),
# then split the ordered interval list into <=max_shards contiguous shard TSVs.
task MakeIntervals {
  input {
    File ref_fasta_fai
    Int interval_size
    Int max_shards
    String bcftools_docker
  }
  command <<<
    set -euo pipefail
    awk -v W=~{interval_size} 'BEGIN{OFS="\t"}
      { len=$2; s=1; while (s<=len){ e=s+W-1; if(e>len)e=len; print $1,s,e; s=e+1 } }' \
      ~{ref_fasta_fai} > intervals.tsv

    N=$(wc -l < intervals.tsv)
    S=~{max_shards}; if [ "$N" -lt "$S" ]; then S=$N; fi
    PER=$(( (N + S - 1) / S ))                 # ceil(N/S) intervals per shard
    split -d -a 3 -l "$PER" intervals.tsv shard_
    for f in shard_[0-9][0-9][0-9]; do mv "$f" "$f.tsv"; done
    echo "intervals=$N shards=$(ls shard_*.tsv | wc -l) per_shard=$PER" 1>&2
  >>>
  runtime {
    docker: bcftools_docker
    memory: "2 GB"
    disks: "local-disk 10 HDD"
  }
  output {
    File intervals_tsv = "intervals.tsv"
    Array[File] shard_tsvs = glob("shard_*.tsv")
  }
}

# Derive the GenomicsDB sample-name-map from the gVCF headers (mirrors
# broadinstitute/long-read-pipelines SRJointGenotyping.CreateSampleNameMap):
# stream each gVCF (localization_optional -> paths stay gs://), read its one
# sample name with `bcftools query -l`, validate, and write name<TAB>gs://path.
# Runs ONCE; the same map feeds every shard's GenomicsDBImport, which streams the
# gVCFs straight from those gs:// paths (no localization of thousands of files).
task CreateSampleNameMap {
  input {
    Array[File] gvcfs
    Array[File] gvcf_indices
    String prefix
    # lr-basic carries bcftools + gcloud (needed to mint GCS_OAUTH_TOKEN for streaming)
    String docker = "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    Int memory_gb = 3
    Int disk_gb = 20
  }
  parameter_meta {
    gvcfs: { localization_optional: true }
    gvcf_indices: { localization_optional: true }
  }
  command <<<
    set -euxo pipefail
    gvcf_list=~{write_lines(gvcfs)}
    idx_list=~{write_lines(gvcf_indices)}
    # one index per gVCF (the .tbi must sit beside each gVCF for GATK streaming)
    [ "$(wc -l < "$gvcf_list")" -eq "$(wc -l < "$idx_list")" ] || { echo "gvcfs/indices count mismatch" >&2; exit 1; }

    out="~{prefix}.sample_name_map.tsv"; : > "$out"
    GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token); export GCS_OAUTH_TOKEN

    n=0
    while read -r gvcf; do
      bcftools query -l "$gvcf" > sn.txt
      [ "$(wc -l < sn.txt)" -eq 1 ] || { echo "GVCF must contain exactly one sample: $gvcf" >&2; exit 1; }
      grep -qi 'unnamedsample' sn.txt && { echo "GVCF sample name is 'unnamedsample': $gvcf" >&2; exit 1; }
      printf '%s\t%s\n' "$(cat sn.txt)" "$gvcf" >> "$out"
      n=$((n+1))
      # periodically refresh the token so long lists don't outlive it
      if [ "$n" -ge 50 ]; then GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token); export GCS_OAUTH_TOKEN; n=0; fi
    done < "$gvcf_list"
  >>>
  runtime {
    docker: docker
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
  }
  output {
    File sample_name_map = "~{prefix}.sample_name_map.tsv"
  }
}

# One SHARD (GATK only): loop the shard's intervals, and for each run the exact
# Pf7 per-window path (GenomicsDBImport -ip pad -> GenotypeGVCFs Pf7 flags ->
# SelectVariants excise pad), then MergeVcfs the shard's interval VCFs into one.
# ponytail: one JVM spin per interval (~235 intervals/shard at default 10 shards).
#           Correct & Pf7-exact; if the JVM startup cost dominates, batch multiple
#           -L windows into a single GenomicsDBImport/GenotypeGVCFs call.
task GenotypeShard {
  input {
    File shard_tsv
    Int pad
    File sample_name_map          # name<TAB>gs://gvcf ; gVCFs streamed by GATK (not localized)
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    String gatk_docker
    Int memory_gb = 10
    Int disk_gb = 100
  }
  command <<<
    set -euo pipefail

    mkdir -p out
    i=0
    while IFS=$'\t' read -r contig start end; do
      [ -z "${contig:-}" ] && continue
      i=$((i+1)); iv=$(printf 'iv_%05d' "$i")

      rm -rf gdb
      gatk --java-options "-Xmx~{memory_gb - 2}g" GenomicsDBImport \
        --sample-name-map ~{sample_name_map} \
        -L "$contig:$start-$end" -ip ~{pad} \
        --genomicsdb-workspace-path gdb --batch-size 50 --reader-threads 5 --merge-input-intervals

      pad_start=$(( start - ~{pad} )); if [ "$pad_start" -lt 1 ]; then pad_start=1; fi
      gatk --java-options "-Xmx~{memory_gb - 2}g" GenotypeGVCFs \
        -R ~{ref_fasta} -V gendb://gdb \
        -L "$contig:$pad_start-$end" \
        --only-output-calls-starting-in-intervals \
        --use-new-qual-calculator \
        --annotation-group StandardAnnotation --annotation-group AS_StandardAnnotation \
        -O padded.vcf.gz

      # excise the 500bp-before padding -> pure window (Pf7 used `bcftools view -r`)
      gatk SelectVariants -R ~{ref_fasta} -V padded.vcf.gz \
        -L "$contig:$start-$end" -O "out/$iv.vcf.gz"
    done < ~{shard_tsv}

    # merge this shard's interval VCFs into one shard VCF (MergeVcfs sorts + indexes)
    ls -1 out/*.vcf.gz | sort > vcf.list
    gatk --java-options "-Xmx~{memory_gb - 2}g" MergeVcfs -I vcf.list -O shard.vcf.gz
  >>>
  runtime {
    docker: gatk_docker
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
  }
  output {
    File shard_vcf = "shard.vcf.gz"
    File shard_vcf_index = "shard.vcf.gz.tbi"
  }
}

# concat scattered shard VCFs into the cohort VCF (genome order)
task GatherVcfs {
  input {
    String cohort_name
    Array[File] interval_vcfs
    Array[File] interval_vcf_indices
    String bcftools_docker
    Int disk_gb = 50
  }
  command <<<
    set -euo pipefail
    printf '%s\n' ~{sep=' ' interval_vcfs} > vcf.list
    bcftools concat -a -f vcf.list -Oz -o ~{cohort_name}.genotyped.vcf.gz
    bcftools index -t ~{cohort_name}.genotyped.vcf.gz
  >>>
  runtime {
    docker: bcftools_docker
    memory: "4 GB"
    disks: "local-disk " + disk_gb + " HDD"
  }
  output {
    File cohort_vcf = "~{cohort_name}.genotyped.vcf.gz"
    File cohort_vcf_index = "~{cohort_name}.genotyped.vcf.gz.tbi"
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
# core; non-core kept unconditionally; merge -> final multisample VCF.
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
