version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/VariantUtils.wdl"
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ONTMethylation {

    meta {
        description: "ONT Methylation pipeline"
    }
    parameter_meta {
        gcs_fast5_dir: "GCS directory containing fast5 files"
        ref_map_file: "Reference map file"
        variants: "VCF file containing variants"
        variants_tbi: "Tabix index for VCF file"
        participant_name: "Participant name"
        prefix: "Prefix for output files"
        gcs_out_root_dir: "GCS directory to write output files"
    }

    input {
        String gcs_fast5_dir

        File ref_map_file
        File variants
        File variants_tbi

        String participant_name
        String prefix

        String gcs_out_root_dir
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTMethylation/~{prefix}"

    call Utils.ListFilesOfType { input: gcs_dir = gcs_fast5_dir, suffixes = [ ".fast5" ] }
    call Utils.ChunkManifest { input: manifest = ListFilesOfType.manifest, manifest_lines_per_chunk = 30 }

    scatter (manifest_chunk in ChunkManifest.manifest_chunks) {
        call Megalodon {
            input:
                fast5_files = read_lines(manifest_chunk),
                ref_fasta   = ref_map['fasta'],
                variants    = variants,
                SM          = participant_name
        }
    }

    call MergeVariantDBs { input: dbs = Megalodon.per_read_variant_calls_db }
    call MergeModifiedBaseCallDBs { input: dbs = Megalodon.per_read_modified_base_calls_db }

    call Utils.MergeFastqs { input: fastqs = Megalodon.basecalls_fastq, prefix = "basecalls" }

    call Utils.Cat as CatModifiedBases5mC {
        input:
            files = Megalodon.modified_bases_5mC,
            has_header = false,
            out = "modified_bases.5mC.bed"
    }

    call Utils.Cat as CatMappingSummaries {
        input:
            files = Megalodon.mappings_summary,
            has_header = true,
            out = "mappings_summary.txt"
    }

    call Utils.Cat as CatSequencingSummaries {
        input:
            files = Megalodon.sequencing_summary,
            has_header = true,
            out = "sequencing_summary.txt"
    }

    call Utils.MergeBams as MergeMappings { input: bams = Megalodon.mappings_bam }
    call Utils.MergeBams as MergeModMappings { input: bams = Megalodon.mod_mappings_bam }
    call Utils.MergeBams as MergeVarMappings { input: bams = Megalodon.variant_mappings_bam }

    call WhatsHapFilter { input: variants = variants, variants_tbi = variants_tbi }
    call IndexVariants { input: variants = WhatsHapFilter.whatshap_filt_vcf }

    call Utils.MakeChrIntervalList {
        input:
            ref_dict = ref_map['dict'],
            filter = ['GL', 'JH']
    }

    scatter (c in MakeChrIntervalList.chrs) {
        String contig = c[0]

        call Utils.SubsetBam {
            input:
                bam = MergeVarMappings.merged_bam,
                bai = MergeVarMappings.merged_bai,
                locus = contig,
                prefix = contig
        }

        call VariantUtils.SubsetVCF {
            input:
                vcf_gz = IndexVariants.vcf_gz,
                vcf_tbi = IndexVariants.vcf_tbi,
                locus = contig,
                prefix = contig
        }

        call PhaseVariants {
             input:
                variants = SubsetVCF.subset_vcf,
                variants_tbi = SubsetVCF.subset_tbi,
                variant_mappings_bam = SubsetBam.subset_bam,
                variant_mappings_bai = SubsetBam.subset_bai,
                ref_fasta = ref_map['fasta'],
                ref_fai   = ref_map['fai'],
                chr = contig
        }

        call Haplotag {
            input:
                variants_phased = PhaseVariants.phased_vcf_gz,
                variants_phased_tbi = PhaseVariants.phased_vcf_tbi,
                variant_mappings_bam = MergeVarMappings.merged_bam,
                variant_mappings_bai = MergeVarMappings.merged_bai
        }
    }

    call VariantUtils.MergePerChrCalls {
        input:
            vcfs = PhaseVariants.phased_vcf_gz,
            ref_dict = ref_map['dict'],
            prefix = "phased.merged"
    }

    call Utils.MergeBams as MergeHaplotagBams { input: bams = Haplotag.variant_mappings_haplotagged_bam }

    call ExtractHaplotypeReads {
         input:
            haplotagged_bam = MergeHaplotagBams.merged_bam,
            phased_variants_vcf = MergePerChrCalls.vcf,
            phased_variants_tbi = MergePerChrCalls.tbi,
            per_read_variant_calls_db = MergeVariantDBs.per_read_variant_calls_db,
            per_read_modified_base_calls_db = MergeModifiedBaseCallDBs.per_read_modified_base_calls_db
    }

    call CallHaploidVariants as CallHaplotype1Variants {
        input:
            per_read_variant_calls_db = MergeVariantDBs.per_read_variant_calls_db,
            per_read_modified_base_calls_db = MergeModifiedBaseCallDBs.per_read_modified_base_calls_db,
            read_ids = ExtractHaplotypeReads.haplotype_1_read_ids,
            suffix = 1
    }

    call CallHaploidVariants as CallHaplotype2Variants {
        input:
            per_read_variant_calls_db = MergeVariantDBs.per_read_variant_calls_db,
            per_read_modified_base_calls_db = MergeModifiedBaseCallDBs.per_read_modified_base_calls_db,
            read_ids = ExtractHaplotypeReads.haplotype_2_read_ids,
            suffix = 2
    }

    # Finalize
    String vdir = outdir + "/variants/phased"
    String adir = outdir + "/alignments/haplotagged"

    File vcf = MergePerChrCalls.vcf
    File tbi = MergePerChrCalls.tbi
    File bam = MergeHaplotagBams.merged_bam
    File bai = MergeHaplotagBams.merged_bai

    call FF.FinalizeToFile as FinalizePhasedVcf { input: outdir = vdir, file = vcf, name = "~{participant_name}.phased.vcf.gz" }
    call FF.FinalizeToFile as FinalizePhasedTbi { input: outdir = vdir, file = tbi, name = "~{participant_name}.phased.vcf.gz.tbi" }

    call FF.FinalizeToFile as FinalizeModMappedBam { input: outdir = adir, file = MergeModMappings.merged_bam, name = "~{participant_name}.mod_mapped.bam" }
    call FF.FinalizeToFile as FinalizeModMappedBai { input: outdir = adir, file = MergeModMappings.merged_bai, name = "~{participant_name}.mod_mapped.bam.bai" }

    call FF.FinalizeToFile as FinalizeHaplotaggedBam { input: outdir = adir, file = bam, name = "~{participant_name}.haplotagged.bam" }
    call FF.FinalizeToFile as FinalizeHaplotaggedBai { input: outdir = adir, file = bai, name = "~{participant_name}.haplotagged.bam.bai" }

    call FF.FinalizeToFile as FinalizeHaplotype1Vcf { input: outdir = vdir, file = CallHaplotype1Variants.haploid_vcf, name = "~{participant_name}.aggregated.haplotype1.vcf.gz" }
    call FF.FinalizeToFile as FinalizeHaplotype1Tbi { input: outdir = vdir, file = CallHaplotype1Variants.haploid_tbi, name = "~{participant_name}.aggregated.haplotype1.vcf.gz.tbi" }
    call FF.FinalizeToFile as FinalizeHaplotype2Vcf { input: outdir = vdir, file = CallHaplotype2Variants.haploid_vcf, name = "~{participant_name}.aggregated.haplotype2.vcf.gz" }
    call FF.FinalizeToFile as FinalizeHaplotype2Tbi { input: outdir = vdir, file = CallHaplotype2Variants.haploid_tbi, name = "~{participant_name}.aggregated.haplotype2.vcf.gz.tbi" }

    output {
        File phased_vcf = FinalizePhasedVcf.gcs_path
        File phased_tbi = FinalizePhasedTbi.gcs_path

        File mod_mapped_bam = FinalizeModMappedBam.gcs_path
        File mod_mapped_bai = FinalizeModMappedBai.gcs_path

        File haplotagged_bam = FinalizeHaplotaggedBam.gcs_path
        File haplotagged_bai = FinalizeHaplotaggedBai.gcs_path

        File haplotype1_vcf = FinalizeHaplotype1Vcf.gcs_path
        File haplotype1_tbi = FinalizeHaplotype1Tbi.gcs_path
        File haplotype2_vcf = FinalizeHaplotype2Vcf.gcs_path
        File haplotype2_tbi = FinalizeHaplotype2Tbi.gcs_path
    }
}

task Megalodon {
    input {
        Array[File] fast5_files
        File ref_fasta
        File variants
        String SM

        Array[String] zones = ["us-central1-c", "us-central1-f", "us-central1-a", "us-central1-b"]

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4 * ceil(size(fast5_files, "GB") + size([ref_fasta, variants], "GB"))

    command <<<
        set -euxo pipefail

        num_cores=$(grep -c '^processor' /proc/cpuinfo | awk '{ print $1 - 1 }')
        dir=$(dirname ~{fast5_files[0]})

        for fast5 in $dir/*.fast5; do
            BN="$(basename "$fast5" .fast5)"

            TMP_DIR=tmp/$BN
            mkdir -p $TMP_DIR

            mkdir f5
            mv $fast5 f5/

            megalodon f5 \
                --guppy-params "-d /rerio/basecall_models" \
                --guppy-config res_dna_r941_prom_modbases_5mC_CpG_v001.cfg \
                --outputs basecalls mappings mod_mappings mods variant_mappings \
                --reference ~{ref_fasta} \
                --mod-motif m CG 0 \
                --variant-filename ~{variants} \
                --devices cuda:0 \
                --processes $num_cores \
                --guppy-server-path /usr/bin/guppy_basecall_server \
                --output-directory $TMP_DIR \
                --overwrite \
                --suppress-progress-bars \
                --suppress-queues-status

            rm -rf f5
        done

        mkdir megalodon_results

        cat tmp/*/basecalls.fastq > megalodon_results/basecalls.fastq

        ((samtools dict ~{ref_fasta}) && (echo -e '@RG\tID:mappings\tSM:~{SM}')) > mappings.header.sam
        samtools cat -o megalodon_results/mappings.bam tmp/*/mappings.bam
        samtools addreplacerg -O BAM -r '@RG\tID:mappings\tSM:~{SM}' megalodon_results/mappings.bam | samtools reheader mappings.header.sam - > megalodon_results/mappings.rh.bam
        samtools sort megalodon_results/mappings.rh.bam > megalodon_results/mappings.sorted.bam
        samtools index megalodon_results/mappings.sorted.bam

        ((samtools dict ~{ref_fasta}) && (echo -e '@RG\tID:mod_mappings\tSM:~{SM}')) > mod_mappings.header.sam
        samtools cat -o megalodon_results/mod_mappings.bam tmp/*/mod_mappings.bam
        samtools addreplacerg -O BAM -r '@RG\tID:mod_mappings\tSM:~{SM}' megalodon_results/mod_mappings.bam | samtools reheader mod_mappings.header.sam - > megalodon_results/mod_mappings.rh.bam
        samtools sort megalodon_results/mod_mappings.rh.bam > megalodon_results/mod_mappings.sorted.bam
        samtools index megalodon_results/mod_mappings.sorted.bam

        ((samtools dict ~{ref_fasta}) && (echo -e '@RG\tID:variant_mappings\tSM:~{SM}')) > variant_mappings.header.sam
        samtools cat -o megalodon_results/variant_mappings.bam tmp/*/variant_mappings.bam
        samtools addreplacerg -O BAM -r '@RG\tID:variant_mappings\tSM:~{SM}' megalodon_results/variant_mappings.bam | samtools reheader variant_mappings.header.sam - > megalodon_results/variant_mappings.rh.bam
        samtools sort megalodon_results/variant_mappings.rh.bam > megalodon_results/variant_mappings.sorted.bam
        samtools index megalodon_results/variant_mappings.sorted.bam

        megalodon_extras merge variants tmp/*
        megalodon_extras merge modified_bases tmp/*

        cat tmp/*/modified_bases.5mC.bed > megalodon_results/modified_bases.5mC.bed

        ((head -1 $(ls tmp/*/mappings.summary.txt | head -1)) \
            && (tail -q -n +2 tmp/*/mappings.summary.txt)) \
            > megalodon_results/mappings.summary.txt

        ((head -1 $(ls tmp/*/sequencing_summary.txt | head -1)) \
            && (tail -q -n +2 tmp/*/sequencing_summary.txt)) \
            > megalodon_results/sequencing_summary.txt
    >>>

    output {
        File basecalls_fastq = "megalodon_results/basecalls.fastq"

        File mappings_bam = "megalodon_results/mappings.sorted.bam"
        File mappings_bai = "megalodon_results/mappings.sorted.bam.bai"

        File mod_mappings_bam = "megalodon_results/mod_mappings.sorted.bam"
        File mod_mappings_bai = "megalodon_results/mod_mappings.sorted.bam.bai"

        File variant_mappings_bam = "megalodon_results/variant_mappings.sorted.bam"
        File variant_mappings_bai = "megalodon_results/variant_mappings.sorted.bam.bai"

        File modified_bases_5mC = "megalodon_results/modified_bases.5mC.bed"
        File mappings_summary = "megalodon_results/mappings.summary.txt"
        File sequencing_summary = "megalodon_results/sequencing_summary.txt"

        File per_read_modified_base_calls_db = "megalodon_merge_mods_results/per_read_modified_base_calls.db"
        File per_read_variant_calls_db = "megalodon_merge_vars_results/per_read_variant_calls.db"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             50,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-megalodon:2.3.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        gpuType:                "nvidia-tesla-p100"
        gpuCount:               1
        nvidiaDriverVersion:    "418.152.00"
        zones:                  "~{sep=' ' zones}"
        cpuPlatform:            "Intel Haswell"
    }
}

task MergeVariantDBs {
    input {
        Array[File] dbs

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(dbs, "GB")) + 1

    command <<<
        set -x

        DIRS=$(find /cromwell_root/ -name '*.db' -exec dirname {} \; | tr '\n' ' ')

        megalodon_extras merge variants $DIRS
    >>>

    output {
        File per_read_variant_calls_db = "megalodon_merge_vars_results/per_read_variant_calls.db"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             48,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-megalodon:2.3.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MergeModifiedBaseCallDBs {
    input {
        Array[File] dbs

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(dbs, "GB")) + 1

    command <<<
        set -x

        DIRS=$(find /cromwell_root/ -name '*.db' -exec dirname {} \; | tr '\n' ' ')

        megalodon_extras merge modified_bases $DIRS
    >>>

    output {
        File per_read_modified_base_calls_db = "megalodon_merge_mods_results/per_read_modified_base_calls.db"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             48,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-megalodon:2.3.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task WhatsHapFilter {
    input {
        File variants
        File variants_tbi

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 100*ceil(size([variants, variants_tbi], "GB")) + 1

    command <<<
        set -x

        gunzip -c ~{variants} > variants.sorted.vcf

        # filter whatshap incompatible variants and create indices
        megalodon_extras \
            phase_variants whatshap_filter \
            variants.sorted.vcf \
            variants.sorted.whatshap_filt.vcf \
            --filtered-records whatshap_filt.txt
    >>>

    output {
        File whatshap_filt_vcf = "variants.sorted.whatshap_filt.vcf"
        File whatshap_filt_txt = "whatshap_filt.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-megalodon:2.3.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task IndexVariants {
    input {
        File variants

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(variants, "GB")) + 1

    command <<<
        set -euxo pipefail

        mv ~{variants} variants.vcf

        bgzip variants.vcf
        tabix variants.vcf.gz
    >>>

    output {
        File vcf_gz = "variants.vcf.gz"
        File vcf_tbi = "variants.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-whatshap:1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task PhaseVariants {
    input {
        File variants
        File variants_tbi
        File variant_mappings_bam
        File variant_mappings_bai
        File ref_fasta
        File ref_fai
        String chr

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size([variants, variants_tbi, variant_mappings_bam, variant_mappings_bai], "GB")) + 1

    command <<<
        set -euxo pipefail

        # run whatshap with produced mappings and variants
        whatshap phase \
            --reference ~{ref_fasta} \
            --distrust-genotypes \
            --ignore-read-groups \
            --chromosome ~{chr} \
            -o variants.phased.~{chr}.vcf \
            ~{variants} \
            ~{variant_mappings_bam}

        # assign haplotypes to reads
        bgzip variants.phased.~{chr}.vcf
        tabix variants.phased.~{chr}.vcf.gz
    >>>

    output {
        File phased_vcf_gz = "variants.phased.~{chr}.vcf.gz"
        File phased_vcf_tbi = "variants.phased.~{chr}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-whatshap:0.13"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Haplotag {
    input {
        File variants_phased
        File variants_phased_tbi
        File variant_mappings_bam
        File variant_mappings_bai

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size([variants_phased, variants_phased_tbi, variant_mappings_bam, variant_mappings_bai], "GB")) + 1

    command <<<
        set -euxo pipefail

        whatshap \
            haplotag ~{variants_phased} \
            ~{variant_mappings_bam} \
            -o variant_mappings.haplotagged.bam

        find . -type f -exec ls -lah {} \;
    >>>

    output {
        File variant_mappings_haplotagged_bam = "variant_mappings.haplotagged.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-whatshap:1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ExtractHaplotypeReads {
    input {
        File haplotagged_bam
        File phased_variants_vcf
        File phased_variants_tbi
        File per_read_variant_calls_db
        File per_read_modified_base_calls_db

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size([haplotagged_bam, phased_variants_vcf, phased_variants_tbi, per_read_variant_calls_db, per_read_modified_base_calls_db], "GB")) + 1

    command <<<
        set -x

        nproc=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        mkdir out_dir
        mv ~{per_read_variant_calls_db} out_dir
        mv ~{per_read_modified_base_calls_db} out_dir

        # extract haplotype reads and call haploid variants
        megalodon_extras \
            phase_variants extract_haplotype_reads \
            ~{haplotagged_bam} \
            out_dir/variant_mappings
    >>>

    output {
        File haplotype_1_read_ids = "out_dir/variant_mappings.haplotype_1_read_ids.txt"
        File haplotype_2_read_ids = "out_dir/variant_mappings.haplotype_2_read_ids.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-megalodon:2.3.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CallHaploidVariants {
    input {
        File per_read_variant_calls_db
        File per_read_modified_base_calls_db
        File read_ids
        Int suffix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size([per_read_variant_calls_db, per_read_modified_base_calls_db], "GB")) + 1

    command <<<
        set -x

        nproc=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        mkdir out_dir
        mv ~{per_read_variant_calls_db} out_dir
        mv ~{per_read_modified_base_calls_db} out_dir

        # call haploid variants
        megalodon_extras \
            aggregate run \
            --megalodon-directory out_dir --output-suffix haplotype_~{suffix} \
            --read-ids-filename ~{read_ids} \
            --outputs variants --haploid --processes $nproc

        tree -h
    >>>

    output {
        File haploid_vcf = "out_dir/variants.haplotype_~{suffix}.sorted.vcf.gz"
        File haploid_tbi = "out_dir/variants.haplotype_~{suffix}.sorted.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-megalodon:2.3.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
