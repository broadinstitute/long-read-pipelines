version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/ONTUtils.wdl" as ONTUtils
import "tasks/VariantUtils.wdl"
import "tasks/Guppy.wdl" as Guppy
import "tasks/Finalize.wdl" as FF

workflow ONTMethylation {
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
    call Utils.ChunkManifest { input: manifest = ListFilesOfType.manifest, manifest_lines_per_chunk = 2 }

    scatter (manifest_chunk in [ ChunkManifest.manifest_chunks[0], ChunkManifest.manifest_chunks[1] ]) {
        call Megalodon {
            input:
                fast5_files = read_lines(manifest_chunk),
                ref_fasta   = ref_map['fasta'],
                variants    = variants,
                SM          = participant_name
        }
    }

#    call Merge as MergeVariantDBs {
#        input:
#            dbs = Megalodon.per_read_variant_calls_db,
#            merge_type = "variants",
#            runtime_attr_override = { 'mem_gb': 48 }
#    }
#
#    call Merge as MergeModifiedBaseCallDBs {
#        input:
#            dbs = Megalodon.per_read_modified_base_calls_db,
#            merge_type = "modified_bases"
#    }
#
#    call Utils.MergeFastqGzs { input: fastq_gzs = Megalodon.basecalls_fastq, prefix = "basecalls" }

    call Utils.MergeBams as MergeMappings { input: bams = Megalodon.mappings_bam }
    call Utils.MergeBams as MergeModMappings { input: bams = Megalodon.mod_mappings_bam }
    call Utils.MergeBams as MergeVarMappings { input: bams = Megalodon.variant_mappings_bam }

#    call Utils.Cat as CatModifiedBases5mC {
#         input:
#            files = Megalodon.modified_bases_5mC,
#            has_header = false,
#            out = "modified_bases.5mC.bed"
#    }
#
#    call Utils.Cat as CatMappingSummaries {
#        input:
#            files = Megalodon.mappings_summary,
#            has_header = true,
#            out = "mappings_summary.txt"
#    }
#
#    call Utils.Cat as CatSequencingSummaries {
#        input:
#            files = Megalodon.sequencing_summary,
#            has_header = true,
#            out = "sequencing_summary.txt"
#    }

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
                chr = contig
        }
    }

    call VariantUtils.MergePerChrCalls {
        input:
            vcfs = PhaseVariants.phased_vcf_gz,
            ref_dict = ref_map['dict'],
            prefix = "phased.merged"
    }

    call Haplotag {
        input:
            variants_phased = MergePerChrCalls.vcf,
            variants_phased_tbi = MergePerChrCalls.tbi,
            variant_mappings_bam = MergeVarMappings.merged_bam,
            variant_mappings_bai = MergeVarMappings.merged_bai
    }

#    output {
#        #String gcs_basecall_dir = Guppy.gcs_dir
#    }
}

task Megalodon {
    input {
        Array[File] fast5_files
        File ref_fasta
        File variants
        String SM

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
                --sort-mappings \
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

        (samtools view -H $(ls tmp/*/*.bam | head -1) | grep -v '^@RG') && (echo -e '@RG\tID:mappings\tSM:~{SM}') > mappings.header.sam
        (samtools view -H $(ls tmp/*/*.bam | head -1) | grep -v '^@RG') && (echo -e '@RG\tID:mod_mappings\tSM:~{SM}') > mod_mappings.header.sam
        (samtools view -H $(ls tmp/*/*.bam | head -1) | grep -v '^@RG') && (echo -e '@RG\tID:var_mappings\tSM:~{SM}') > variant_mappings.header.sam

        cat tmp/*/basecalls.fastq > megalodon_results/basecalls.fastq

        samtools merge megalodon_results/mappings.bam tmp/*/mappings.bam
        samtools reheader mappings.header.sam megalodon_results/mappings.bam > megalodon_results/mappings.rh.bam
        samtools sort megalodon_results/mappings.rh.bam > megalodon_results/mappings.sorted.bam
        samtools index megalodon_results/mappings.sorted.bam

        samtools merge megalodon_results/mod_mappings.bam tmp/*/mod_mappings.bam
        samtools reheader mod_mappings.header.sam megalodon_results/mod_mappings.bam > megalodon_results/mod_mappings.rh.bam
        samtools sort megalodon_results/mod_mappings.rh.bam > megalodon_results/mod_mappings.sorted.bam
        samtools index megalodon_results/mod_mappings.sorted.bam

        samtools merge megalodon_results/variant_mappings.bam tmp/*/variant_mappings.bam
        samtools reheader variant_mappings.header.sam megalodon_results/variant_mappings.bam > megalodon_results/variant_mappings.rh.bam
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
        boot_disk_gb:       30,
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
        zones:                  ["us-central1-c", "us-central1-f", "us-east1-b", "us-east1-c", "us-west1-a", "us-west1-b"]
        cpuPlatform:            "Intel Haswell"
    }
}

task Merge {
    input {
        Array[File] dbs
        String merge_type

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(dbs, "GB")) + 1

    command <<<
        set -x

        DIRS=$(find /cromwell_root/ -name '*.db' -exec dirname {} \; | tr '\n' ' ')

        megalodon_extras merge ~{merge_type} $DIRS
    >>>

    output {
        File? per_read_variant_calls_db = "megalodon_merge_vars_results/per_read_variant_calls.db"
        File? per_read_modified_base_calls_db = "megalodon_merge_mods_results/per_read_modified_base_calls.db"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
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
        boot_disk_gb:       10,
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
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-whatshap:0.1.1"
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
        String variant_mappings_bam
        File variant_mappings_bai
        String chr

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size([variants, variants_tbi, variant_mappings_bam, variant_mappings_bai], "GB")) + 1

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -hb ~{variant_mappings_bam} ~{chr} > ~{chr}.bam
        samtools index ~{chr}.bam

        # run whatshap with produced mappings and variants
        whatshap phase \
            --distrust-genotypes \
            --ignore-read-groups \
            --chromosome ~{chr} \
            -o variants.phased.~{chr}.vcf \
            ~{variants} \
            ~{chr}.bam

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
        mem_gb:             72,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
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

        tree -h
    >>>

    output {
        File variant_mappings_haplotagged_bam = "variant_mappings.haplotagged.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-whatshap:0.1.1"
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
