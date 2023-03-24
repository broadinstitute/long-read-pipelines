version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/Guppy.wdl" as Guppy
import "tasks/Finalize.wdl" as FF

workflow ONTMethylation {
    input {
        String gcs_fast5_dir

        File ref_map_file

        String participant_name
        String prefix

        String gcs_out_root_dir
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTMethylation/~{prefix}"

    String fast5_dir = sub(gcs_fast5_dir,"/$", "") + "/" + participant_name

    call Utils.ListFilesOfType { input: gcs_dir = fast5_dir, suffixes = [ ".fast5" ] }
    call Utils.ChunkManifest { input: manifest = ListFilesOfType.manifest, manifest_lines_per_chunk = 30 }

    scatter (manifest_chunk in ChunkManifest.manifest_chunks) {
        call Megalodon {
            input:
                fast5_files = read_lines(manifest_chunk),
                ref_fasta   = ref_map['fasta'],
                SM          = participant_name
        }
    }

    call Utils.MergeFastqGzs { input: fastq_gzs = Megalodon.basecalls_fastq, prefix = "basecalls" }

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

    # Finalize
    String adir = outdir + "/alignments"

    call FF.FinalizeToFile as FinalizeModMappedBam { input: outdir = adir, file = MergeModMappings.merged_bam, name = "~{participant_name}.mod_mapped.bam" }
    call FF.FinalizeToFile as FinalizeModMappedBai { input: outdir = adir, file = MergeModMappings.merged_bai, name = "~{participant_name}.mod_mapped.bam.bai" }

    call FF.FinalizeToFile as FinalizeMappedBam { input: outdir = adir, file = MergeMappings.merged_bam, name = "~{participant_name}.mapped.bam" }
    call FF.FinalizeToFile as FinalizeMappedBai { input: outdir = adir, file = MergeMappings.merged_bai, name = "~{participant_name}.mapped.bam.bai" }

    call FF.FinalizeToFile as FinalizeModifiedBases5mC { input: outdir = outdir, file = CatModifiedBases5mC.combined, name = "~{participant_name}.modified_bases.5mc.bed" }
    call FF.FinalizeToFile as FinalizeMappingSummary { input: outdir = outdir, file = CatMappingSummaries.combined, name = "~{participant_name}.mappings_summary.txt" }
    call FF.FinalizeToFile as FinalizeSequencingSummary { input: outdir = outdir, file = CatSequencingSummaries.combined, name = "~{participant_name}.sequencing_summary.txt" }
    call FF.FinalizeToFile as FinalizeBasecalls { input: outdir = outdir, file = MergeFastqGzs.merged_fastq_gz, name = "~{participant_name}.merged_fastq.gz" }

    output {
        File mod_mapped_bam = FinalizeModMappedBam.gcs_path
        File mod_mapped_bai = FinalizeModMappedBai.gcs_path
        File mapped_bam = FinalizeMappedBam.gcs_path
        File mapped_bai = FinalizeMappedBai.gcs_path

        File basecalls = FinalizeBasecalls.gcs_path
    }
}

task Megalodon {
    input {
        Array[File] fast5_files
        File ref_fasta
        String SM

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4 * ceil(size(fast5_files, "GB") + size(ref_fasta, "GB"))

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
                --outputs basecalls mappings mod_mappings mods \
                --reference ~{ref_fasta} \
                --mod-motif m CG 0 \
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

        samtools cat -o megalodon_results/mappings.bam tmp/*/mappings.bam
        samtools sort megalodon_results/mappings.bam > megalodon_results/mappings.sorted.bam
        samtools index megalodon_results/mappings.sorted.bam

        samtools cat -o megalodon_results/mod_mappings.bam tmp/*/mod_mappings.bam
        samtools sort megalodon_results/mod_mappings.bam > megalodon_results/mod_mappings.sorted.bam
        samtools index megalodon_results/mod_mappings.sorted.bam

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

        File modified_bases_5mC = "megalodon_results/modified_bases.5mC.bed"
        File mappings_summary = "megalodon_results/mappings.summary.txt"
        File sequencing_summary = "megalodon_results/sequencing_summary.txt"

        File per_read_modified_base_calls_db = "megalodon_merge_mods_results/per_read_modified_base_calls.db"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             50,
        disk_gb:            disk_size,
        boot_disk_gb:       30,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/ymostovoy/lr-megalodon:2.5.0"
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
        megalodon_extras per_read_text modified_bases megalodon_merge_mods_results/

    >>>

    output {
        File per_read_modified_base_calls_db = "megalodon_merge_mods_results/per_read_modified_base_calls.db"
        File per_read_modified_base_calls_txt = "megalodon_merge_mods_results/per_read_modified_base_calls.txt"

    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             48,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/ymostovoy/lr-megalodon:2.5.0"
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
