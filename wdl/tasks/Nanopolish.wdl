version 1.0

import "Structs.wdl"

workflow PolishAssembly {
    input {
        String fast5_gcs_dir
        File combined_read_fasta
        File sequencing_summary
        File draft_assembly_fasta
        File output_dir
    }

    call NanopolishIndex {
        input:
            fast5_gcs_dir = fast5_gcs_dir,
            combined_read_fasta = combined_read_fasta,
            sequencing_summary = sequencing_summary,
            draft_assembly_fasta = draft_assembly_fasta
    }

    scatter (nanopolish_range in NanopolishIndex.nanopolish_job_ranges) {
        call NanopolishVariants {
            input:
                fast5_tar = NanopolishIndex.fast5_tar,
                draft_assembly_fasta = draft_assembly_fasta,
                draft_alignment_bam = NanopolishIndex.draft_alignment_bam,
                draft_alignment_bai = NanopolishIndex.draft_alignment_bai,
                combined_read_fasta = combined_read_fasta,
                nanopolish_range = nanopolish_range
        }
    }

    call MergeVcfs {
        input:
            vcfs = flatten(NanopolishVariants.vcfs),
            draft_assembly_fasta = draft_assembly_fasta
    }

    output {
        File polished_assembly = MergeVcfs.polished_assembly
    }
}

task NanopolishIndex {
    input {
        String fast5_gcs_dir
        File combined_read_fasta
        File sequencing_summary
        File draft_assembly_fasta

        RuntimeAttr? runtime_attr_override
    }

    String fast5_dir = sub(fast5_gcs_dir, "/$", "")

    command <<<
        set -euxo pipefail

        mkdir fast5
        gsutil -m cp ~{fast5_dir}/* fast5/
        cp ~{sequencing_summary} fast5/sequencing_summary.txt

        nanopolish index -d fast5/ -s fast5/sequencing_summary.txt ~{combined_read_fasta}

        minimap2 -ax map-ont -t 8 ~{draft_assembly_fasta} ~{combined_read_fasta} | samtools sort -o draft_alignment.sorted.bam
        samtools index draft_alignment.sorted.bam

        # python3 /nanopolish/scripts/nanopolish_makerange.py ~{draft_assembly_fasta} > nanopolish_intervals.txt

        tar -cf fast5.tar fast5

        echo "tig00000001:250000-300200" >> nanopolish_job_1_ranges.txt
        echo "tig00000611:800000-850200" >> nanopolish_job_2_ranges.txt
        echo "tig00000618:0-30854" >> nanopolish_job_2_ranges.txt
    >>>

    output {
        File fast5_tar = "fast5.tar"
        File draft_alignment_bam = "draft_alignment.sorted.bam"
        File draft_alignment_bai = "draft_alignment.sorted.bam.bai"

        Array[File] nanopolish_job_ranges =  glob("nanopolish_job_*_ranges.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             30,
        disk_gb:            250, # this should be dynamic
        boot_disk_gb:       10,
        preemptible_tries:  0, # change
        max_retries:        0, # change
        docker:             "quay.io/broad-long-read-pipelines/lr-nanopolish:0.3.0"
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

task NanopolishVariants {
    input {
        File fast5_tar
        File draft_assembly_fasta
        File draft_alignment_bam
        File draft_alignment_bai
        File combined_read_fasta

        File nanopolish_range
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        ls -ahl

        cat ~{nanopolish_range} | parallel --results nanopolish.results -P 4 \
            nanpolish variants --consensus -o polished.{1}.vcf \
              -w {1} \
              -r ~{combined_read_fasta} \
              -b ~{draft_alignment_bam} \
              -g ~{draft_assembly_fasta} \
              -t 2 \
              --min-candidate-frequency 0.1
    >>>

    output {
        Array[File] vcfs = glob("polished.*.vcf")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            250,
        boot_disk_gb:       10,
        preemptible_tries:  0, #change
        max_retries:        0, #change
        docker:             "quay.io/broad-long-read-pipelines/lr-nanopolish:0.1.0"
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

task MergeVcfs {
    input {
        Array[File] vcfs
        File draft_assembly_fasta

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        ls -ahl

        nanopolish vcf2fasta -g ~{draft_assembly_fasta} ~{sep=" " vcfs} > polished_assembly.fasta
    >>>

    output {
       File polished_assembly = "polished_assembly.fasta"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             20,
        disk_gb:            250,
        boot_disk_gb:       10,
        preemptible_tries:  0, #change
        max_retries:        0, #change
        docker:             "quay.io/broad-long-read-pipelines/lr-nanopolish:0.1.0"
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