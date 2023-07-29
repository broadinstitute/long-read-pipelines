version 1.0

##########################################################################################
# Workflow that runs Nanopolish to polish an ONT genome assembly.
# Quite computationally expensive so there's a parallelization factor parameter.
##########################################################################################

import "../../structs/Structs.wdl"

workflow PolishAssembly {
    input {
        String fast5_dir
        File reads_fasta
        File sequencing_summary # From basecalled output directory
        File draft_assembly_fasta
        Int parallel_instances
    }

    parameter_meta {
        fast5_dir: "Original fast5 reads used to basecall reads_fasta and assemble draft_assembly_fasta"
        reads_fasta: "fastq reads used to assemble draft_assembly_fasta"
        sequencing_summary: "Generated by the basecalling step, look for sequencing_summary.txt in the basecalling output dir"
        draft_assembly_fasta: "Draft assembly to be polished"
        parallel_instances: "Parallelization factor; number of instances to spin up"
    }

    call NanopolishIndex {
        input:
            fast5_dir = fast5_dir,
            reads_fasta = reads_fasta,
            sequencing_summary = sequencing_summary,
            draft_assembly_fasta = draft_assembly_fasta,
            parallel_instances = parallel_instances
    }

    scatter (nanopolish_range in NanopolishIndex.nanopolish_job_ranges) {
        call NanopolishVariants {
            input:
                fast5_and_indexes = NanopolishIndex.fast5_and_indexes,
                draft_assembly_fasta = draft_assembly_fasta,
                draft_assembly_fai = NanopolishIndex.draft_assembly_fai,
                draft_alignment_bam = NanopolishIndex.draft_alignment_bam,
                draft_alignment_bai = NanopolishIndex.draft_alignment_bai,
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
        String fast5_dir
        File reads_fasta
        File sequencing_summary
        File draft_assembly_fasta
        Int parallel_instances

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        fast5_dir:              "GCS path to \"directory\" holding the fast5 files"
        reads_fasta:            "Basecalled reads"
        sequencing_summary:     "summary file produced by the sequencer"
        draft_assembly_fasta:   "FASTA holding the draft assembly"
        parallel_instances:     "how many ways to split the work"
    }

    # The 16 multiplier is to estimate the space needed for the fast5 files
    Int disk_size = 2 * ceil((16 * size(reads_fasta, "GB")) + size(draft_assembly_fasta, "GB") + size(sequencing_summary, "GB"))
    String fast5_dir_cleaned = sub(fast5_dir, "/$", "")
    String draft_basename = basename(draft_assembly_fasta)

    command <<<
        set -euxo pipefail

        mkdir fast5
        gsutil -m cp "~{fast5_dir_cleaned}/*" fast5/
        cp ~{sequencing_summary} fast5/sequencing_summary.txt
        cp ~{reads_fasta} fast5/reads.fasta

        cp ~{draft_assembly_fasta} .
        samtools faidx ~{draft_basename}

        nanopolish index -d fast5/ -s fast5/sequencing_summary.txt fast5/reads.fasta

        minimap2 -ax map-ont -t 8 ~{draft_basename} ~{reads_fasta} | samtools sort -o draft_alignment.sorted.bam
        samtools index draft_alignment.sorted.bam

        python3 /nanopolish/scripts/nanopolish_makerange.py ~{draft_basename} > nanopolish_intervals.txt

        split -n r/~{parallel_instances} -d -a 3 --additional-suffix=".txt" nanopolish_intervals.txt nanopolish_job_
    >>>

    output {
        Array[File] fast5_and_indexes = glob("fast5/*")
        File draft_assembly_fai = "~{draft_basename}.fai"
        File draft_alignment_bam = "draft_alignment.sorted.bam"
        File draft_alignment_bai = "draft_alignment.sorted.bam.bai"

        Array[File] nanopolish_job_ranges =  glob("nanopolish_job_*.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-nanopolish:0.3.0"
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
        Array[File] fast5_and_indexes
        File draft_assembly_fasta
        File draft_assembly_fai
        File draft_alignment_bam
        File draft_alignment_bai

        File nanopolish_range
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        fast5_and_indexes:      ""
        draft_assembly_fasta:   ""
        draft_assembly_fai:     ""
        draft_alignment_bam:    ""
        draft_alignment_bai:    ""
    }

    Int disk_size = 4 * ceil(size(fast5_and_indexes, "GB"))

    command <<<
        set -euxo pipefail

        fast5_dir=$(dirname ~{fast5_and_indexes[0]})
        mv $fast5_dir fast5

        cp ~{draft_assembly_fasta} draft_assembly.fasta
        cp ~{draft_assembly_fai} draft_assembly.fasta.fai

        cat ~{nanopolish_range} | parallel --results nanopolish.results -P 4 \
            nanopolish variants --consensus -o polished.{1}.vcf \
              --max-haplotypes=5000 \
              -w {1} \
              -r fast5/reads.fasta \
              -b ~{draft_alignment_bam} \
              -g draft_assembly.fasta \
              -t 2 \
              --min-candidate-frequency 0.1
    >>>

    output {
        Array[File] vcfs = glob("polished.*.vcf")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        # We should try to find a way to support this to reduce costs.
        # Right now, it is almost always preempted with the amount of
        # time it takes to localize the fast5 files
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-nanopolish:0.3.0"
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

    parameter_meta {
        vcfs:                   ""
        draft_assembly_fasta:   ""
    }

    Int disk_size = 2 * ceil(size(vcfs, "GB") + size(draft_assembly_fasta, "GB"))

    command <<<
        set -euxo pipefail

        nanopolish vcf2fasta -g ~{draft_assembly_fasta} ~{sep=" " vcfs} > polished_assembly.fasta
    >>>

    output {
       File polished_assembly = "polished_assembly.fasta"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-nanopolish:0.3.0"
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