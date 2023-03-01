version 1.0

import "Structs.wdl"

task Quantify {
    input {
        File aligned_bam
        File aligned_bai
        File gtf
        Boolean keep_retained_introns = false
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size([aligned_bam, aligned_bai, gtf], "GB"))

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        stringtie \
            -Lv -p $num_core ~{true='-i' false='' keep_retained_introns} \
            -G ~{gtf} \
            -o ~{prefix}.gtf \
            -A ~{prefix}.gene_abund.out \
            ~{aligned_bam}
    >>>

    output {
        File st_gtf = "~{prefix}.gtf"
        File st_abund = "~{prefix}.gene_abund.out"
    }

    # TODO: Debug memory here.  Is getting seg fault...  More memory doesn't help.
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-stringtie2:2.2.1"
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

task ExtractTranscriptSequences {
    input {
        File ref_fasta
        File ref_fasta_fai
        File gtf
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([ref_fasta, ref_fasta_fai, gtf], "GB"))

    command <<<
        set -euxo pipefail

        gffread -w ~{prefix}.fa -g ~{ref_fasta} ~{gtf}
        samtools faidx ~{prefix}.fa
        samtools dict ~{prefix}.fa > ~{prefix}.dict
    >>>

    output {
        File transcripts_fa = "~{prefix}.fa"
        File transcripts_fai = "~{prefix}.fa.fai"
        File transcripts_dict = "~{prefix}.dict"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-stringtie2:2.2.1"
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

task CompareTranscriptomes {
    input {
        File guide_gtf
        File new_gtf
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([guide_gtf, new_gtf], "GB"))

    # TODO: Look into the comment below and figure out why this is happening.
    command <<<
        set -euxo pipefail

        gffcompare -R -r ~{guide_gtf} -o ~{prefix} ~{new_gtf}

        dir=$(dirname ~{new_gtf})
        mv $dir/*.refmap ~{prefix}.refmap
        mv $dir/*.tmap ~{prefix}.tmap

        if [ ! -e ~{prefix}.stats ] ; then
            # Sometimes the stats file hasn't been named `.stats` and just has the prefix name.
            # I have no idea why.  It's very weird.  But we don't use this for our output anyway, so
            # for now we'll try to use the prefix file itself if it exists:
            mv ~{prefix} ~{prefix}.stats
        fi

        tree -h
    >>>

    output {
        File annotated_gtf = "~{prefix}.annotated.gtf"
        File loci = "~{prefix}.loci"
        File tracking = "~{prefix}.tracking"
        File refmap = "~{prefix}.refmap"
        File tmap = "~{prefix}.tmap"
        File stats = "~{prefix}.stats"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-stringtie2:2.2.1"
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