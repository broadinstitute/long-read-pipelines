version 1.0

import "../../structs/Structs.wdl"

task GenerateBinFiles {
    meta {
        desciption:
        ""
    }

    parameter_meta {

    }
    input {
        File reads
        String prefix = "out"

        String zones
        Int? memory_override
        RuntimeAttr? runtime_attr_override
    }
    output {
        File resouce_monitor_log = "resources.log"
        File log = "hifiasm.log"

        File ec_bin = "~{prefix}.ec.bin"
        File ovlp_reverse_bin = "~{prefix}.ovlp.reverse.bin"
        File ovlp_source_bin = "~{prefix}.ovlp.source.bin"
    }

    command <<<
    set -euxo pipefail

        export MONITOR_MOUNT_POINT="/cromwell_root/"
        bash /opt/vm_local_monitoring_script.sh &> resources.log &
        job_id=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')

        time \
        hifiasm \
            --bin-only \
            -o ~{prefix} \
            -t~{num_cpus} \
            ~{reads} \
            2>&1 | tee hifiasm.log

        if ps -p "${job_id}" > /dev/null; then kill "${job_id}"; fi
        tree -h .
    >>>

    #########################
    Int min_memory = 192 # this min_memory magic number is purely empirical
    Int proposed_memory = 2 * ceil(size(reads, "GiB"))
    Int memory = if proposed_memory < min_memory then min_memory else if proposed_memory > 512 then 512 else proposed_memory
    Int n = memory / 4  # this might be an odd number
    Int num_cpus_proposal = if (n/2)*2 == n then n else n+1  # a hack because WDL doesn't have modulus operator
    Int num_cpus = if num_cpus_proposal > 96 then 96 else num_cpus_proposal

    Int min_disk = 50
    Int half_reads_sz = (1 + ceil(size(reads, "GiB")))/2
    Int proposed_disk = 10 + 8 * half_reads_sz # a trick to do 3 times the reads file size, yet make sure it is even
    Int disk_size = if proposed_disk < min_disk then min_disk else proposed_disk

    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             select_first([memory_override, memory]),
        disk_gb:            disk_size,
        preemptible_tries:  if (size(reads, "GiB") < 24) then 1 else 0,  # a herustic to use preemptible instances for assembling shallow genomes
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-hifiasm:0.19.5"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones: zones
    }
}

task AssembleForHaplotigs {
    parameter_meta {
        bin_files_are_fake: "when true (i.e. this step took in bin files), the bin files are not real, so don't use them"
    }

    input {
        File reads
        String prefix = "out"
        Array[File]? bin_files

        String zones
        RuntimeAttr? runtime_attr_override
    }

    Boolean mv_bin_files = defined(bin_files)
    Array[File] local_bin_files = select_first([bin_files, [reads]])

    command <<<
        set -euxo pipefail

        export MONITOR_MOUNT_POINT="/cromwell_root/"
        bash /opt/vm_local_monitoring_script.sh &> resources.log &
        job_id=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')

        if ~{mv_bin_files} ; then mv ~{sep=" " local_bin_files} . ; ls; fi

        time \
        hifiasm \
            -o ~{prefix} \
            -t~{num_cpus} \
            ~{reads} \
            2>&1 | tee hifiasm.log

        if ps -p "${job_id}" > /dev/null; then kill "${job_id}"; fi
        tree -h .

        # GFA graph to contigs, primary
        # outputs generated this way has "bp" in their names
        time \
        awk '/^S/{print ">"$2; print $3}' \
            ~{prefix}.bp.p_ctg.gfa \
            > ~{prefix}.bp.p_ctg.fa

        ls "~{prefix}.bp.hap"*".p_ctg.gfa"

        # GFA graph to contigs, for each haplotig set
        for haplotype_gfa in ~{prefix}.bp.hap*.p_ctg.gfa; do
            filename=$(basename -- "${haplotype_gfa}")
            haplotype="${filename%.*}"
            time \
            awk '/^S/{print ">"$2; print $3}' \
                "${haplotype_gfa}" \
                > "${haplotype}".fa
        done

        for ff in ./*.p_ctg.fa; do
            bgzip -@4 --index "${ff}" &
        done
        wait
    >>>

    output {
        File resouce_monitor_log = "resources.log"
        File log = "hifiasm.log"

        File hap1_gfa = "~{prefix}.bp.hap1.p_ctg.gfa"
        File hap1_tigs = "~{prefix}.bp.hap1.p_ctg.fa.gz"
        File hap1_tig_gzi = "~{prefix}.bp.hap1.p_ctg.fa.gz.gzi"
        File hap1_lowQ_bed = "~{prefix}.bp.hap1.p_ctg.lowQ.bed"

        File hap2_gfa = "~{prefix}.bp.hap2.p_ctg.gfa"
        File hap2_tigs = "~{prefix}.bp.hap2.p_ctg.fa.gz"
        File hap2_tig_gzi = "~{prefix}.bp.hap2.p_ctg.fa.gz.gzi"
        File hap2_lowQ_bed = "~{prefix}.bp.hap2.p_ctg.lowQ.bed"

        # these are saved, but the one with alt contigs genearted will be preferred for now
        File primary_gfa = "~{prefix}.bp.p_ctg.gfa"
        File primary_fa = "~{prefix}.bp.p_ctg.fa.gz"
        File primary_fa_gzi = "~{prefix}.bp.p_ctg.fa.gz.gzi"

        # raw unitig graph
        File raw_unitig_graph = "~{prefix}.bp.r_utg.gfa"
        File raw_unitig_lowQ_bed = "~{prefix}.bp.r_utg.lowQ.bed"

        # if bin files are provided, then don't delocalize them as that wastes time
        Boolean bin_files_are_fake = ! defined(bin_files)
        File? ec_bin = if defined(bin_files) then "hifiasm.log" else "~{prefix}.ec.bin"
        File? ovlp_reverse_bin = if defined(bin_files) then "hifiasm.log" else "~{prefix}.ovlp.reverse.bin"
        File? ovlp_source_bin = if defined(bin_files) then "hifiasm.log" else "~{prefix}.ovlp.source.bin"
    }

    #########################
    Int min_memory = 192 # this min_memory magic number is purely empirical
    Int proposed_memory = ceil(size(reads, "GiB"))
    Int memory = if proposed_memory < min_memory then min_memory else if proposed_memory > 512 then 512 else proposed_memory
    Int n = memory / 6  # this might be an odd number
    Int num_cpus_proposal = if (n/2)*2 == n then n else n+1  # a hack because WDL doesn't have modulus operator
    Int num_cpus = if num_cpus_proposal > 96 then 96 else num_cpus_proposal

    Int min_disk = 50
    Int half_reads_sz = (1 + ceil(size(reads, "GiB")))/2
    Int inflation_factor = if defined(bin_files) then 10 else 6  # bin files sizes together are ~ 10 + size of raw reads
    Int proposed_disk = inflation_factor * half_reads_sz
    Int disk_size = if proposed_disk < min_disk then min_disk else proposed_disk

    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             memory,
        disk_gb:            disk_size,
        preemptible_tries:  if defined(bin_files) then 1 else 0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-hifiasm:0.19.5"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones: zones
    }
}

task AssembleForAltContigs {

    parameter_meta {
        bin_files:
        "If provided, hifiasm will continue the assembly in the primary-vs-alt mode by reusing the bin files previously generated (likely in the hap mode)."
    }

    input {
        File reads
        String prefix = "out"
        String zones

        Array[File]? bin_files

        RuntimeAttr? runtime_attr_override
    }

    Boolean mv_bin_files = defined(bin_files)
    Array[File] local_bin_files = select_first([bin_files, [reads]])

    command <<<
        set -euxo pipefail

        export MONITOR_MOUNT_POINT="/cromwell_root/"
        bash /opt/vm_local_monitoring_script.sh &> resources.log &
        job_id=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')

        if ~{mv_bin_files} ; then mv ~{sep=" " local_bin_files} . ; ls; fi

        time \
        hifiasm \
            -o ~{prefix} \
            -t~{num_cpus} \
            --primary \
            ~{reads} \
            2>&1 | tee hifiasm.log

        if ps -p "${job_id}" > /dev/null; then kill "${job_id}"; fi
        tree -h .

        date
        # tricky, outputs generated this way has no "bp" in their file names
        # GFA graph to contigs, primary
        awk '/^S/{print ">"$2; print $3}' \
            ~{prefix}.p_ctg.gfa \
        > ~{prefix}.p_ctg.fa &

        # GFA graph to contigs, alternate
        awk '/^S/{print ">"$2; print $3}' \
            ~{prefix}.a_ctg.gfa \
        > ~{prefix}.a_ctg.fa &
        wait
        date

        for ff in ./*_ctg.fa; do
            bgzip -@4 --index "${ff}" &
        done
        wait
    >>>

    output {
        File resouce_monitor_log = "resources.log"
        File log = "hifiasm.log"

        File primary_gfa  = "~{prefix}.p_ctg.gfa"
        File primary_tigs = "~{prefix}.p_ctg.fa.gz"
        File primary_fa_gzi = "~{prefix}.p_ctg.fa.gz.gzi"
        File primary_lowQ_bed = "~{prefix}.p_ctg.lowQ.bed"

        File alternate_gfa  = "~{prefix}.a_ctg.gfa"
        File alternate_tigs = "~{prefix}.a_ctg.fa.gz"
        File alternate_tigs_gzi = "~{prefix}.a_ctg.fa.gz.gzi"
        File alternate_lowQ_bed = "~{prefix}.a_ctg.lowQ.bed"
    }

    #########################
    Int min_memory = 192 # this min_memory magic number is purely empirical
    Int proposed_memory = ceil(size(reads, "GiB"))
    Int memory = if proposed_memory < min_memory then min_memory else if proposed_memory > 512 then 512 else proposed_memory
    Int n = memory / 6  # this might be an odd number
    Int num_cpus_proposal = if (n/2)*2 == n then n else n+1  # a hack because WDL doesn't have modulus operator
    Int num_cpus = if num_cpus_proposal > 96 then 96 else num_cpus_proposal

    Int min_disk = 50
    Int half_reads_sz = (1 + ceil(size(reads, "GiB")))/2
    Int inflation_factor = if defined(bin_files) then 10 else 6  # bin files sizes together are ~ 10 + size of raw reads
    Int proposed_disk = inflation_factor * half_reads_sz
    Int disk_size = if proposed_disk < min_disk then min_disk else proposed_disk

    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             memory,
        disk_gb:            disk_size,
        preemptible_tries:  if defined(bin_files) then 3 else 0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-hifiasm:0.19.5"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones: zones
    }
}
