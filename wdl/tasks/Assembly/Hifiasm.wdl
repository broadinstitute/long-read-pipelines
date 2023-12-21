version 1.0

import "../../structs/Structs.wdl"

task AssembleForHaplotigs {
    input {
        File reads
        String prefix = "out"
        String zones

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        export MONITOR_MOUNT_POINT="/cromwell_root/"
        bash /opt/vm_local_monitoring_script.sh &> resources.log &
        job_id=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')

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
        awk '/^S/{print ">"$2; print $3}' \
            ~{prefix}.bp.p_ctg.gfa \
            > ~{prefix}.bp.p_ctg.fa

        ls "~{prefix}.bp.hap"*".p_ctg.gfa"

        # GFA graph to contigs, for each haplotig set
        for haplotype_gfa in ~{prefix}.bp.hap*.p_ctg.gfa; do
            filename=$(basename -- "${haplotype_gfa}")
            haplotype="${filename%.*}"
            awk '/^S/{print ">"$2; print $3}' \
                "${haplotype_gfa}" \
                > "${haplotype}".fa
        done

        for ff in ./*.p_ctg.fa; do
            bgzip -@2 -k --index "${ff}"
        done
    >>>

    output {
        File resouce_monitor_log = "resources.log"
        File log = "hifiasm.log"

        File hap1_gfa = "~{prefix}.bp.hap1.p_ctg.gfa"
        File hap1_tigs = "~{prefix}.bp.hap1.p_ctg.fa.gz"
        File hap1_tig_gzi = "~{prefix}.bp.hap1.p_ctg.fa.gz.gzi"

        File hap2_gfa = "~{prefix}.bp.hap2.p_ctg.gfa"
        File hap2_tigs = "~{prefix}.bp.hap2.p_ctg.fa.gz"
        File hap2_tig_gzi = "~{prefix}.bp.hap2.p_ctg.fa.gz.gzi"

        # these are saved, but the one with alt contigs genearted will be preferred for now
        File primary_gfa = "~{prefix}.bp.p_ctg.gfa"
        File primary_fa = "~{prefix}.bp.p_ctg.fa.gz"
        File primary_fa_gzi = "~{prefix}.bp.p_ctg.fa.gz.gzi"
    }

    #########################
    Int proposed_memory = 4 * ceil(size(reads, "GB"))
    Int memory = if proposed_memory < 96 then 96 else proposed_memory  # this 96 magic number is purely empirical
    Int n = memory / 4  # this might be an odd number
    Int num_cpus_proposal = if (n/2)*2 == n then n else n+1  # a hack because WDL doesn't have modulus operator
    Int num_cpus = if num_cpus_proposal > 96 then 96 else num_cpus_proposal

    Int min_disk = 75
    Int proposed_disk = 5 * ceil(size(reads, "GB"))
    Int disk_size = if proposed_disk < min_disk then min_disk else proposed_disk

    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             memory,
        disk_gb:            disk_size,
        preemptible_tries:  0,  # don't bother with premption on GCP with resource request + the expected long run, you'll end up waiting longer & paying more because GCP will prempt you anyway
        max_retries:        1,
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
    input {
        File reads
        String prefix = "out"
        String zones

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        export MONITOR_MOUNT_POINT="/cromwell_root/"
        bash /opt/vm_local_monitoring_script.sh &> resources.log &
        job_id=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')

        time \
        hifiasm \
            -o ~{prefix} \
            -t~{num_cpus} \
            --primary \
            ~{reads} \
            2>&1 | tee hifiasm.log

        if ps -p "${job_id}" > /dev/null; then kill "${job_id}"; fi
        tree -h .

        # tricky, outputs generated this way has no "bp" in their file names
        # GFA graph to contigs, primary
        awk '/^S/{print ">"$2; print $3}' \
            ~{prefix}.p_ctg.gfa \
        > ~{prefix}.p_ctg.fa

        # GFA graph to contigs, alternate
        awk '/^S/{print ">"$2; print $3}' \
            ~{prefix}.a_ctg.gfa \
        > ~{prefix}.a_ctg.fa

        for ff in ./*_ctg.fa; do
            bgzip -@2 -k --index "${ff}"
        done
    >>>

    output {
        File resouce_monitor_log = "resources.log"
        File log = "hifiasm.log"

        File primary_gfa  = "~{prefix}.p_ctg.gfa"
        File primary_tigs = "~{prefix}.p_ctg.fa.gz"
        File primary_fa_gzi = "~{prefix}.p_ctg.fa.gz.gzi"

        File alternate_gfa  = "~{prefix}.a_ctg.gfa"
        File alternate_tigs = "~{prefix}.a_ctg.fa.gz"
        File alternate_tigs_gzi = "~{prefix}.a_ctg.fa.gz.gzi"
    }

    #########################
    Int proposed_memory = 4 * ceil(size(reads, "GB"))
    Int memory = if proposed_memory < 96 then 96 else if proposed_memory > 512 then 512 else proposed_memory # this 96 magic number is purely empirical
    Int n = memory / 4  # this might be an odd number
    Int num_cpus_proposal = if (n/2)*2 == n then n else n+1  # a hack because WDL doesn't have modulus operator
    Int num_cpus = if num_cpus_proposal > 96 then 96 else num_cpus_proposal

    Int min_disk = 75
    Int proposed_disk = 5 * ceil(size(reads, "GB"))
    Int disk_size = if proposed_disk < min_disk then min_disk else proposed_disk

    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             memory,
        disk_gb:            disk_size,
        preemptible_tries:  0, # don't bother with premption on GCP with resource request + the expected long run, you'll end up waiting longer & paying more because GCP will prempt you anyway
        max_retries:        1,
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
