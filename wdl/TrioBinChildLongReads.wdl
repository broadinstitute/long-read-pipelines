version 1.0

import "CollectParentsKmerStats.wdl" as stats

# A workflow that performs triobinning of child long reads given parental short reads
workflow TrioBinChildLongReads {
    input{

        String workdir_name

        Int? kmerSize

        String father_short_reads_bucket
        String mother_short_reads_bucket

        String child_long_reads_bucket
        # currently the following only takes one of [pacbio-raw, nanopore-raw]
        String long_read_platform

        File vm_local_monitoring_script

        # these numbers will be used to request VMs on which all meryl jobs are run, in stages
        Int meryl_operations_threads_est = 8

        Int child_read_assign_threads_est = 64
        Int child_read_assign_memoryG_est = 64

        Boolean? run_with_debug = false
    }

    call stats.CollectParentsKmerStats as CollectParentsKmerStats {

        input:
            workdir_name = workdir_name,
            kmerSize = kmerSize,
            father_short_reads_bucket = father_short_reads_bucket,
            mother_short_reads_bucket = mother_short_reads_bucket,
            meryl_operations_threads_est = meryl_operations_threads_est,
            run_with_debug = run_with_debug
    }

    call AssignChildLongReads {
        input:
            workdir_name = workdir_name,

            meryl_subtract_father = CollectParentsKmerStats.Father_haplotype_merylDB,
            meryl_subtract_mother = CollectParentsKmerStats.Mother_haplotype_merylDB,
            meryl_stats_father = CollectParentsKmerStats.Father_reads_statistics,
            meryl_stats_mother = CollectParentsKmerStats.Mother_reads_statistics,

            child_long_reads_bucket = child_long_reads_bucket,
            long_read_platform = long_read_platform,

            child_read_assign_threads_est = child_read_assign_threads_est,
            child_read_assign_memoryG_est = child_read_assign_memoryG_est,
            vm_local_monitoring_script = vm_local_monitoring_script,
            run_with_debug = run_with_debug
    }

    output {
        File reads_assigned_to_father = AssignChildLongReads.reads_assigned_to_father
        File reads_assigned_to_mother = AssignChildLongReads.reads_assigned_to_mother
        File unassigned_reads = AssignChildLongReads.unassigned_reads
    }
}

###############################################################

# actually assign child long reads
task AssignChildLongReads {
    input{

        String workdir_name

        String child_long_reads_bucket
        # currently the following only takes one of [pacbio-raw, nanopore-raw]
        String long_read_platform

        Array[File] meryl_subtract_father
        Array[File] meryl_subtract_mother

        File meryl_stats_father
        File meryl_stats_mother

        Int child_read_assign_threads_est
        Int child_read_assign_memoryG_est

        File vm_local_monitoring_script

        Boolean run_with_debug = false

        RuntimeAttr? runtime_attr_override
    }

    String extra_args = if (run_with_debug) then "-debug" else " "
    String resource_script_name = basename(vm_local_monitoring_script)

    command <<<
        set -euo pipefail

        ##############################
        # parallel localize the input reads (remove trailing slash first to be safe)
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        child_path=$(echo ~{child_long_reads_bucket} | sed 's:/*$::')
        a=$(gsutil ls "${child_path}/"*.fastq.gz | wc -l)
        if [[ $a == 0 ]]; then
          echo "no reads in ~{child_long_reads_bucket}" && exit 1
        fi
        echo "==================================="
        echo "BEGIN LOCALIZING CHILD LONG READS"
        date -u
        mkdir child && gsutil -mq cp "${child_path}/"*.fastq.gz child/ && echo "Child localized"
        date -u
        echo "DONE LOCALIZING CHILD LONG READS"
        echo "==================================="

        ##########
        # prep files from previous stages
        mkdir -p workdir/canu-logs workdir/canu-scripts workdir/haplotype/0-kmers/
        mkdir -p workdir/haplotype/0-kmers/haplotype-Father.meryl \
                 workdir/haplotype/0-kmers/haplotype-Mother.meryl
        for ff in `ls ~{sep=' ' meryl_subtract_father}`; do
            mv $ff workdir/haplotype/0-kmers/haplotype-Father.meryl/
        done
        for ff in `ls ~{sep=' ' meryl_subtract_mother}`; do
            mv $ff workdir/haplotype/0-kmers/haplotype-Mother.meryl/
        done

        mv ~{meryl_stats_father} workdir/haplotype/0-kmers/
        mv ~{meryl_stats_mother} workdir/haplotype/0-kmers/

        # we need to have these success files to fool canu
        touch workdir/haplotype/0-kmers/meryl-count.success
        touch workdir/haplotype/0-kmers/meryl-merge.success
        touch workdir/haplotype/0-kmers/meryl-subtract.success

        ##########
        echo "==================================="
        date -u
        machine_mem=`echo $(($(getconf _PHYS_PAGES) * $(getconf PAGE_SIZE) / (1024 * 1024 * 1024)))`
        hap_mem=$(($machine_mem - 2))
        echo "Limiting hap memory to: ${hap_mem} G"
        export MONITOR_MOUNT_POINT="/cromwell_root"
        bash ~{vm_local_monitoring_script} &> resources.log &
        job_id=$(ps -aux | grep -F ~{resource_script_name} | head -1 | awk '{print $2}')
        canu \
            -haplotype \
            -p ~{workdir_name} \
            -d /cromwell_root/workdir/ \
            genomeSize=3.1G \
            beginConfigAt=hap \
            stopAfter=haplotype \
            hapMemory=${hap_mem} \
            -hapNames "Father Mother" \
            -~{long_read_platform} /cromwell_root/child/*.fastq.gz \
            ~{extra_args} ||
            cat workdir/haplotype/*.out
        kill $job_id || true
        du -sh workdir/haplotype/*
        date -u
        echo "==================================="
        ##########

        mv workdir/haplotype/haplotype.log .
        mv workdir/haplotype/splitHaplotype.000001.out .
        mv workdir/haplotype/splitHaplotype.sh .
        mv workdir/haplotype/haplotype-*.fasta.gz .

        # save logs and scripts
        tar -czf canu-logs.tar.gz workdir/canu-logs
        tar -czf canu-scripts.tar.gz workdir/canu-scripts
    >>>

    output {
        File logs = "canu-logs.tar.gz"
        File scripts = "canu-scripts.tar.gz"

        File resource_monitoring_log = "resources.log"

        File assignment_log = "haplotype.log"
        File assignment_job_log = "splitHaplotype.000001.out"
        File assignment_script = "splitHaplotype.sh"

        File reads_assigned_to_father = "haplotype-Father.fasta.gz"
        File reads_assigned_to_mother = "haplotype-Mother.fasta.gz"
        File unassigned_reads = "haplotype-unknown.fasta.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          child_read_assign_threads_est,
        mem_gb:             child_read_assign_memoryG_est,
        disk_gb:            300,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/canu:v1.9_wdl_patch_varibale_k"
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

