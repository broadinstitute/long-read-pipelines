version 1.0

import "Structs.wdl"

# A workflow that performs triobinning of child long reads given parental short reads
workflow TrioBinChildLongReads {
    input{

        String assembly_name

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

    call ParentalReadsRepartitionAndMerylConfigure {
        input:
            assembly_name = assembly_name,
            father_short_reads_bucket = father_short_reads_bucket,
            mother_short_reads_bucket = mother_short_reads_bucket,
            meryl_operations_threads_est = meryl_operations_threads_est,
            run_with_debug = run_with_debug
    }

    call PrintMerylMemory {
        input:
            meryl_memory_file = ParentalReadsRepartitionAndMerylConfigure.meryl_memory
    }

    scatter (pair in ParentalReadsRepartitionAndMerylConfigure.batch_id_and_parental_short_reads_tar) {
        call MerylCount {
            input:
                assembly_name = assembly_name,
                batch_id_hold_file = pair.left,
                parental_reads_for_this_batch = pair.right,
                meryl_count_script = ParentalReadsRepartitionAndMerylConfigure.count_script,
                meryl_operations_threads_est = meryl_operations_threads_est,
                meryl_memory_in_GB = PrintMerylMemory.meryl_memory_in_GB
        }
    }

    call MerylMergeAndSubtract {
        input:
            meryl_merge_script = ParentalReadsRepartitionAndMerylConfigure.merge_script,
            meryl_subtract_script = ParentalReadsRepartitionAndMerylConfigure.subtract_script,
            meryl_count_batches = MerylCount.count_output,
            meryl_operations_threads_est = meryl_operations_threads_est,
            meryl_memory_in_GB = PrintMerylMemory.meryl_memory_in_GB,
            run_with_debug = run_with_debug
    }

    call AssignChildLongReads {
        input:
            assembly_name = assembly_name,
            child_long_reads_bucket = child_long_reads_bucket,
            long_read_platform = long_read_platform,
            meryl_subtract_father = MerylMergeAndSubtract.subtract_output_father,
            meryl_subtract_mother = MerylMergeAndSubtract.subtract_output_mother,
            meryl_stats_father = MerylMergeAndSubtract.merge_stats_father,
            meryl_stats_mother = MerylMergeAndSubtract.merge_stats_mother,
            child_read_assign_threads_est = child_read_assign_threads_est,
            child_read_assign_memoryG_est = child_read_assign_memoryG_est,
            vm_local_monitoring_script = vm_local_monitoring_script,
            run_with_debug = run_with_debug
    }
}

###############################################################
# repartition the parental short reads, IO bound
task ParentalReadsRepartitionAndMerylConfigure {
    input{

        String assembly_name

        String father_short_reads_bucket
        String mother_short_reads_bucket

        Int meryl_operations_threads_est

        Boolean run_with_debug = false

        RuntimeAttr? runtime_attr_override
    }

    String extra_args = if (run_with_debug) then "-debug" else " "

    command <<<
        set -euo pipefail

        # parallel localize the input reads (remove trailing slash first to be safe)
        father_path=$(echo ~{father_short_reads_bucket} | sed 's:/*$::')
        mother_path=$(echo ~{mother_short_reads_bucket} | sed 's:/*$::')
        a=$(gsutil ls "${father_path}/"*.fastq.gz | wc -l)
        b=$(gsutil ls "${mother_path}/"*.fastq.gz | wc -l)
        if [[ $a == 0 ]]; then
          echo "no reads in ~{father_short_reads_bucket}" && exit 1
        elif [[ $b == 0 ]]; then
          echo "no reads in ~{mother_short_reads_bucket}" && exit 1
        fi
        echo "==================================="
        echo "BEGIN LOCALIZING PARENTAL READS"
        date -u
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        mkdir father && gsutil -mq cp "${father_path}/"*.fastq.gz father/ && echo "Father localized"
        mkdir mother && gsutil -mq cp "${mother_path}/"*.fastq.gz mother/ && echo "Mother localized"
        date -u
        echo "DONE LOCALIZING PARENTAL READS"
        echo "==================================="

        ##########
        # re-partition parental reads and
        # very importantly, we don't stop at that, we stop after parental kmer stats configure
        # (i.e. meryl-count, -merge, and -subtract) because
        # the custom docker allows us to output a batch-specific tar.gz of parental read files
        # but we don't use the configuration scripts in later stages
        echo "==================================="
        date -u
        canu \
            -haplotype \
            -p ~{assembly_name} \
            -d /cromwell_root/workdir/ \
            genomeSize=3.1G \
            stopAfter=parent-kmer-stat-conf \
            -haplotypeFather /cromwell_root/father/*.fastq.gz \
            -haplotypeMother /cromwell_root/mother/*.fastq.gz \
            ~{extra_args}
        date -u
        tree workdir # helps debugging, in case something went wrong
        # tar the repartitioned reads according to batches
        cd workdir/haplotype/0-kmers/
        for dict in *.dict; do
            new_name=$(echo $dict | sed 's/dict//')
            cat $dict | tar -czf $new_name"batch.tar.gz" -T - &
        done
        wait
        cd -
        echo "==================================="
        ##########

        # move configured shell scripts up for delocalization, then sed replace the thread configuration
        # recall that memory is set purely based on file count
        mv workdir/haplotype/0-kmers/meryl-count.memory .
        mv workdir/haplotype/0-kmers/*.sh .
        th_cnt=~{meryl_operations_threads_est}
        for script in *.sh; do
            sed -i -E "s/threads=[0-9]+/threads=$th_cnt/g" $script;
        done

        # move parental reads up for delocalization
        mv workdir/haplotype/0-kmers/*.dict .
        mv workdir/haplotype/0-kmers/*.batch.tar.gz .
        # grep and sort, then generate an array of flat files that just hold the "$hap-bathid"
        for batch_id in `grep -Eo "output=\"(Father|Mother)-[0-9]+\"" meryl-count.sh | awk -F '=' '{print $2}' | sed 's/"//g'`; do
            echo $batch_id > "$batch_id.txt"
        done
        ls *.txt

        # save logs and scripts
        tar -czf canu-logs.tar.gz workdir/canu-logs
        tar -czf canu-scripts.tar.gz workdir/canu-scripts
    >>>

    output {

        File logs = "canu-logs.tar.gz"
        File scripts = "canu-scripts.tar.gz"

        Array[File] batch_ids = glob("*.txt")
        Array[File] repartitioned_parental_reads_per_batch = glob("*.batch.tar.gz")
        Array[Pair[String, File]] batch_id_and_parental_short_reads_tar = zip(batch_ids, repartitioned_parental_reads_per_batch)

        File meryl_memory = "meryl-count.memory"
        File count_script = "meryl-count.sh"
        File merge_script = "meryl-merge.sh"
        File subtract_script = "meryl-subtract.sh"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            600,
        boot_disk_gb:       10,
        preemptible_tries:  0, # explicitly turn this off as we don't save that much for the disk, and pre-emption kills us
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/canu:v1.9_wdl_patch"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL" # LOCAL because this task is mostly IO operation
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task PrintMerylMemory {
    input {
        File meryl_memory_file
    }

    command <<<
        cat ~{meryl_memory_file}
    >>>

    output {
        Int meryl_memory_in_GB = read_int(stdout())
    }

    runtime {
        cpu: 1
        memory: "1 GiB"
        docker: "ubuntu:18.04"
    }
}

# mery-count on one batch
task MerylCount {
    input{

        String assembly_name

        File batch_id_hold_file
        File parental_reads_for_this_batch

        File meryl_count_script

        Int meryl_operations_threads_est
        Int meryl_memory_in_GB

        RuntimeAttr? runtime_attr_override
    }

    String postfix = basename(batch_id_hold_file, ".txt")

    Int half_meryl_memory = ceil(meryl_memory_in_GB/2)

    command <<<
        set -euo pipefail

        mkdir -p workdir/canu-logs workdir/canu-scripts
        mkdir -p workdir/haplotype/0-kmers/reads-Father workdir/haplotype/0-kmers/reads-Mother

        cp ~{meryl_count_script} workdir/haplotype/0-kmers/
        echo ~{meryl_memory_in_GB} > workdir/haplotype/0-kmers/meryl-count.memory

        echo "==================================="
        echo "BEGIN UNPACKING REPARTITIONED PARENTAL READS TO THE DESIRED LOCATIONS"
        date -u
        df -h
        tar xzfv ~{parental_reads_for_this_batch} -C workdir/haplotype/0-kmers/
        rm ~{parental_reads_for_this_batch} # save some disk space
        df -h
        date -u
        tree workdir
        echo "DONE UNPACKING REPARTITIONED PARENTAL READS TO THE DESIRED LOCATIONS"
        echo "==================================="

        ##########
        # run the script
        # this is essentially the command by canu::Execution::submitOrRunParallelJob
        echo "==================================="
        echo "BEGIN KMER COUNTING"
        date -u
        n=$(cat ~{batch_id_hold_file} | awk -F '-' '{print $2}' | sed 's/^0*//')
        log_name="meryl-count."~{postfix}".out"
        echo "Dealing with batch: ${n}, with log name: ${log_name}"
        cd workdir/haplotype/0-kmers/ && chmod +x meryl-count.sh
        ./meryl-count.sh ${n} > ${log_name} 2>&1 || cat ${log_name}
        tail -n 5 ${log_name}
        date -u
        tar -czf reads-~{postfix}.meryl.tar.gz reads-~{postfix}.meryl
        du -sh reads-~{postfix}.meryl.tar.gz
        cd -
        df -h
        date -u
        echo "DONE  KMER COUNTING"
        echo "==================================="
        ##########

        mv workdir/haplotype/0-kmers/reads-~{postfix}.meryl.tar.gz .
        mv workdir/haplotype/0-kmers/meryl-count.~{postfix}.out .

        # save logs and scripts
        tar -czf canu-logs.tar.gz workdir/canu-logs
        tar -czf canu-scripts.tar.gz workdir/canu-scripts
    >>>

    output {

        File logs = "canu-logs.tar.gz"
        File scripts = "canu-scripts.tar.gz"

        File count_log = "meryl-count.${postfix}.out"

        File count_output = "reads-~{postfix}.meryl.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          meryl_operations_threads_est / 2,
        mem_gb:             if (16 < half_meryl_memory) then half_meryl_memory else 16,
        disk_gb:            20,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/canu:v1.9_wdl_patch"
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

# merge the batches of the parental reads, generate one result per parent
# then subtract the k-mers
task MerylMergeAndSubtract {
    input{

        File meryl_merge_script
        File meryl_subtract_script

        Array[File] meryl_count_batches

        Int meryl_operations_threads_est
        Int meryl_memory_in_GB

        Boolean run_with_debug = false

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        echo "==================================="
        date -u
        echo "BEGIN UNPACKING SCATTERED meryl-count RESULTS"
        mkdir -p workdir/canu-logs workdir/canu-scripts workdir/haplotype/0-kmers/
        cp ~{meryl_merge_script} \
           ~{meryl_subtract_script} \
           workdir/haplotype/0-kmers/
        echo ~{meryl_memory_in_GB} > workdir/haplotype/0-kmers/meryl-count.memory
        for zipped_count in `ls ~{sep=' ' meryl_count_batches}`; do
            out_log_prefix=$(basename ${zipped_count})
            (tar xzfv ${zipped_count} -C workdir/haplotype/0-kmers/ > ${out_log_prefix}".unpacking.log" && rm -rf ${zipped_count}) &
        done
        wait
        echo "disk use after unpacking:"
        df -h
        date -u
        echo "DONE  UNPACKING SCATTERED meryl-count RESULTS"
        echo "==================================="

        ##########
        # run the shell scripts without canu
        echo "==================================="
        date -u
        cd workdir/haplotype/0-kmers/
        chmod +x meryl-merge.sh && chmod +x meryl-subtract.sh
        # merge
        ./meryl-merge.sh 1 > meryl-merge.000001.out 2>&1 &
        ./meryl-merge.sh 2 > meryl-merge.000002.out 2>&1 &
        wait
        if [[ ~{run_with_debug} == true ]]; then
            cat meryl-merge.000001.out
            cat meryl-merge.000002.out
        fi
        echo "disk use after merge operation:"
        df -h
        du -sh *
        rm -rf reads-Father-*.meryl.tar.gz # delete these count files after merge
        rm -rf reads-Mother-*.meryl.tar.gz # delete these count files after merge
        date -u
        echo "-----------------------------------"
        # subtract
        ./meryl-subtract.sh 1 > meryl-subtract.000001.out 2>&1 &
        ./meryl-subtract.sh 2 > meryl-subtract.000002.out 2>&1 &
        wait
        if [[ ~{run_with_debug} == true ]]; then
            cat meryl-subtract.000001.out
            cat meryl-subtract.000002.out
        fi
        echo "disk use after subtract operation:"
        df -h
        du -sh *
        date -u
        cd -
        echo "==================================="
        ##########

        mv workdir/haplotype/0-kmers/*.out .
        mv workdir/haplotype/0-kmers/reads-*.statistics .

        mv workdir/haplotype/0-kmers/haplotype-Father.meryl .
        mv workdir/haplotype/0-kmers/haplotype-Mother.meryl .

        # save logs and scripts
        tar -czf canu-logs.tar.gz workdir/canu-logs
        tar -czf canu-scripts.tar.gz workdir/canu-scripts
    >>>

    output {
        File logs = "canu-logs.tar.gz"
        File scripts = "canu-scripts.tar.gz"
        Array [File] unpacking_logs = glob("*.unpacking.log")

        Array[File] merge_log = glob("meryl-merge.*.out")
        File merge_stats_father = "reads-Father.statistics"
        File merge_stats_mother = "reads-Mother.statistics"
        # Array[File] merge_output_mother = glob("workdir/haplotype/0-kmers/reads-Mother.meryl/*") # we take these out for now because
        # Array[File] merge_output_father = glob("workdir/haplotype/0-kmers/reads-Father.meryl/*") # these are usually deleted in the canu pipeline

        Array[File] subtract_log = glob("meryl-subtract.*.out")
        Array[File] subtract_output_father = glob("haplotype-Father.meryl/*")
        Array[File] subtract_output_mother = glob("haplotype-Mother.meryl/*")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2 * meryl_operations_threads_est + 2,
        mem_gb:             2 * meryl_memory_in_GB + 2, # choosing this specification so that two parallel jobs can be executed at the same time
        disk_gb:            1500, # we strongly recommend you NOT lower this number, as we've seen it typically gets closer to 1.2T peak usage, and local SSD's increase by unit of 375GB
        boot_disk_gb:       10,
        preemptible_tries:  0, # explicitly turn off as this takes a long time
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/canu:v1.9_wdl_patch"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

# actually assign child long reads
task AssignChildLongReads {
    input{

        String assembly_name

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
            -p ~{assembly_name} \
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
        docker:             "quay.io/broad-long-read-pipelines/canu:v1.9_wdl_patch"
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

