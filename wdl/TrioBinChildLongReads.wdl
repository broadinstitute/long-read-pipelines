version 1.0

import "Structs.wdl"

# A workflow that performs triobinning of child long reads given parental short reads
workflow TrioBinChildLongReads {
    input{

        String assembly_name

        String father_short_reads_bucket
        String mother_short_reads_bucket

        String child_long_reads_bucket

        # these numbers will be used to request VMs on which and run all meryl jobs are run, in stages, including the configuration job
        Int meryl_operations_threads_est = 8
        Int meryl_operations_memoryG_est = 32

        Int child_read_assign_threads_est = 40
        Int child_read_assign_memoryG_est = 256

        Boolean? run_with_debug = false
    }

    call ParentalReadsRepartition {
        input:
            assembly_name = assembly_name,
            father_short_reads_bucket = father_short_reads_bucket,
            mother_short_reads_bucket = mother_short_reads_bucket,
            run_with_debug = run_with_debug
    }

    call MerylConfigure {
        input:
            assembly_name = assembly_name,
            father_and_mother_read_files_count_formated = ParentalReadsRepartition.father_and_mother_read_files_count_formated,
            meryl_operations_threads_est = meryl_operations_threads_est,
            meryl_operations_memoryG_est = meryl_operations_memoryG_est
    }

    scatter (batchIdFile in MerylConfigure.batch_ids) {
        call MerylCount {
            input:
                assembly_name = assembly_name,
                batch_id_hold_file = batchIdFile,
                meryl_count_script = MerylConfigure.count_script,
                meryl_count_memory = MerylConfigure.meryl_memory,
                repartitioned_father_reads = ParentalReadsRepartition.repartitioned_father_reads,
                repartitioned_mother_reads = ParentalReadsRepartition.repartitioned_mother_reads,
                meryl_operations_threads_est = meryl_operations_threads_est,
                meryl_operations_memoryG_est = meryl_operations_memoryG_est
        }
    }

    call MerylMergeAndSubtract {
        input:
            meryl_merge_script = MerylConfigure.merge_script,
            meryl_subtract_script = MerylConfigure.subtract_script,
            meryl_count_memory = MerylConfigure.meryl_memory,
            meryl_count_batches = MerylCount.count_output,
            meryl_operations_threads_est = meryl_operations_threads_est,
            meryl_operations_memoryG_est = meryl_operations_memoryG_est,
            run_with_debug = run_with_debug
    }

    call AssignChildLongReads {
        input:
            assembly_name = assembly_name,
            child_long_reads_bucket = child_long_reads_bucket,
            meryl_subtract_father = MerylMergeAndSubtract.subtract_output_father,
            meryl_subtract_mother = MerylMergeAndSubtract.subtract_output_mother,
            meryl_stats_father = MerylMergeAndSubtract.merge_stats_father,
            meryl_stats_mother = MerylMergeAndSubtract.merge_stats_mother,
            child_read_assign_threads_est = child_read_assign_threads_est,
            child_read_assign_memoryG_est = child_read_assign_memoryG_est,
            run_with_debug = run_with_debug
    }
}

###############################################################
# TODO: serial localization is done with Cromwell 45, slowing down the pipeline signnificantly.
#       Try and see if Cromwell 47 helps.
###############################################################
# repartition the parental short reads, IO bound
task ParentalReadsRepartition {
    input{

        String assembly_name

        String father_short_reads_bucket
        String mother_short_reads_bucket

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
          echo "no reads in ~{father_short_reads_bucket}"
        elif [[ $b == 0 ]]; then
          echo "no reads in ~{mother_short_reads_bucket}"
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
        # re-partition parental reads and stop there
        echo "==================================="
        date -u
        canu \
            -haplotype \
            -p ~{assembly_name} \
            -d /cromwell_root/workdir/ \
            genomeSize=3.1G \
            stopAfter=parental-reads-repartition \
            -haplotypeFather /cromwell_root/father/*.fastq.gz \
            -haplotypeMother /cromwell_root/mother/*.fastq.gz \
            ~{extra_args}
        date -u
        tree workdir # helps debugging, in case something went wrong
        echo "==================================="
        ##########

        # save logs and scripts
        tar -czf canu-logs.tar.gz workdir/canu-logs
        tar -czf canu-scripts.tar.gz workdir/canu-scripts

        # pack up re-partitioned reads
        tar -czf reads-Father.tar.gz workdir/haplotype/0-kmers/reads-Father &
        tar -czf reads-Mother.tar.gz workdir/haplotype/0-kmers/reads-Mother &
        wait
        du -sh father mother reads-Father.tar.gz reads-Mother.tar.gz

        # a trick to be used later
        a=$(ls workdir/haplotype/0-kmers/reads-Father/*.fasta.gz | wc -l | awk '{print $1}')
        b=$(ls workdir/haplotype/0-kmers/reads-Mother/*.fasta.gz | wc -l | awk '{print $1}')
        echo $a"_"$b > parental_read_files_count.txt
    >>>

    output {

        File logs = "canu-logs.tar.gz"
        File scripts = "canu-scripts.tar.gz"

        File repartitioned_father_reads = "reads-Father.tar.gz"
        File repartitioned_mother_reads = "reads-Mother.tar.gz"

        File father_and_mother_read_files_count_formated = "parental_read_files_count.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            600,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/canu:v1.9_wdl"
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

# configures the meryl-count, -merge, and -subtract shell scripts
# a good estimate on the machine configuration is critial for good performance on count/merge/subtract operations
task MerylConfigure {
    input{

        String assembly_name

        # this is a specially coded string of format "a_b", where
        # a is the number of re-partitioned file count for father, and b is for mother
        # we are going to do a little trick below
        File father_and_mother_read_files_count_formated

        Int meryl_operations_threads_est
        Int meryl_operations_memoryG_est

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        mkdir -p workdir/canu-logs workdir/canu-scripts
        mkdir -p workdir/haplotype/0-kmers/reads-Father workdir/haplotype/0-kmers/reads-Mother

        # here we do a trick, because we know that the configuration subroutine only takes a look at how many fasta files there are
        # instead of how large they are (they are assumed to be of fixed size--the reason for repartition)
        # we are simply going to just touch-create those files, to "fool" the configuration subroutine
        father_count=$(cat ~{father_and_mother_read_files_count_formated} | awk -F '_' '{print $1}')
        mother_count=$(cat ~{father_and_mother_read_files_count_formated} | awk -F '_' '{print $2}')
        echo "Father read files count: $father_count"
        echo "Mother read files count: $mother_count"
        for i in `seq $father_count`; do
            if [[ $i -le 999 ]]; then
                padded=$(printf %03d $i) # necessary, as the canu code format things in that order
            else
                padded=$(printf %04d $i)
            fi
            touch workdir/haplotype/0-kmers/reads-Father/reads-Father-$padded.fasta.gz;
        done
        for i in `seq $mother_count`; do
            if [[ $i -le 999 ]]; then
                padded=$(printf %03d $i) # necessary, as the canu code format things in that order
            else
                padded=$(printf %04d $i)
            fi
            touch workdir/haplotype/0-kmers/reads-Mother/reads-Mother-$padded.fasta.gz;
        done
        # this avoids re-doing the repartition work exactly because
        # accompanying the repartitioned reads, there are these "*.success" files
        # fooling canu again
        touch workdir/haplotype/0-kmers/reads-Father/reads-Father.success
        touch workdir/haplotype/0-kmers/reads-Mother/reads-Mother.success

        ##########
        echo "==================================="
        date -u
        canu \
            -haplotype \
            -p ~{assembly_name} \
            -d /cromwell_root/workdir/ \
            genomeSize=3.1G \
            beginConfigAt=meryl \
            stopAfter=meryl-configure \
            -haplotypeFather /cromwell_root/workdir/haplotype/0-kmers/reads-Father/*.fasta.gz \
            -haplotypeMother /cromwell_root/workdir/haplotype/0-kmers/reads-Mother/*.fasta.gz \
            -debug # explicitly turn on debugging here because we want to see more of the configuration process, for easier debugging if something goes wrong
        date -u
        echo "==================================="
        ##########

        # save result
        cp workdir/haplotype/0-kmers/*.sh .
        cp workdir/haplotype/0-kmers/meryl-count.memory .

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

        File meryl_memory = "meryl-count.memory"
        File count_script = "meryl-count.sh"
        File merge_script = "meryl-merge.sh"
        File subtract_script = "meryl-subtract.sh"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          meryl_operations_threads_est + 2,
        mem_gb:             meryl_operations_memoryG_est + 2,
        disk_gb:            50, # such a small size because we don't actuall have any large files to deal with in this stage
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/canu:v1.9_wdl"
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

# mery-count on one batch
task MerylCount {
    input{

        String assembly_name

        File batch_id_hold_file
        File meryl_count_script
        File meryl_count_memory

        # tar.gz files, to avoid the super slow localization process
        File repartitioned_father_reads
        File repartitioned_mother_reads

        Int meryl_operations_threads_est
        Int meryl_operations_memoryG_est

        RuntimeAttr? runtime_attr_override
    }

    String postfix = basename(batch_id_hold_file, ".txt")

    command <<<
        set -euo pipefail

        mkdir -p workdir/canu-logs workdir/canu-scripts
        mkdir -p workdir/haplotype/0-kmers/reads-Father workdir/haplotype/0-kmers/reads-Mother

        cp ~{meryl_count_script} ~{meryl_count_memory} workdir/haplotype/0-kmers/

        # move reads to the correct location
        echo "==================================="
        echo "BEGIN UNPACKING REPARTITIONED PARENTAL READS TO THE DESIRED LOCATIONS"
        date -u
        tar xzf ~{repartitioned_father_reads} -C /cromwell_root/ &
        tar xzf ~{repartitioned_mother_reads} -C /cromwell_root/ &
        wait
        date -u
        echo "DONE UNPACKING REPARTITIONED PARENTAL READS TO THE DESIRED LOCATIONS"
        echo "==================================="

        ##########
        # run the script
        echo "==================================="
        date -u
        n=$(cat ~{batch_id_hold_file} | awk -F '-' '{print $2}' | sed 's/^0*//')
        log_name="meryl-count."~{postfix}".out"
        echo "Dealing with batch: ${n}, with log name: ${log_name}"
        cd workdir/haplotype/0-kmers/ && chmod +x meryl-count.sh # && sed -i '2iset -eu' meryl-count.sh # turn off this usual safeguard because the shell script has some non-best practice behavior
        ./meryl-count.sh ${n} > ${log_name} 2>&1 || cat ${log_name} # this is essentially the command by canu::Execution::submitOrRunParallelJob
        cd -
        date -u
        echo "==================================="
        ##########

        # rm to make space for tar
        rm -rf ~{repartitioned_father_reads} \
               ~{repartitioned_mother_reads}
        tar -czf reads-~{postfix}.meryl.tar.gz workdir/haplotype/0-kmers/reads-~{postfix}.meryl
        du -sh reads-~{postfix}.meryl.tar.gz

        # save logs and scripts
        tar -czf canu-logs.tar.gz workdir/canu-logs
        tar -czf canu-scripts.tar.gz workdir/canu-scripts
    >>>

    output {

        File logs = "canu-logs.tar.gz"
        File scripts = "canu-scripts.tar.gz"

        File count_log = "workdir/haplotype/0-kmers/meryl-count.${postfix}.out"

        File count_output = "reads-~{postfix}.meryl.tar.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          meryl_operations_threads_est + 2,
        mem_gb:             meryl_operations_memoryG_est + 2,
        disk_gb:            500,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/canu:v1.9_wdl"
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

# merge the batches of the parental reads, generate one result per parent
# then subtract the k-mers
task MerylMergeAndSubtract {
    input{

        File meryl_merge_script
        File meryl_subtract_script
        File meryl_count_memory

        Array[File] meryl_count_batches

        Int meryl_operations_threads_est
        Int meryl_operations_memoryG_est

        Boolean run_with_debug = false

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        echo "==================================="
        date -u
        echo "BEGIN UNPACKING SCATTERED meryl-count RESULTS"
        mkdir -p workdir/canu-logs workdir/canu-scripts workdir/haplotype/0-kmers/
        cp ~{meryl_merge_script} ~{meryl_subtract_script} ~{meryl_count_memory} workdir/haplotype/0-kmers/
        for zipped_count in `ls ~{sep=' ' meryl_count_batches}`; do
            tar xzf ${zipped_count} -C workdir/haplotype/0-kmers/ &
        done
        wait
        for zipped_count in `ls ~{sep=' ' meryl_count_batches}`; do
            rm -rf ${zipped_count}
        done
        date -u
        echo "DONE  UNPACKING SCATTERED meryl-count RESULTS"
        echo "==================================="

        ##########
        # run the shell scripts without canu
        echo "==================================="
        date -u
        cd workdir/haplotype/0-kmers/
        chmod +x meryl-merge.sh && chmod +x meryl-subtract.sh
        # sed -i '2iset -eu' meryl-merge.sh && sed -i '2iset -eu' meryl-subtract.sh && \ # turn off this usual safeguard because the shell script has some non-best practice behavior
        mv workdir/haplotype/0-kmers/reads-*.meryl . && rm -rf workdir
        ./meryl-merge.sh 1 > meryl-merge.000001.out 2>&1 &
        ./meryl-merge.sh 2 > meryl-merge.000002.out 2>&1 &
        wait
        if [[ ~{run_with_debug} == true ]]; then
            cat meryl-merge.000001.out
            cat meryl-merge.000002.out
        fi
        date -u
        ./meryl-subtract.sh 1 > meryl-subtract.000001.out 2>&1 &
        ./meryl-subtract.sh 2 > meryl-subtract.000002.out 2>&1 &
        wait
        if [[ ~{run_with_debug} == true ]]; then
            cat meryl-subtract.000001.out
            cat meryl-subtract.000002.out
        fi
        cd -
        date -u
        du -sh workdir
        echo "==================================="
        ##########

        # save logs and scripts
        tar -czf canu-logs.tar.gz workdir/canu-logs
        tar -czf canu-scripts.tar.gz workdir/canu-scripts
    >>>

    output {
        File logs = "canu-logs.tar.gz"
        File scripts = "canu-scripts.tar.gz"

        Array[File] merge_log = glob("workdir/haplotype/0-kmers/meryl-merge.*.out")
        File merge_stats_father = "workdir/haplotype/0-kmers/reads-Father.statistics"
        File merge_stats_mother = "workdir/haplotype/0-kmers/reads-Mother.statistics"
        # Array[File] merge_output_mother = glob("workdir/haplotype/0-kmers/reads-Mother.meryl/*") # we take these out for now because
        # Array[File] merge_output_father = glob("workdir/haplotype/0-kmers/reads-Father.meryl/*") # these are usually deleted in the canu pipeline

        Array[File] subtract_log = glob("workdir/haplotype/0-kmers/meryl-subtract.*.out")
        Array[File] subtract_output_father = glob("workdir/haplotype/0-kmers/haplotype-Father.meryl/*")
        Array[File] subtract_output_mother = glob("workdir/haplotype/0-kmers/haplotype-Mother.meryl/*")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2 * meryl_operations_threads_est + 2,
        mem_gb:             2 * meryl_operations_memoryG_est + 2,# choosing this specification so that two parallel jobs can be executed at the same time
        disk_gb:            1000,
        boot_disk_gb:       10,
        preemptible_tries:  0, # explicitly turn off as this takes a long time
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/canu:v1.9_wdl"
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

        Array[File] meryl_subtract_father
        Array[File] meryl_subtract_mother

        File meryl_stats_father
        File meryl_stats_mother

        Int child_read_assign_threads_est
        Int child_read_assign_memoryG_est

        Boolean run_with_debug = false

        RuntimeAttr? runtime_attr_override
    }

    String extra_args = if (run_with_debug) then "-debug" else " "

    command <<<
        set -euo pipefail

        ##############################
        # parallel localize the input reads (remove trailing slash first to be safe)
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        child_path=$(echo ~{child_long_reads_bucket} | sed 's:/*$::')
        a=$(gsutil ls "${child_path}/"*.fastq.gz | wc -l)
        if [[ $a == 0 ]]; then
          echo "no reads in ~{child_long_reads_bucket}"
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
        canu \
            -haplotype \
            -p ~{assembly_name} \
            -d /cromwell_root/workdir/ \
            genomeSize=3.1G \
            beginConfigAt=hap \
            stopAfter=haplotype \
            hapMemory=${hap_mem} \
            -hapNames "Father Mother" \
            -pacbio-raw /cromwell_root/child/*.fastq.gz \
            ~{extra_args} ||
            cat workdir/haplotype/*.out
        du -sh workdir/haplotype/*
        date -u
        echo "==================================="
        ##########

        # save logs and scripts
        tar -czf canu-logs.tar.gz workdir/canu-logs
        tar -czf canu-scripts.tar.gz workdir/canu-scripts
    >>>

    output {
        File logs = "canu-logs.tar.gz"
        File scripts = "canu-scripts.tar.gz"

        File assignment_log = "workdir/haplotype/haplotype.log"
        File assignment_job_log = "workdir/haplotype/splitHaplotype.000001.out"
        File assignment_script = "workdir/haplotype/splitHaplotype.sh"

        File reads_assigned_to_father = "workdir/haplotype/haplotype-Father.fasta.gz"
        File reads_assigned_to_mother = "workdir/haplotype/haplotype-Mother.fasta.gz"
        File unassigned_reads = "workdir/haplotype/haplotype-unknown.fasta.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          child_read_assign_threads_est,
        mem_gb:             child_read_assign_memoryG_est,
        disk_gb:            500,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/canu:v1.9_wdl"
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

