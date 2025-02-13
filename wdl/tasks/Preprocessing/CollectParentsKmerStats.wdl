version 1.0

import "../../structs/Structs.wdl"

workflow CollectParentsKmerStats {

    meta {
        description: "A workflow that performs trio-binning of child long reads given parental (short) reads. Based on the trio-canu publication https://www.nature.com/articles/nbt.4277. This holds the sub-workflow for part one: collect k-mer stats given parental (short) reads"
    }
    parameter_meta {
        workdir_name:                 "name of working directory"
        genome_size:                  "an esimate on genome size of the specicies (affects k-value picking)"
        kmerSize:                     "[optional] force specifying k-value in collecting k-mer stats on parents"

        father_short_reads_bucket:    "GCS bucket path holding FASTA/FASTQ of (short) reads of paternal origin"
        mother_short_reads_bucket:    "GCS bucket path holding FASTA/FASTQ of (short) reads of maternal origin"

        meryl_operations_threads_est: "[default-valued] estimate on how many threads to allocate to k-mer stats collection step"
        run_with_debug:               "[optional] whether to run in debug mode (takes significantly more disk space and more logs); defaults to false"
    }

    input{

        String workdir_name

        String genome_size
        Int? kmerSize

        String father_short_reads_bucket
        String mother_short_reads_bucket

        Int meryl_operations_threads_est = 8

        Boolean? run_with_debug = false
    }

    ############################################################################
    # we based the implementation of this workflow on a forked canu v1.9
    # the original canu 1.9 is not cloud-friendly
    call ParentalReadsRepartitionAndMerylConfigure {
        input:
            workdir_name = workdir_name,
            genome_size = genome_size,
            kmerSize = kmerSize,
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
                workdir_name = workdir_name,
                batch_id_hold_file = pair.left,
                parental_reads_for_this_batch = pair.right,
                meryl_count_script = ParentalReadsRepartitionAndMerylConfigure.count_script,
                kmerSize = kmerSize,
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

    output {
        Array[File] Father_haplotype_merylDB = MerylMergeAndSubtract.subtract_output_father
        Array[File] Mother_haplotype_merylDB = MerylMergeAndSubtract.subtract_output_mother

        File Father_reads_statistics = MerylMergeAndSubtract.merge_stats_father
        File Mother_reads_statistics = MerylMergeAndSubtract.merge_stats_mother
    }
}

###############################################################

# repartition the parental short reads for easier batch-processing by canu/meryl itself
# note that this step is IO bound
task ParentalReadsRepartitionAndMerylConfigure {
    input{

        String workdir_name

        String genome_size
        Int? kmerSize

        String father_short_reads_bucket
        String mother_short_reads_bucket

        Int meryl_operations_threads_est

        Boolean? run_with_debug = false

        RuntimeAttr? runtime_attr_override
    }

    String debug_option = if select_first([run_with_debug, false]) then "-debug" else " "
    String kmer_option = if (defined(kmerSize)) then ("-triobinK " + select_first([kmerSize])) else " "
    String extra_args = kmer_option + debug_option

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
            -p ~{workdir_name} \
            -d /cromwell_root/workdir/ \
            genomeSize=~{genome_size} \
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

        # move configured shell scripts up for delocalization
        mv workdir/haplotype/0-kmers/meryl-count.memory .
        mv workdir/haplotype/0-kmers/*.sh .
        # then sed replace the thread configuration:
        # 1. recall that memory was set purely based on file count hence no-need/better-not change
        # 2. the number of threads was configured above using threads available on this VM
        #    and we don't want that.
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
        boot_disk_gb:       25,
        preemptible_tries:  0, # explicitly turn this off as we don't save that much for the disk, and pre-emption kills us
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/canu:v1.9_wdl_patch_varibale_k"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD" # If this is too slow, revert to LOCAL (LOCAL because this task is mostly IO operation)
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

# a hackish step to simply print out memory configuration from the above step
# this value is used for configuring the VMs for actual batch Meryl count processing
#   as well as memory allocation for the merging and subtracting step
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

        String workdir_name

        File batch_id_hold_file
        File parental_reads_for_this_batch

        File meryl_count_script

        Int? kmerSize
        Int meryl_operations_threads_est
        Int meryl_memory_in_GB

        RuntimeAttr? runtime_attr_override
    }

    String postfix = basename(batch_id_hold_file, ".txt")

    Int emperical_memory_lower_limit = 18
    Int memory_to_use = if (meryl_memory_in_GB < emperical_memory_lower_limit) then emperical_memory_lower_limit else meryl_memory_in_GB

    Int emperical_thread_cout = 4 # based on monitor, the task is not CPU intensive (but memory intensive), and higher available CPU improved runtime only marginally

    Int disk_space_gb = if(defined(kmerSize)) then 100 else 50

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
        echo "----------"
        echo "tail log files"
        tail -n 5 ${log_name}
        echo "----------"
        date -u
        echo "Done counting, now compressing for delocalization..."
        tar --use-compress-program=pigz -cf reads-~{postfix}.meryl.tar.gz reads-~{postfix}.meryl
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
        cpu_cores:          emperical_thread_cout,
        mem_gb:             memory_to_use,
        disk_gb:            disk_space_gb,
        boot_disk_gb:       25,
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

# merge the batches of the parental reads, generate one result per parent
# then subtract the k-mers
task MerylMergeAndSubtract {
    input{

        File meryl_merge_script
        File meryl_subtract_script

        Array[File] meryl_count_batches

        Int meryl_operations_threads_est
        Int meryl_memory_in_GB

        Boolean? run_with_debug = false

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
        touch this.is.father.db && mv this.is.father.db haplotype-Father.meryl/
        touch this.is.mother.db && mv this.is.mother.db haplotype-Mother.meryl/

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
        cpu_cores:          2 * meryl_operations_threads_est + 6, # a bit more threads, a bit more concurrency for decompression at the beginning
        mem_gb:             3 * meryl_memory_in_GB + 6, # choosing this specification so that two parallel jobs can be executed at the same time
        disk_gb:            3000, # we strongly recommend you NOT change this number: 1) we've seen close to full disk peak usage, and 2) local SSD's increase by unit of 375GB, this is the maximum
        boot_disk_gb:       25,
        preemptible_tries:  0, # explicitly turn off as this takes a long time
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/canu:v1.9_wdl_patch_varibale_k"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD" # If this is too slow, revert to LOCAL
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
