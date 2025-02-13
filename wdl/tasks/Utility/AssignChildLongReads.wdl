version 1.0

import "../../structs/Structs.wdl"

workflow AssignChildLongReadsGivenParentalKmerStats {

    meta {
        description: "A workflow that performs trio-binning of child long reads given parental (short) reads. Based on the trio-canu publication 'De novo assembly of haplotype-resolved genomes with trio binning' https://www.nature.com/articles/nbt.4277 . This holds the sub-workflow for part two: given the k-mer stats database from part one, classify child long reads.  We separate this out based on two concerns:  1. we can test out using different k-value when collecting parental k-mer states  2. we can collect parental k-mer stats once and classify all children reads (different sibblings, technologies) separately."
    }

    parameter_meta {
        workdir_name:            "name of working directory"
        child_long_reads_bucket: "GCS bucket path holding FASTA/FASTQ of child long reads"
        long_read_platform:      "platform of long read sequencing; currently only one of [pacbio-raw, nanopore-raw] is supported"

		meryl_db_files_father: "Meryl databases files on paternal (short) reads"
		meryl_db_files_mother: "Meryl databases files on maternal (short) reads"
		meryl_stats_father:    "Meryl statistics single file on paternal (short) reads"
		meryl_stats_mother:    "Meryl statistics single file on maternal (short) reads"

        vm_local_monitoring_script:    "GCS file holding a resouce monitoring script that runs locally and collects info for a very specific purpose"
        meryl_operations_threads_est:  "[default-valued] estimate on how many threads to allocate to k-mer stats collection step"
        child_read_assign_threads_est: "[default-valued] estimate on how many threads to allocate to the child longread classification step"
        child_read_assign_memoryG_est: "[default-valued] estimate on how many GB memory to allocate to the child longread classification step"
        run_with_debug:                "[optional] whether to run in debug mode (takes significantly more disk space and more logs); defaults to false"
    }

	input{

        String workdir_name

        String child_long_reads_bucket

        String long_read_platform

        Array[File] meryl_db_files_father
        Array[File] meryl_db_files_mother

        File meryl_stats_father
        File meryl_stats_mother

        File vm_local_monitoring_script

        Int child_read_assign_threads_est = 36
        Int child_read_assign_memoryG_est = 32

        Boolean? run_with_debug = false
    }



	call AssignChildLongReads {
        input:
            workdir_name = workdir_name,

            meryl_db_files_father = meryl_db_files_father,
            meryl_db_files_mother = meryl_db_files_mother,
            meryl_stats_father = meryl_stats_father,
            meryl_stats_mother = meryl_stats_mother,

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

        Array[File] meryl_db_files_father
        Array[File] meryl_db_files_mother

        File meryl_stats_father
        File meryl_stats_mother

        Int child_read_assign_threads_est
        Int child_read_assign_memoryG_est

        File vm_local_monitoring_script

        Boolean? run_with_debug = false

        RuntimeAttr? runtime_attr_override
    }

    String extra_args = if (defined(run_with_debug) && run_with_debug) then "-debug" else " "
    String resource_script_name = basename(vm_local_monitoring_script)

    command <<<
        set -euo pipefail

        ##############################
        # parallel localize the input reads (remove trailing slash first to be safe)
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        child_path=$(echo ~{child_long_reads_bucket} | sed 's:/*$::')
        if ! gsutil ls "${child_path}/" | grep -Eq "(fa|fq|fasta|fastq)(.gz)?$"; then
          echo "no reads in ~{child_long_reads_bucket}" && exit 1
        fi
        echo "==================================="
        echo "BEGIN LOCALIZING CHILD LONG READS"
        date -u
        mkdir child
        gsutil ls "${child_path}/" | grep -E "(fa|fq|fasta|fastq)(.gz)?$" | gsutil -mq cp -I child/
        echo "Child localized"
        date -u
        echo "DONE LOCALIZING CHILD LONG READS"
        echo "==================================="
        read_input_format=$(gsutil ls "${child_path}/" | grep -Eo "(fa|fq|fasta|fastq)(.gz)?$" | uniq)

        ##########
        # prep files from previous stages
        mkdir -p workdir/canu-logs workdir/canu-scripts workdir/haplotype/0-kmers/
        mkdir -p workdir/haplotype/0-kmers/haplotype-Father.meryl \
                 workdir/haplotype/0-kmers/haplotype-Mother.meryl
        for ff in `ls ~{sep=' ' meryl_db_files_father}`; do
            mv $ff workdir/haplotype/0-kmers/haplotype-Father.meryl/
        done
        for ff in `ls ~{sep=' ' meryl_db_files_mother}`; do
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
            -~{long_read_platform} /cromwell_root/child/*${read_input_format} \
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
        disk_gb:            500,
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
