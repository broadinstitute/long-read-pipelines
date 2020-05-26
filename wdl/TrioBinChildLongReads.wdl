version 1.0

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.8/wdl/tasks/CollectParentsKmerStats.wdl" as stats

import "https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/lrp_2.1.8/wdl/tasks/AssignChildLongReads.wdl" as asign

# A workflow that performs triobinning of child long reads given parental short reads
workflow TrioBinChildLongReads {
    input{

        String workdir_name

        String genome_size
        Int? kmerSize

        String father_short_reads_bucket
        String mother_short_reads_bucket

        String child_long_reads_bucket
        # currently the following only takes one of [pacbio-raw, nanopore-raw]
        String long_read_platform

        File vm_local_monitoring_script

        # these numbers will be used to request VMs on which all meryl jobs are run, in stages
        Int meryl_operations_threads_est = 8

        Int child_read_assign_threads_est = 36
        Int child_read_assign_memoryG_est = 32

        Boolean? run_with_debug = false
    }

    call stats.CollectParentsKmerStats as CollectParentsKmerStats {

        input:
            workdir_name = workdir_name,
            genome_size = genome_size,
            kmerSize = kmerSize,
            father_short_reads_bucket = father_short_reads_bucket,
            mother_short_reads_bucket = mother_short_reads_bucket,
            meryl_operations_threads_est = meryl_operations_threads_est,
            run_with_debug = run_with_debug
    }

    call asign.AssignChildLongReads as AssignChildLongReads {
        input:
            workdir_name = workdir_name,

            meryl_db_files_father = CollectParentsKmerStats.Father_haplotype_merylDB,
            meryl_db_files_mother = CollectParentsKmerStats.Mother_haplotype_merylDB,
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

