version 1.0

##########################################################################################
## A workflow that performs trio-binning of child long reads given parental (short) reads.
## Based on the trio-canu publication
##    De novo assembly of haplotype-resolved genomes with trio binning
##    https://www.nature.com/articles/nbt.4277
## We divide the workflow into two parts
##   part one: collect k-mer stats given parental (short) reads
##   part two: given the k-mer stats database from part one, classify child long reads
##########################################################################################

import "../../../tasks/Preprocessing/CollectParentsKmerStats.wdl" as stats

import "../../../tasks/Utility/AssignChildLongReads.wdl" as asign

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

    parameter_meta {
        workdir_name:                  "name of working directory"
        genome_size:                   "an esimate on genome size of the specicies (affects k-value picking)"
        kmerSize:                      "[optional] force specifying k-value in collecting k-mer stats on parents"
        father_short_reads_bucket:     "GCS bucket path holding FASTA/FASTQ of (short) reads of paternal origin"
        mother_short_reads_bucket:     "GCS bucket path holding FASTA/FASTQ of (short) reads of maternal origin"
        child_long_reads_bucket:       "GCS bucket path holding FASTA/FASTQ of child long reads"
        long_read_platform:            "platform of long read sequencing; currently only one of [pacbio-raw, nanopore-raw] is supported"
        vm_local_monitoring_script:    "GCS file holding a resouce monitoring script that runs locally and collects info for a very specific purpose"
        meryl_operations_threads_est:  "[default-valued] estimate on how many threads to allocate to k-mer stats collection step"
        child_read_assign_threads_est: "[default-valued] estimate on how many threads to allocate to the child longread classification step"
        child_read_assign_memoryG_est: "[default-valued] estimate on how many GB memory to allocate to the child longread classification step"
        run_with_debug:                "[optional] whether to run in debug mode (takes significantly more disk space and more logs); defaults to false"
    }

    ############################################################################
    # we divide the workflow into two parts
    # part one: collect k-mer stats given parental (short) reads
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

    # part two: given the k-mer stats database from part one, classify child long reads
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

    # most important outputs are classified child longreads
    output {
        File reads_assigned_to_father = AssignChildLongReads.reads_assigned_to_father
        File reads_assigned_to_mother = AssignChildLongReads.reads_assigned_to_mother
        File unassigned_reads = AssignChildLongReads.unassigned_reads
    }
}
