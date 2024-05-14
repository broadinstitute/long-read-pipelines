version 1.0

import "../../../structs/Structs.wdl"

workflow CoverageFromMosdepthFiles {

    input {
        File ref_fasta_fai
        Array[File] contig_files
        Array[String]? contigs
        Int cov_threshold = 5
    }

    call GetChrsFromFai { input: ref_fasta_fai = ref_fasta_fai }
    Array[String] select_contigs = select_first([contigs, GetChrsFromFai.chrs])

    call MosDepthGenomeCoverage {
        input:
            ref_fasta_fai = ref_fasta_fai,
            contigs = select_contigs,
            contig_files = contig_files,
            cov_threshold = cov_threshold
    }

    output {
        Float genome_cvg_at_thresh = MosDepthGenomeCoverage.genome_pct_cov_at_threshold
    }
}

task MosDepthGenomeCoverage {
    meta {
        description: "Takes coverage information from a list of Mosdepth files to give genome-wide coverage metrics"
    }
    
    parameter_meta {
        ref_fasta_fai: "Index for reference genome"
        contigs: "The contigs you are interested in to calculate genome-wide coverage metrics"
        contig_files: "Files to consider when calculating genome-wide coverage"
        cov_threshold: "Coverage threshold to report"
    }

    input {
        File ref_fasta_fai
        Array[String] contigs
        Array[File] contig_files
        Int cov_threshold = 5

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        cut -f1,2 ~{ref_fasta_fai} > all_contigs_w_lens.txt
        
        python3 <<CODE
        import numpy as np
        import os

        contigs_to_use = np.loadtxt("~{write_lines(contigs)}", dtype='str').tolist()
        files_to_use = np.loadtxt("~{write_lines(contig_files)}", dtype='str').tolist()
        
        # Get lengths of every contig
        contig_lens = {}
        with open("all_contigs_w_lens.txt", 'r') as f:
            for line in f:
                contig_name = line.split('\t')[0]
                contig_len = int(line.split('\t')[1])
                contig_lens[contig_name] = contig_len

        total_length = 0
        total_weighted_coverage = 0
        contig_cov_files = {}

        # Initialize contig_cov_files dictionary
        filepaths = []
        for filepath in files_to_use:
            if (any(c in filepath for c in contigs_to_use)):
                for contig in contigs_to_use:
                    if contig in filepath:
                        contig_cov_files[contig] = filepath
        
        # Loop through all contigs
        for contig, contig_file in contig_cov_files.items():
            contig_cvg = -1.0
            with open(contig_file) as f:
                # The Mosdepth file should have a three-column tsv of the following order: contig_name, coverage, percent_at_coverage
                for line in f:
                    if (contig in line) and ("\t~{cov_threshold}\t" in line):
                        contig_cvg = float(line.split('\t')[2])
                        break
            # if contig_cvg is still -1 (i.e. in the file, there was no line that satisfied this cov_threshold), set contig_cvg to 0
            if contig_cvg < 0:
                contig_cvg = 0.0

            # Rolling weighted average
            contig_len = contig_lens[contig]
            total_weighted_coverage += contig_len * contig_cvg
            total_length += contig_len
        
        # Write it into a file
        with open("genome_wide_coverage.txt", "w") as file:
            file.write(f"{total_weighted_coverage / total_length}")
        
        CODE
    >>>

    output {
        Float genome_pct_cov_at_threshold = read_float("genome_wide_coverage.txt")
    }

    #######################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-10x:0.1.18"
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

task GetChrsFromFai {
    input {
        File ref_fasta_fai

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        cut -f1 ~{ref_fasta_fai}
    >>>

    output {
        Array[String] chrs = read_lines(stdout())
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-10x:0.1.18"
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