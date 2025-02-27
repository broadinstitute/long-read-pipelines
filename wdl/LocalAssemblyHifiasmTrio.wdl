version 1.0

import "tasks/Utils.wdl" as Utils
import "tasks/Hifiasm.wdl" as Hifiasm
import "tasks/CallAssemblyVariants.wdl" as Align

workflow LocalAssembly {
    input {
        Array[String]+ loci
        String aligned_bam
        File   aligned_bai
        String prefix
#        Boolean add_unaligned_reads = false
        File maternal_bam
        File maternal_bai
        File paternal_bam
        File paternal_bai

        File ref_map_file
    }

    parameter_meta {
        loci:          "Loci to assemble. At least one is required. Reads from all loci will be merged for assembly. Format: [\"chr1:1000-2000\", \"chr1:5000-10000\"]"
        aligned_bam:   "aligned file"
        aligned_bai:   "index file"
        prefix:        "prefix for output files"

#        add_unaligned_reads: "set to true to include unaligned reads in the assembly (default: false)"

        ref_map_file:  "table indicating reference sequence and auxillary file locations"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    # get parental bam subsets
    scatter (locus in loci) {
        call Utils.SubsetBam as Paternal_subset {
            input:
                bam = paternal_bam,
                bai = paternal_bai,
                locus = locus
        }
    }

        scatter (locus in loci) {
        call Utils.SubsetBam as Maternal_subset {
            input:
                bam = maternal_bam,
                bai = maternal_bai,
                locus = locus
        }
    }

    # merge parental bam subsets
    if (length(loci) > 1)  {
        call Utils.MergeBams as MergePaternal {
            input:
                bams = Paternal_subset.subset_bam,
                prefix = "merged_paternal"
        }

        call Utils.MergeBams as MergeMaternal {
            input:
                bams = Maternal_subset.subset_bam,
                prefix = "merged_maternal"
        }
    }

    File paternal_subset_bam = select_first([MergePaternal.merged_bam, Paternal_subset.subset_bam[0]])
    File maternal_subset_bam = select_first([MergeMaternal.merged_bam, Maternal_subset.subset_bam[0]])

    # convert parental bam subsets to fastq
    call Utils.BamToFastq as PaternalFastq {
        input:
            bam = paternal_subset_bam,
            prefix = "paternal"
    }

    call Utils.BamToFastq as MaternalFastq {
        input:
            bam = maternal_subset_bam,
            prefix = "maternal"
    }

    # count kmers from parents
    call CountKmers as PaternalKmers {
        input:
            fq = PaternalFastq.reads_fq,
            prefix = "paternal"
    }

    call CountKmers as MaternalKmers {
        input:
            fq = MaternalFastq.reads_fq,
            prefix = "maternal"
    }

    # get child bam subset
    scatter (locus in loci) {
        call Utils.SubsetBam {
            input:
                bam = aligned_bam,
                bai = aligned_bai,
                locus = locus
        }
    }

    if (length(loci) > 1)  {
        call Utils.MergeBams {
            input:
                bams = SubsetBam.subset_bam,
                prefix = "merged"
        }
    }

    File subset_bam = select_first([MergeBams.merged_bam, SubsetBam.subset_bam[0]])

    call Utils.BamToFastq {
        input:
            bam = subset_bam,
            prefix = prefix
    }

    call Assemble_trio {
        input:
            reads = BamToFastq.reads_fq,
            prefix = prefix,
            mat_yak = MaternalKmers.yak,
            pat_yak = PaternalKmers.yak
    }

    call Align.AlignAsPAF as align_hap1 {
        input:
            asm_fasta = Assemble_trio.h1_fa,
            ref_fasta = ref_map['fasta'],
            prefix = "~{prefix}_h1"
    }

        call Align.AlignAsPAF as align_hap2 {
        input:
            asm_fasta = Assemble_trio.h2_fa,
            ref_fasta = ref_map['fasta'],
            prefix = "~{prefix}_h2"
    }

    output {
        File h1_fa = Assemble_trio.h1_fa
        File h2_fa = Assemble_trio.h2_fa
        File h1_aligned = align_hap1.paf
        File h2_aligned = align_hap2.paf
    }
}

task CountKmers {
    input {
        File fq
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 2,
        disk_gb: 5*ceil(size(fq, "GB")) + 20,
        cpu_cores: 16,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10,
        docker: "us.gcr.io/broad-dsp-lrma/lr-hifiasm:0.20.0"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         runtime_default.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      runtime_default.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       runtime_default.max_retries])
        docker:                 select_first([runtime_attr.docker,            runtime_default.docker])
    }

    command <<<
        set -euxo pipefail

        yak count -k31 -b37 -t16 -o ~{prefix}.yak ~{fq}
    >>>

    output {
        File yak = "~{prefix}.yak"
    }
}

task Assemble_trio {
    input {
        File reads
        String prefix = "out"
        File mat_yak
        File pat_yak
        Int num_cpus = 32

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        reads:    "reads (in fasta or fastq format, compressed or uncompressed)"
        prefix:   "prefix to apply to assembly output filenames"
        num_cpus: "number of CPUs to parallelize over"
        mat_yak:  "maternal yak file"
        pat_yak:  "paternal yak file"
    }

    Int disk_size = ceil(size(mat_yak, "GB")) + ceil(size(pat_yak, "GB")) + 10 * ceil(size(reads, "GB"))

    command <<<
        set -euxo pipefail

        hifiasm -o ~{prefix} -t~{num_cpus} -1 ~{pat_yak} -2 ~{mat_yak} ~{reads}
        ls
        awk '/^S/{print ">"$2; print $3}' ~{prefix}.hap1.p_ctg.gfa > ~{prefix}.hap1.p_ctg.fa
        awk '/^S/{print ">"$2; print $3}' ~{prefix}.hap2.p_ctg.gfa > ~{prefix}.hap2.p_ctg.fa
    >>>

    output {
        File h1_fa = "~{prefix}.hap1.p_ctg.fa"
        File h2_fa = "~{prefix}.hap2.p_ctg.fa"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             150,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-hifiasm:0.20.0"
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
