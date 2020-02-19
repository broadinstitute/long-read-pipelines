version 1.0

import "Structs.wdl"
import "AssembleTarget.wdl" as ASM

workflow AssembleMHC {
    input {
        Array[File] mother_fastqs
        Array[File] father_fastqs

        File child_corrected_bam
        File child_corrected_bai
        File child_remaining_bam
        File child_remaining_bai

        String platform

        File ref_fasta
        File ref_dict
        File ref_fasta_fai
        File ref_fasta_amb
        File ref_fasta_ann
        File ref_fasta_bwt
        File ref_fasta_pac
        File ref_fasta_sa

        Array[String] mhc_region

        String prefix
    }

    call BwaMem as BwaMemMother {
        input:
            end1 = mother_fastqs[0],
            end2 = mother_fastqs[1],
            ref_fasta = ref_fasta,
            ref_dict = ref_dict,
            ref_fasta_fai =  ref_fasta_fai,
            ref_fasta_amb = ref_fasta_amb,
            ref_fasta_ann = ref_fasta_ann,
            ref_fasta_bwt = ref_fasta_bwt,
            ref_fasta_pac = ref_fasta_pac,
            ref_fasta_sa = ref_fasta_sa,
            aligned_name = "mother.bam"
    }

    call BwaMem as BwaMemFather {
        input:
            end1 = father_fastqs[0],
            end2 = father_fastqs[1],
            ref_fasta = ref_fasta,
            ref_dict = ref_dict,
            ref_fasta_fai =  ref_fasta_fai,
            ref_fasta_amb = ref_fasta_amb,
            ref_fasta_ann = ref_fasta_ann,
            ref_fasta_bwt = ref_fasta_bwt,
            ref_fasta_pac = ref_fasta_pac,
            ref_fasta_sa = ref_fasta_sa,
            aligned_name = "father.bam"
    }

    call ASM.SelectReadsFromRegion as SelectMHCFromMother {
        input:
            bam = BwaMemMother.aligned_bam,
            bai = BwaMemMother.aligned_bai,
            region = mhc_region
    }

    call ASM.SelectReadsFromRegion as SelectMHCFromFather {
        input:
            bam = BwaMemFather.aligned_bam,
            bai = BwaMemFather.aligned_bai,
            region = mhc_region
    }

    # select corrected reads from MHC
    call ASM.SelectReadsFromRegion as SelectMHCFromCorrected {
        input:
            bam = child_corrected_bam,
            bai = child_corrected_bai,
            region = mhc_region
    }

    # select remaining reads from MHC
    call ASM.SelectReadsFromRegion as SelectMHCFromRemaining {
        input:
            bam = child_remaining_bam,
            bai = child_remaining_bai,
            region = mhc_region
    }

#    call AssembleReadsWithTrioCanu {
#        input:
#            mother = CombineMother.reads,
#            father = CombineFather.reads,
#            child  = SelectMHCFromCorrected.reads,
#
#            target_size = "5m",
#            platform = "PACBIO",
#            is_corrected = true,
#            prefix = prefix
#    }
#
#    # correct/trim remaining reads
#    call ASM.CorrectAndTrimReadsWithCanu as CorrectMHCFromRemaining {
#        input:
#            reads = SelectMHCFromRemaining.reads,
#            target_size = "5m",
#            platform = platform,
#            is_corrected = false,
#            prefix = "mhc"
#    }
#
#    # combine corrected and remaining reads
#    call ASM.CombineReads as CombineReads {
#        input:
#            reads = [ SelectMHCFromCorrected.reads, CorrectMHCFromRemaining.trimmed_reads ]
#    }
#
#    # assemble combined reads
#    call ASM.AssembleReadsWithCanu as AssembleMHCFromCombined {
#        input:
#            reads = CombineReads.reads,
#            target_size = "5m",
#            platform = platform,
#            prefix = prefix
#    }
#
#    # align assembly of combined reads
#    call ASM.AlignContigs as AlignedMHCFromCombined {
#        input:
#            contigs = AssembleMHCFromCombined.contigs_fasta,
#            ref_fasta = ref_fasta,
#            SM = "test",
#            ID = "combined",
#            PL = "PACBIO",
#            is_corrected = true,
#            prefix = prefix
#    }
#
#    # call haploid variants on combined assembly
#    call ASM.CallHaploidVariants as CallHaploidVariantsFromCombined {
#        input:
#            bam = AlignedMHCFromCombined.aligned_bam,
#            bai = AlignedMHCFromCombined.aligned_bai,
#            ref_fasta = ref_fasta,
#            prefix = prefix
#    }
#
#    output {
#        File report          = AssembleMHCFromCombined.report
#
#        File contigs_fasta   = AssembleMHCFromCombined.contigs_fasta
#        File unassembled     = AssembleMHCFromCombined.unassembled
#        File unitigs_fasta   = AssembleMHCFromCombined.unitigs_fasta
#
#        File contigs_layout  = AssembleMHCFromCombined.contigs_layout
#        File unitigs_layout  = AssembleMHCFromCombined.unitigs_layout
#        File unitigs_bed     = AssembleMHCFromCombined.unitigs_bed
#
#        File contigs_gfa     = AssembleMHCFromCombined.contigs_gfa
#        File unitigs_gfa     = AssembleMHCFromCombined.unitigs_gfa
#
#        File aligned_bam     = AlignedMHCFromCombined.aligned_bam
#        File aligned_bai     = AlignedMHCFromCombined.aligned_bai
#        File calls           = CallHaploidVariantsFromCombined.calls
#    }
}

task BwaMem {
    input {
        File end1
        File end2
        File ref_fasta
        File ref_dict
        File ref_fasta_fai
        File ref_fasta_amb
        File ref_fasta_ann
        File ref_fasta_bwt
        File ref_fasta_pac
        File ref_fasta_sa
        String aligned_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(end1, "GB") + size(end2, "GB") + size(ref_fasta, "GB") + size(ref_dict, "GB") + size(ref_fasta_fai, "GB") + size(ref_fasta_amb, "GB") + size(ref_fasta_ann, "GB") + size(ref_fasta_bwt, "GB") + size(ref_fasta_pac, "GB") + size(ref_fasta_sa, "GB"))
    Int cpus = 16

    command <<<
        set -euxo pipefail

        bwa mem -M -t ~{cpus} ~{ref_fasta} ~{end1} ~{end2} | samtools sort -@~{cpus} -o ~{aligned_name} -
        samtools index ~{aligned_name}
    >>>

    output {
        File aligned_bam = "~{aligned_name}"
        File aligned_bai = "~{aligned_name}.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             30,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-asm:0.01.11"
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

task AssembleReadsWithTrioCanu {
    input {
        File mother
        File father
        File child

        String target_size
        String platform
        Boolean? is_corrected
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Boolean correct = select_first([is_corrected, false])
    String data_type = "-" + (if platform == "PACBIO" then "pacbio" else "nanopore") + "-" + (if correct then "corrected" else "raw")

    Int disk_size = 4*ceil(size(mother, "GB") + size(father, "GB") + size(child, "GB"))
    Int cpus = 8

    command <<<
        set -euxo pipefail

        # assemble haplotypes
        canu \
            -p ~{prefix} -d ./ \
            genomeSize=~{target_size} \
            -haplotypeA ~{mother} \
            -haplotypeB ~{father} \
            ~{data_type} ~{child}
    >>>

    output {
        File report          = "~{prefix}.report"

        File contigs_fasta   = "~{prefix}.contigs.fasta"
        File unassembled     = "~{prefix}.unassembled.fasta"
        File unitigs_fasta   = "~{prefix}.unitigs.fasta"

        File contigs_layout  = "~{prefix}.contigs.layout"
        File unitigs_layout  = "~{prefix}.unitigs.layout"
        File unitigs_bed     = "~{prefix}.unitigs.bed"

        File contigs_gfa     = "~{prefix}.contigs.gfa"
        File unitigs_gfa     = "~{prefix}.unitigs.gfa"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             12,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-asm:0.01.11"
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
