version 1.0


##########################################################################################
## A workflow
##########################################################################################
import "tasks/Utils.wdl" as Utils
import "tasks/Hifiasm.wdl" as HA
import "tasks/AlignReads.wdl" as AR
import "tasks/Quast.wdl" as Quast
import "tasks/CallAssemblyVariants.wdl" as  CallAssemblyVariants

task RG_Parsing {
    input{
        File bam
    }

    parameter_meta {
        bam: "GCS path to raw subread bam"
    }

    command <<<
        set -euxo pipefail
        samtools view -H ~{bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^ID:' | sed 's/ID://g' > ID.txt
        samtools view -H ~{bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g' > SM.txt
        samtools view -H ~{bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^PL:' | sed 's/PL://g' > PL.txt
        samtools view -H ~{bam} | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^PU:' | sed 's/PU://g' > PU.txt
    >>>

    output{
        String ID = read_string("ID.txt")
        String SM = read_string("SM.txt")
        String PL = read_string("PL.txt")
        String PU = read_string("PU.txt")
    }

    runtime {
        cpu:                    1
        memory:                 3 + " GiB"
        disks:                  "local-disk " +  10 + " HDD"
        bootDiskSizeGb:         10
        preemptible:            3
        maxRetries:             2
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}


workflow MitochondriaProcessing{
    input{
        File bam
        File bai
        String locus
        String prefix
        File ref_fasta
        File ref_fai
        String participant_name

    }

    parameter_meta{
        bam:        "GCS path to raw subread bam"
        bai:        "index for bam file"
        locus:      "genomic locus to select"
        prefix:     "prefix for output bam and bai file names"
        ref_fasta:  "chrM reference fasta"
        map_preset: "preset to be used for minimap2 parameter '-x'"
        ref_fai:    "index of fa"
    }


    #String RG = "@RG\\tID:~{ID}\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}"


    call Utils.SubsetBam as SubsetBam {input: bam = bam, bai = bai, locus=locus} #task_var = workflow_var
    call RG_Parsing as Parsing {input: bam = SubsetBam.subset_bam}
    call Utils.BamToFastq as BamToFastq {input: bam = SubsetBam.subset_bam, prefix = prefix}

    String ID = Parsing.ID
    String SM = Parsing.SM
    String PL = Parsing.PL
    String PU = Parsing.PU
    String RG = "@RG\\tID:~{ID}\\tSM:~{SM}\\tPL:~{PL}\\tPU:~{PU}"

    call HA.Hifiasm as Hifiasm {input: reads = BamToFastq.reads_fq, prefix = prefix}
    call AR.Minimap2 as Minimap2 {input: reads = [Hifiasm.fa], ref_fasta = ref_fasta, map_preset = "map-hifi", RG = RG}
    call Quast.Quast as Quast {input: ref = ref_fasta, assemblies = [Hifiasm.fa]}


    call CallAssemblyVariants.CallAssemblyVariants as  CallAssemblyVariants {input: asm_fasta = Hifiasm.fa,
                                                                ref_fasta = ref_fasta,
                                                                participant_name = participant_name,
                                                                prefix = prefix}

    output{ File chrM_bam = SubsetBam.subset_bam
            File chrM_bam_bai = SubsetBam.subset_bai
            File reads_fq = BamToFastq.reads_fq
            File gfa = Hifiasm.gfa
            File fa = Hifiasm.fa
            File chrM_aligned_bam = Minimap2.aligned_bam
            File chrM_aligned_bai = Minimap2.aligned_bai

            File report_html = Quast.report_html
            File report_txt = Quast.report_txt
            File report_pdf = Quast.report_pdf
            Array[File] quast_plots = Quast.plots

            File paf = CallAssemblyVariants.paf
            File paftools_vcf = CallAssemblyVariants.paftools_vcf}
}
















