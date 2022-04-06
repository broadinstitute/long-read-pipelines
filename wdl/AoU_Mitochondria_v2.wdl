version 1.0


import "tasks/Utils.wdl" as Utils
import "tasks/Hifiasm.wdl" as HA
import "tasks/AlignReads.wdl" as AR
import "tasks/Quast.wdl" as Quast
import "tasks/CallAssemblyVariants.wdl" as  CallAssemblyVariants
import "tasks/Canu.wdl" as Canu

task RG_Parsing {

    meta {
        description: "This workflow parses"
    }

    input{
        File bam
    }

    parameter_meta {
        bam: {
            description: "GCS path to raw subread bam",
            localization_optional: true
        }
    }

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

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

    meta {
    description:
    "This workflow subsets mitochondrial reads from whole genome bam file and calls variants from mitochondrial reads."
    }

    input{
        File bam
        File bai
        String locus
        String prefix
        File ref_fasta
        File ref_fai
        File reads_fq
        String participant_name
        Int genome_size
        Float correct_error_rate = 0.045
        Float trim_error_rate = 0.045
        Float assemble_error_rate = 0.045
        String technology
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

    call RG_Parsing as Parsing {input: bam = bam}
    String RG = "@RG\\tID:~{Parsing.ID}\\tSM:~{Parsing.SM}\\tPL:~{Parsing.PL}\\tPU:~{Parsing.PU}"

    call Canu.Canu {
        input:
            reads = reads_fq,
            genome_size = genome_size,
            correct_error_rate = correct_error_rate,
            trim_error_rate = trim_error_rate,
            assemble_error_rate = assemble_error_rate,
            prefix = prefix,
            technology = technology}
    call AR.Minimap2 as Minimap2 {input: reads = [Canu.assemble_fa], ref_fasta = ref_fasta, map_preset = "map-hifi", RG = RG}
    call Quast.Quast as Quast {input: ref = ref_fasta, assemblies = [Canu.assemble_fa]}
    call CallAssemblyVariants.CallAssemblyVariants as  CallAssemblyVariants {input: asm_fasta = Canu.assemble_fa,
                                                                ref_fasta = ref_fasta,
                                                                participant_name = participant_name,
                                                                prefix = prefix}
    output{
            File chrM_aligned_bam = Minimap2.aligned_bam
            File chrM_aligned_bai = Minimap2.aligned_bai
            File canu_contig = Canu.assemble_fa
            File canu_assemble_log = Canu.assemble_log
            File report_html = Quast.report_html
            File report_txt = Quast.report_txt
            File report_pdf = Quast.report_pdf
            Array[File] quast_plots = Quast.plots
            File paf = CallAssemblyVariants.paf
            File paftools_vcf = CallAssemblyVariants.paftools_vcf}
}
















