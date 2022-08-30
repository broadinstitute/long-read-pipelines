version 1.0

import "tasks/Structs.wdl"
import "tasks/AlignReads.wdl" as AR
import "AoU_Mitochondria_Canu_filteredReads.wdl" as AoU
import "tasks/Quast.wdl" as Quast
import "tasks/CallAssemblyVariants.wdl" as  CallAssemblyVariants
import "tasks/Canu.wdl" as Canu
import "tasks/Clair_mito.wdl" as Clair_Mito

workflow Trim_Contigs {
    input {
        File assembly_fasta
        File reads
        File bam
        String preset
#        String prefix

    }

    call Filter_Contigs {
        input:
        assembly_fasta = assembly_fasta
    }

    call Self_Align {
        input:
        filtered_contigs = Filter_Contigs.filtered_contigs
    }

    call AoU.RG_Parsing as Parsing {input: bam = bam}
    String RG = "@RG\\tID:~{Parsing.ID}\\tSM:~{Parsing.SM}\\tPL:~{Parsing.PL}\\tPU:~{Parsing.PU}"


    scatter (pair in zip(Self_Align.trimmed_contigs, Self_Align.trimmed_contigs_idx)) {
        call AR.Minimap2 as Minimap2 {
            input:
                reads = [reads],
                ref_fasta = pair.left,
                map_preset = "map-hifi",
                RG = RG
        }
        call Clair_Mito.Clair as Clair_Mito {
            input:
                bam = Minimap2.aligned_bam,
                bai = Minimap2.aligned_bai,
                ref_fasta = pair.left,
                ref_fasta_fai = pair.right,
                preset = preset

        }
    }

    scatter (vcf in Clair_Mito.pileup_vcf) {
    call Count_VCF {
        input:
            single_pileup_vcf = vcf
        }
    }




    output{
        Array[File] trimmed_candidate_contigs = Self_Align.trimmed_contigs
        Array[File] trimmed_cadidate_fai = Self_Align.trimmed_contigs_idx
        Array[File] aligned_bam = Minimap2.aligned_bam
        Array[File] aligned_bai = Minimap2.aligned_bai
        Array[File?] full_alignment_vcf = Clair_Mito.full_alignment_vcf
        Array[File?] gvcf = Clair_Mito.gvcf
#        Array[File] gvcf_tbi = Clair_Mito.gvcf_tbi
        Array[File?] merged_vcf = Clair_Mito.vcf
        Array[File] pileup_vcf = Clair_Mito.pileup_vcf
        Array[File] pileup_vcf_tbi = Clair_Mito.pileup_vcf_tbi
    }
}


task Filter_Contigs {
    input {
        File assembly_fasta
    }

    parameter_meta {
        assembly_fasta: "Hifiasm Assembly Fasta File"
    }

    command <<<
        set -euxo pipefail
        awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >16000 && length($2) < 30000 {print ">"$0}' ~{assembly_fasta} > filtered_contigs.fasta

    >>>

    output {
        File filtered_contigs = "filtered_contigs.fasta"
    }

    #########################
    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}


task Self_Align {
    input {
        File filtered_contigs
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        filtered_contigs: "Filtered contigs based on genome length"

    }



    command <<<
        python <<CODE
        import mappy as mp
        import subprocess
        with open("~{filtered_contigs}", "r") as f:
            assembly = f.readlines()
            for i in assembly:
                if i.startswith(">"):
                    contig_name = i.strip().split('>')[1]
                else:
                    seqlen = len(i)
                    split_pos = seqlen // 2
                    left_half = i[0:split_pos]
                    right_half = i[split_pos::]
                    left_header = contig_name+'_left'
                    right_header = contig_name+'_right'


                    with open(left_header+'.fa','w') as lfa:
                        lfa.write('>'+left_header)
                        lfa.write('\n')
                        lfa.write(left_half)
                        lfa.close()

                    with open(right_header+'.fa','w') as rfa:
                        rfa.write('>'+right_header)
                        rfa.write('\n')
                        rfa.write(right_half)
                        rfa.close()

                    aligner = mp.Aligner(left_header+'.fa')
                    if not aligner: raise Exception("ERROR: failed to load/build index")
                    for name, seq, qual in mp.fastx_read(right_header+'.fa'):
                        for hit in aligner.map(seq):
                            if hit:
                                soft_clipped = hit.q_st
                                trimmed_seq = left_half + right_half[0:hit.q_st]

                                with open(contig_name+'_trimmed.fa', 'w') as tfa:
                                    tfa.write('>'+contig_name+'_trimmed')
                                    tfa.write('\n')
                                    tfa.write(trimmed_seq)
                                    tfa.close()

                                subprocess.run(["samtools","faidx",contig_name+"_trimmed.fa"])
        CODE
    >>>

    output {
        Array[File] trimmed_contigs = glob("*_trimmed.fa")
        Array[File] trimmed_contigs_idx = glob("*_trimmed.fa.fai")

    }

#     Int disk_size = 2*ceil(size(filtered_contigs, "GB"))

    ########################
    #########################
#    RuntimeAttr default_attr = object {
#        cpu_cores:          1,
#        mem_gb:             8,
#        disk_gb:            disk_size,
#        boot_disk_gb:       10,
#        preemptible_tries:  1,
#        max_retries:        0,
#        docker:             "us.gcr.io/broad-dsp-lrma/lr-c3poa:2.2.2"
#    }
#    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
#    runtime {
#        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
#        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
#        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
#        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
#        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
#        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
#        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
#    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-c3poa:2.2.2"
    }

}


task Count_VCF {

    input{
        File single_pileup_vcf

    }

    command <<<
        set -euxo pipefail
        zcat < ~{single_pileup_vcf} > pileup_unzip.vcf
        n_snp = $(grep -v "^#" pileup_unzip.vcf | awk -F "\t" '{a=length($4); if (a==1) print $4}' | grep -c '[A-Za-z]')

    >>>



    output {
        # tsv file? fasta file?


    }

    #########################
    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }

}



#task Align_to_Candidates {
#
#    input{
#        Array[File] trimmed_contigs
#        File assembly_fasta
#    }
#
#    parameter_meta {
#        trimmed_contigs: "trimmed contigs"
#    }
#
#    command <<<
#        minimap2  -aYL --MD -t 8 ~{ref} ~{con}
#    >>>
#
#    output {
#        Array[File] candidates_alignment_results = glob("*.bam")
#    }
#
#    runtime {
#        disks: ""
#        docker: ""
#    }
#}


