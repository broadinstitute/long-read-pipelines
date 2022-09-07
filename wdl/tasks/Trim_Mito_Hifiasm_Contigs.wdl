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

        call Count_Variants {
            input:
                vcfgz = Clair_Mito.vcf
            }

    }

    call Find_Min {
        input:
            variant_count = Count_Variants.count
    }

    scatter (txt in Find_Min.min_idx) {
        Int idx = read_int(txt)
        File custom_reference = Self_Align.trimmed_contigs[idx]
#        String test = idx
    }


    output{
        Array[File] trimmed_candidate_contigs = Self_Align.trimmed_contigs
        Array[File] trimmed_cadidate_fai = Self_Align.trimmed_contigs_idx
        Array[File] aligned_bam = Minimap2.aligned_bam
        Array[File] aligned_bai = Minimap2.aligned_bai
        Array[File?] full_alignment_vcf = Clair_Mito.full_alignment_vcf
        Array[File?] merged_vcf = Clair_Mito.vcf
        Array[Int] variant_counts = Count_Variants.count
        Array[File] min_idx = Find_Min.min_idx
        Array[File] custom_reference_all = custom_reference
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


    ###################
    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-c3poa:2.2.2"
    }

}

task Count_Variants {

    input{
        File? vcfgz
        RuntimeAttr? runtime_attr_override
    }

    command <<<

        zcat < ~{vcfgz} > merge_output.vcf
        counts=$(grep -v "^#" merge_output.vcf | awk -F "\t" '{a=length($4); if (a==1) print $4}' | grep -c '[A-Za-z]')
        if [[ $counts -eq 0 ]];
        then
            echo  0 > counts.txt
        else
            echo $counts > counts.txt
        fi
    >>>


    output {
        Int count = read_int("counts.txt")
    }

    #########################

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }

}

task Find_Min {
    input {
        Array[Int] variant_count
    }

    Int n = length(variant_count)-1

    command <<<
        set -eux
        seq 0 ~{n} > indices.txt
        echo "~{sep='\n' variant_count}" > v_counts.txt
        paste -d' ' indices.txt v_counts.txt > counts_index.txt
        sort -k2 -n counts_index.txt > sorted_counts.txt


        min_val=$(head -n 1 sorted.txt | awk '{print $2}')

        for idx in $(awk -F " " '{if ($2=="$min_val") print $1}' sorted_counts.txt); do
            echo $idx > ${idx%.*}.txt
        done
        >>>



    output {
        Array[File] min_idx = glob("^[0-9]*[1-9][0-9]*.txt")
    }

    ###################
    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}











