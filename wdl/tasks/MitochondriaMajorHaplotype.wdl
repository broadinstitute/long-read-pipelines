version 1.0

import "tasks/Structs.wdl"
import "tasks/AlignReads.wdl" as AR
import "tasks/CallAssemblyVariants.wdl" as  CallAssemblyVariants
import "tasks/Canu.wdl" as Canu
import "tasks/Clair_mito.wdl" as Clair_Mito
import "tasks/Finalize.wdl" as Finalize

workflow SelectContigs {
    input {
        File assembly_fasta
        File reads
        File bam
        String preset
        String outdir
    }

    call FilterContigs {
        input:
            assembly_fasta = assembly_fasta
    }

    call SelfAlign {
        input:
            filtered_contigs = FilterContigs.filtered_contigs
    }

    String RG = "@RG\\tID:NA\\tSM:HG00514.2"

    scatter (pair in zip(SelfAlign.trimmed_contigs, SelfAlign.trimmed_contigs_idx)) {
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

        call CountVariants {
            input:
                vcfgz = Clair_Mito.vcf
            }

    }

    call FindMin {
        input:
            variant_count = CountVariants.count,
            contigs = SelfAlign.trimmed_contigs
    }

    call Finalize.FinalizeToDir {
        input:
            files = FindMin.picked_tigs,
            outdir = outdir
    }


    output{
        Array[File] trimmed_candidate_contigs = SelfAlign.trimmed_contigs
        Array[File] trimmed_cadidate_fai = SelfAlign.trimmed_contigs_idx
        Array[File] aligned_bam = Minimap2.aligned_bam
        Array[File] aligned_bai = Minimap2.aligned_bai
        Array[File?] full_alignment_vcf = Clair_Mito.full_alignment_vcf
        Array[File?] merged_vcf = Clair_Mito.vcf
        Array[Int] variant_counts = CountVariants.count
        Array[File] picked_tigs = FindMin.picked_tigs
    }
}


task FilterContigs {
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


task SelfAlign {
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

task CountVariants {

    input{
        File? vcfgz
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        vcfgz: ".vcf.gz file"
    }

    command <<<
        zcat < ~{vcfgz} > merge_output.vcf
        grep -v "^#" merge_output.vcf | awk -F "\t" '{ref=length($4); alt=length($5); if (ref==1 && alt==1) print $4}' | wc -l | awk '{print $1}' > counts.txt
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

task FindMin {
    input {
        Array[Int] variant_count
        Array[String] contigs
    }

    parameter_meta {
        variant_count: "number of SNPs from contigs"
        contigs: "gs paths of contigs fasta files"
    }

    Int n = length(variant_count)-1

    command <<<
        set -eux
        echo "~{sep='\n' variant_count}" > v_counts.txt
        echo "~{sep='\n' contigs}" > fasta_path.txt

        paste -d' ' fasta_path.txt v_counts.txt > vcount_contigs.txt
        sort -k2 -n vcount_contigs.txt > sorted_counts.txt

        min_val=$(head -n 1 sorted_counts.txt | awk '{print $2}')
        awk -v mm=${min_val} '{if($2==mm) print $1}' sorted_counts.txt > picked_tigs.txt
        >>>


    output {
        Array[String] picked_tigs = read_lines("picked_tigs.txt")
    }

    ###################
    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}











