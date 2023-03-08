version 1.0


# Refines the genotypes of a population.
#
workflow LRRegenotype {
    input {
        File merged_vcf_gz
        File bam_addresses
        Int use_lrcaller
        Int lrcaller_genotyping_model
        Int use_cutesv
        File reference_fa
        File reference_fai
        Int n_nodes
        Int n_cpus
        Int bam_size_gb
    }
    parameter_meta {
        merged_vcf_gz: "The output of the merging step, whose genotypes must be refined."
        bam_addresses: "File containing a list of bucket addresses."
        lrcaller_genotyping_model: "Genotyping model: 1=AD, 2=VA, 3=JOINT, 4=PRESENCE, 5=VA_OLD. See https://github.com/DecodeGenetics/LRcaller/blob/144204b0edf95c55f31fdf0fbd4f02f4a36edfc1/algo.hpp#L1179"
        n_nodes: "Use this number of nodes to regenotype in parallel."
        n_cpus: "Lower bound on the number of CPUs per regenotype node."
        bam_size_gb: "Upper bound on the size of a single BAM."
    }
    
    call GetVcfToGenotype {
        input:
            merged_vcf_gz = merged_vcf_gz,
            reference_fa = reference_fa
    }
    call GetChunks {
        input:
            bam_addresses = bam_addresses,
            n_chunks = n_nodes
    }
    scatter(chunk_file in GetChunks.chunks) {
        call RegenotypeChunk { 
            input:
                chunk = chunk_file,
                vcf_to_genotype = GetVcfToGenotype.vcf_to_genotype,
                use_lrcaller = use_lrcaller,
                lrcaller_genotyping_model = lrcaller_genotyping_model,
                use_cutesv = use_cutesv,
                reference_fa = reference_fa,
                reference_fai = reference_fai,
                n_cpus = n_cpus,
                bam_size_gb = bam_size_gb
        }
    }
    call PasteGenotypedChunks {
        input:
            vcf_to_genotype = GetVcfToGenotype.vcf_to_genotype,
            genotypes = RegenotypeChunk.genotypes,
            format = RegenotypeChunk.format[0],
            use_lrcaller = use_lrcaller,
            use_cutesv = use_cutesv
    }
    output {
        File vcf_gz = PasteGenotypedChunks.vcf_gz
        File tbi = PasteGenotypedChunks.tbi
    }
}


# Outputs a .vcf file that equals $merged_vcf_gz$ but lacks any sample
# information.
#
task GetVcfToGenotype {
    input {
        File merged_vcf_gz
        File reference_fa
    }
    parameter_meta {
        merged_vcf_gz: "The output of the merging step, whose genotypes must be refined."
    }
    
    Int disk_size_gb = 20*ceil(size(merged_vcf_gz, "GB")) + ceil(size(reference_fa, "GB"))
    Int ram_size_gb = 2*ceil(size(reference_fa, "GB"))

    command <<<
        set -euxo pipefail
        
        VCF_FILE=~{merged_vcf_gz}
        VCF_FILE=${VCF_FILE%.gz}
        gunzip ~{merged_vcf_gz}
        grep '#' ${VCF_FILE} > vcf_header.txt
        N_ROWS=$(wc -l < vcf_header.txt)
        head -n $(( ${N_ROWS} - 1 )) vcf_header.txt > variants_only.vcf
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> variants_only.vcf
        tail -n +$((${N_ROWS} + 1)) ${VCF_FILE} | cut -f 1,2,3,4,5,6,7,8 >> variants_only.vcf
    >>>
    
    output {
        File vcf_to_genotype = "variants_only.vcf"
    }
    runtime {
        docker: "ubuntu:latest"
        memory: ram_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


# Creates an array of balanced lists of BAM addresses.
#
task GetChunks {
    input {
        File bam_addresses
        Int n_chunks
    }
    parameter_meta {
        bam_addresses: "File containing a list of bucket addresses."
    }
    
    command <<<
        set -euxo pipefail
        
        N_LINES=$(wc -l < ~{bam_addresses})
        N_LINES_PER_CHUNK=$(( (${N_LINES} + ~{n_chunks} - 1) / ~{n_chunks} ))
        split -d -l ${N_LINES_PER_CHUNK} ~{bam_addresses} chunk-
    >>>
    output {
        Array[File] chunks = glob("chunk-*")
    }
    runtime {
        docker: "ubuntu:latest"
    }
}


# Genotypes $vcf_to_genotype$ using every remote address of a BAM in list
# $chunk$. The task returns a file with all and only the genotype columns that
# will have to be added to $vcf_to_genotype$ to create the final VCF (such
# columns include the heading with sample ID).
#
# Remark: caller parameters are taken from "Comprehensive evaluation of
# structural variant genotyping methods based on long-read sequencing data".
#
# Resource analysis for a VCF with just chr21 and 22. Node with 17 physical 
# cores.
#
# TASK     | TIME     | RAM    | CORES USED
# LRcaller | ~15 mins | 1 GB   | 17
# cutesv   | ~2 mins  | 3.5 GB | 8
#
task RegenotypeChunk {
    input {
        File chunk
        File vcf_to_genotype
        Int use_lrcaller
        Int lrcaller_genotyping_model
        Int use_cutesv
        File reference_fa
        File reference_fai
        Int n_cpus
        Int bam_size_gb
    }
    parameter_meta {
        bam_size_gb: "Upper bound on the size of a single BAM."
        lrcaller_genotyping_model: "Genotyping model: 1=AD, 2=VA, 3=JOINT, 4=PRESENCE, 5=VA_OLD. See https://github.com/DecodeGenetics/LRcaller/blob/144204b0edf95c55f31fdf0fbd4f02f4a36edfc1/algo.hpp#L1179"
    }
    
    Int ram_size_gb = n_cpus*8 + 2*( ceil(size(vcf_to_genotype,"GB")) + bam_size_gb + ceil(size(reference_fa,"GB")) )
    Int disk_size_gb = bam_size_gb + 10*ceil(size(vcf_to_genotype,"GB")) + ceil(size(reference_fa,"GB"))

    command <<<
        set -euxo pipefail
        
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        CUTESV_MIN_SUPPORTING_READS="2"  # Since we have <=4 reads per haplotype
        
        touch format.txt genotypes.txt
        i="0"
        while read BAM_FILE; do
            while : ; do
                TEST=$(gsutil -m cp ${BAM_FILE} ${BAM_FILE}.bai . && echo 0 || echo 1)
                if [ ${TEST} -eq 1 ]; then
                    echo "Error downloading file <${BAM_FILE}>. Trying again..."
                    sleep ${GSUTIL_DELAY_S}
                else
                    break
                fi
            done
            if [ ~{use_lrcaller} -eq 1 ]; then
                ${TIME_COMMAND} LRcaller --number_of_threads ${N_THREADS} -fa ~{reference_fa} $(basename ${BAM_FILE}) ~{vcf_to_genotype} genotypes.vcf    
            elif [ ~{use_cutesv} -eq 1 ]; then
                rm -rf ./cutesv_tmp; mkdir ./cutesv_tmp
                ${TIME_COMMAND} cuteSV --threads ${N_THREADS} -Ivcf ~{vcf_to_genotype} --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.8 -mi 500 -md 500 --min_support ${CUTESV_MIN_SUPPORTING_READS} --genotype -L -1 $(basename ${BAM_FILE}) ~{reference_fa} genotypes.vcf ./cutesv_tmp
            fi
            N_LINES=$(grep '#' genotypes.vcf | wc -l)
            rm -f $(basename ${BAM_FILE}) $(basename ${BAM_FILE}).bai
            echo "FORMAT" > format.txt
            tail -n +$(( ${N_LINES} + 1 )) genotypes.vcf | cut -f 9 >> format.txt
            INDIVIDUAL=$(basename -s .bam ${BAM_FILE})
            echo ${INDIVIDUAL} > new_genotypes.txt
            if [ ~{use_lrcaller} -eq 1 ]; then
                GENOTYPE_COLUMN=$(( 9 + ~{lrcaller_genotyping_model}))
            else
                GENOTYPE_COLUMN="10"
            fi
            tail -n +$(( ${N_LINES} + 1 )) genotypes.vcf | cut -f ${GENOTYPE_COLUMN} >> new_genotypes.txt
            rm -f genotypes.vcf
            if [ $i -eq 0 ]; then
                mv new_genotypes.txt genotypes.txt
                i="1"
            else
                paste genotypes.txt new_genotypes.txt > tmp.txt
                mv tmp.txt genotypes.txt; rm -f tmp.txt
            fi
            head genotypes.txt
        done < ~{chunk}
    >>>

    output {
        File format = "format.txt"
        File genotypes = "genotypes.txt"
    }
    runtime {
        docker: "fcunial/lr-genotyping"
        cpu: n_cpus
        memory: ram_size_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


# Adds the new genotypes of every chunk as columns of a new VCF file.
#
task PasteGenotypedChunks {
    input {
        File vcf_to_genotype
        Array[File] genotypes
        File format
        Int use_lrcaller
        Int use_cutesv
    }
    parameter_meta {
        vcf_to_genotype: "The .vcf file that was used as input to $RegenotypeChunk$."
        format: "Any of the format files returned by $RegenotypeChunk$ (they are assumed to be all equal)."
    }

    Int disk_size_gb = 2*( ceil(size(vcf_to_genotype,"GB")) + ceil(size(genotypes,"GB")) )
    String first_genotyped_file = genotypes[0]

    command <<<
        set -euxo pipefail
        
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Building the VCF header (excluding the column title line).
        grep '#' ~{vcf_to_genotype} > header.txt
        N_ROWS=$(wc -l < header.txt)
        head -n $(( ${N_ROWS} - 1 )) header.txt > out_header.txt
        if [ ~{use_lrcaller} -eq 1 ]; then
            echo "##FORMAT=<ID=AD,Number=3,Type=Integer,Description=\"Allelic depths from alignment supporting ref and alt allele and total number of reads\">" >> out_header.txt
            echo "##FORMAT=<ID=VA,Number=3,Type=Integer,Description=\"Allelic depths from bam file supporting ref and alt allele and total number of reads\">" >> out_header.txt
            echo "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PHRED-scaled genotype likelihoods\">" >> out_header.txt
        elif [ ~{use_cutesv} -eq 1 ]; then
            echo "##FILTER=<ID=q5,Description=\"Quality below 5\">" >> out_header.txt
            echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> out_header.txt
            echo "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# High-quality reference reads\">" >> out_header.txt
            echo "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# High-quality variant reads\">" >> out_header.txt
            echo "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"# Phred-scaled genotype likelihoods rounded to the closest integer\">" >> out_header.txt
            echo "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"# Genotype quality\">" >> out_header.txt
        fi
        
        # Building the VCF body (including the column title line).
        tail -n 1 header.txt > out_body.txt
        tail -n +$(( ${N_ROWS} + 1 )) ~{vcf_to_genotype} >> out_body.txt
        
        # Adding the FORMAT column (assumed to be the same for all chunks).
        paste out_body.txt ~{format} > tmp.txt
        mv tmp.txt out_body.txt
        
        # Adding the genotype columns
        for FILE in ~{sep=' ' genotypes}; do
            paste out_body.txt ${FILE} > tmp.txt
            mv tmp.txt out_body.txt
        done
        
        # Building the output VCF
        cat out_header.txt out_body.txt > out.vcf
        bgzip --threads ${N_THREADS} out.vcf
        tabix out.vcf.gz
    >>>

    output {
        File vcf_gz = "out.vcf.gz"
        File tbi = "out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/lr-genotyping"
        cpu: 1
        memory: "8 GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
