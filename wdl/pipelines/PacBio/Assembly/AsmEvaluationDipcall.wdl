version 1.0

workflow runDipcall {
    input {
        File truth_hap1_fasta
        File truth_hap2_fasta
        File eval_hap1_fasta
        File eval_hap2_fasta
        File HG38_reference_fasta
        File HG38_reference_fai
        String docker = "timd1/vcfdist:latest"
        Int verbosity = 1
    }
	call dipcall as truth_dipcall { input:
        assemblyFastaPat = truth_hap1_fasta,
        assemblyFastaMat = truth_hap2_fasta,
        referenceFasta = HG38_reference_fasta,
    }
    call dipcall as eval_dipcall { input:
        assemblyFastaPat = eval_hap1_fasta,
        assemblyFastaMat = eval_hap2_fasta,
        referenceFasta = HG38_reference_fasta,
    }

    call RunVcfdistTask {
        input:
            truth_vcf = truth_dipcall.outputVCF,
            eval_vcf = eval_dipcall.outputVCF,
            bed_file = truth_dipcall.outputBED,
            fasta_file = HG38_reference_fasta,
            docker = docker,
            verbosity = verbosity
    }

    output {
        File prs_tsv = RunVcfdistTask.prs_tsv
        File summary = RunVcfdistTask.summary
        File phasing_tsv = RunVcfdistTask.phasing_tsv
        File switch_flip_tsv = RunVcfdistTask.switch_flip_tsv
    }
}

task dipcall {
    input {
        File assemblyFastaPat
        File assemblyFastaMat
        File referenceFasta
        Boolean isMaleSample = false
        Boolean referenceIsHS38 = true
        Int memSizeGB = 64
        Int threadCount = 16
        Int diskSizeGB = 64
        String dockerImage = "rlorigro/dipcall_asm10:latest"
    }

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace
        PATH="/root/bin/samtools_1.9:$PATH"

        # get output base
        PREFIX=$(basename ~{assemblyFastaPat} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/[._][pm]at\(ernal\)*//')
        mkdir $PREFIX.dipcall

        # prep paternal
        PAT_FILENAME=$(basename -- "~{assemblyFastaPat}")
        if [[ $PAT_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFastaPat} .
            gunzip $PAT_FILENAME
            PAT_FILENAME="${PAT_FILENAME%.gz}"
        else
            ln -s ~{assemblyFastaPat}
        fi

        # prep maternal
        MAT_FILENAME=$(basename -- "~{assemblyFastaMat}")
        if [[ $MAT_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFastaMat} .
            gunzip $MAT_FILENAME
            MAT_FILENAME="${MAT_FILENAME%.gz}"
        else
            ln -s ~{assemblyFastaMat}
        fi

        # prep reference
        REF_FILENAME=$(basename -- "~{referenceFasta}")
        if [[ $REF_FILENAME =~ \.gz$ ]]; then
            cp ~{referenceFasta} .
            gunzip $REF_FILENAME
            REF_FILENAME="${REF_FILENAME%.gz}"
        else
            ln -s ~{referenceFasta}
        fi
        samtools faidx $REF_FILENAME

        # initialize script
        cmd=( /opt/dipcall/dipcall.kit/run-dipcall )

        # male samples need PAR region excluded
        if [[ ~{isMaleSample} == true ]]; then
            if [[ ~{referenceIsHS38} ]]; then
                cmd+=( -x /opt/dipcall/dipcall.kit/hs38.PAR.bed )
            else
                cmd+=( -x /opt/dipcall/dipcall.kit/hs37d5.PAR.bed )
            fi

        fi

        # finalize script
        cmd+=( $PREFIX.dipcall/$PREFIX )
        cmd+=( $REF_FILENAME )
        cmd+=( $PAT_FILENAME )
        cmd+=( $MAT_FILENAME )

        # generate makefile
        "${cmd[@]}" >$PREFIX.mak

        # run dipcall
        make -j 2 -f $PREFIX.mak

        # finalize
        rm $PREFIX.dipcall/*sam.gz
        tar czvf $PREFIX.dipcall.tar.gz $PREFIX.dipcall/
        cp $PREFIX.dipcall/$PREFIX.dip.bed $PREFIX.dipcall.bed
        cp $PREFIX.dipcall/$PREFIX.dip.vcf.gz $PREFIX.dipcall.vcf.gz

        # cleanup
        rm $REF_FILENAME
        rm $MAT_FILENAME
        rm $PAT_FILENAME

        # index the vcf
        /opt/bcftools/bcftools-1.9/bcftools index -t $PREFIX.dipcall.vcf.gz

	>>>
	output {
		File outputTarball = glob("*.dipcall.tar.gz")[0]
		File outputVCF = glob("*.dipcall.vcf.gz")[0]
		File outputTBI = glob("*.dipcall.vcf.gz.tbi")[0]
		File outputBED = glob("*.dipcall.bed")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task RunVcfdistTask {
    input {
        File truth_vcf
        File eval_vcf
        File bed_file
        File fasta_file
        Int verbosity
        String sampleid
        
        String docker
        Int disk_size_gb = ceil(size(truth_vcf, "GiB") + 10)
        Int mem_gb = 16
        Int cpu = 2
        Int preemptible = 1
    }

    command <<<
        vcfdist \
            ~{eval_vcf} \
            ~{truth_vcf} \
            ~{fasta_file} \
            -b ~{bed_file} \
            -v ~{verbosity}

        mv "precision-recall-summary.tsv" "~{sampleid}-precision-recall-summary.tsv"
        mv "phasing-summary.tsv" "~{sampleid}-phasing-summary.tsv"
    >>>

    runtime {
        docker: docker
        disks: "local-disk " + disk_size_gb + " HDD"
        memory: mem_gb + " GiB"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File prs_tsv = "~{sampleid}-precision-recall-summary.tsv"
        File summary = "summary.vcf"
        File precrec_tsv = "precision-recall.tsv"
        File query_tsv = "query.tsv"
        File truth_tsv = "truth.tsv"
        File phasing_tsv = "~{sampleid}-phasing-summary.tsv"
        File switch_flip_tsv = "switchflips.tsv"
    }
}