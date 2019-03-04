
workflow RepeatFinderWorkflow {
    File input_vcf
    String output_prefix

    Int match_score = 2
    Int mismatch_score = 3
    Int indel_score = 5
    Int minscore = 8
    Int max_repeat_unit_size = 50


    call RepeatFinder {
        input:
            input_vcf=input_vcf,
            match_score=match_score,
            mismatch_score=mismatch_score,
            indel_score=indel_score,
            minscore=minscore,
            max_repeat_unit_size=max_repeat_unit_size,
            output_prefix=output_prefix
    }

}

task FilterVcfToStrAlleles {
    File input_vcf

    Int match_score
    Int mismatch_score
    Int indel_score
    Int minscore
    Int max_repeat_unit_size

    Int pm = 80
    Int pi = 10

    Int disk_size = 5

    Int min_allele_size = minscore/match_score    # don't added alleles shorter than (minscore/match_score) to the .fasta file because TRF will filter them out anyway

    command {
        set -xe

        PREFIX=$(basename ${input_vcf} | sed s/.vcf.gz// )

        python3 generate_fasta_from_vcf_alleles.py --min-allele-size ${min_allele_size} ${input_vcf}
        /trf409.linux64 $FASTA_FILENAME ${match_score} ${mismatch_score} ${indel_score} ${pm} ${pi} ${minscore} ${max_repeat_unit_size} -d -h -l 6 -ngs
    }

    output {
        File output_dat="output.dat"
    }

    runtime {
        docker: "weisburd/repeat-finder@sha256:51c19bc6ee455f2c27e72c88dd3d063a6822d1bc9672a398505b398d13338efa"
        memory: "4G"
        disks: "local-disk ${disk_size} SSD"
    }
}
