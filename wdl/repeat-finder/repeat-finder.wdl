
workflow RepeatFinderWorkflow {
    Array[String] ref_fasta_gz_urls
    String output_prefix

    Int match_score = 2
    Int mismatch_score = 3
    Int indel_score = 5
    Int minscore = 8
    Int max_repeat_unit_size = 50
    Int max_repeat_total_length = 1000000000

    scatter (ref_fasta_gz_url in ref_fasta_gz_urls) {
        call RepeatFinder {
            input:
                ref_fasta_gz_url=ref_fasta_gz_url,
                match_score=match_score,
                mismatch_score=mismatch_score,
                indel_score=indel_score,
                minscore=minscore,
                max_repeat_unit_size=max_repeat_unit_size

        }
    }

    call GatherDatToBed {
        input:
            input_dat_files=RepeatFinder.output_dat,
            max_repeat_unit_size=max_repeat_unit_size,
            max_repeat_total_length=max_repeat_total_length,
            output_prefix=output_prefix
    }


}

task RepeatFinder {
    String ref_fasta_gz_url

    Int match_score
    Int mismatch_score
    Int indel_score
    Int minscore
    Int max_repeat_unit_size

    Int pm = 80
    Int pi = 10

    Int disk_size = 10

    command {
        set -xe

        FASTA_FILENAME=$(basename ${ref_fasta_gz_url} | sed s/.gz// )

        wget -q ${ref_fasta_gz_url}
        gunzip "$FASTA_FILENAME".gz
        /trf409.linux64 $FASTA_FILENAME ${match_score} ${mismatch_score} ${indel_score} ${pm} ${pi} ${minscore} ${max_repeat_unit_size} -d -h -l 6 -ngs > output.dat
    }

    output {
        File output_dat="output.dat"
    }

    runtime {
        docker: "weisburd/repeat-finder@sha256:51c19bc6ee455f2c27e72c88dd3d063a6822d1bc9672a398505b398d13338efa"
        memory: "8G"
        disks: "local-disk ${disk_size} SSD"
    }
}

task GatherDatToBed {
    Array[File] input_dat_files
    Int max_repeat_unit_size
    Int max_repeat_total_length

    String output_prefix

    Int disk_size = 20

    command {
        set -xe

        cat ${sep=" " input_dat_files} > ${output_prefix}.dat
        python3 /convert_dat_output_to_bed.py --max-repeat-unit-size ${max_repeat_unit_size} --max-repeat-total-length ${max_repeat_total_length} ${output_prefix}.dat ${output_prefix}.bed
        bedtools sort -i ${output_prefix}.bed > ${output_prefix}.sorted.bed
        /IGVTools/igvtools index ${output_prefix}.sorted.bed
    }

    output {
        File combined_dat="${output_prefix}.dat"
        File output_bed="${output_prefix}.bed"
        File output_sorted_bed="${output_prefix}.sorted.bed"
        File output_sorted_bed_idx="${output_prefix}.sorted.bed.idx"
    }

    runtime {
        docker: "weisburd/repeat-finder@sha256:51c19bc6ee455f2c27e72c88dd3d063a6822d1bc9672a398505b398d13338efa"
        disks: "local-disk ${disk_size} SSD"
        memory: "8G"
    }
}
