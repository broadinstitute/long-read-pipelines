version 1.0

workflow ONTBamCheckBeforeIngestion {
    meta {
        description: "For verification of some format concerns regarding provided ONT BAMs. Project specific."
    }

    input {
        File bam
    }

    call CheckReadgroups { input: bam = bam }

    output {
        Array[String] rgs_without_id = CheckReadgroups.rgs_without_id
        Array[String] rgs_without_sm = CheckReadgroups.rgs_without_sm
        Array[String] rgs_without_pu = CheckReadgroups.rgs_without_pu
        Array[String] rgs_absent_in_body = CheckReadgroups.rgs_absent_in_body
        Array[String] rgs_undef_in_header = CheckReadgroups.rgs_undef_in_header
    }
}

task CheckReadgroups {
    meta {
        description: "Check that readgroup lines are well defined in the header. Also check that reads all have the RG tag, and no readgroups defined in the header are absent in the body. Also check that reads' RG tags are defined in the header."
    }
    parameter_meta {
        bam: {
            localization_optional: true
        }
    }
    input {
        File bam
    }

    Int disk_size = 100 + ceil(size(bam, 'GiB'))

    String base = basename(bam, '.bam')
    String local_bam = "/cromwell_root/~{base}.bam"

    output {
        Array[String] rgs_without_id = read_lines("readgroups.without.ID.txt")
        Array[String] rgs_without_sm = read_lines("readgroups.without.SM.txt")
        Array[String] rgs_without_pu = read_lines("readgroups.without.PU.txt")
        Array[String] rgs_absent_in_body = read_lines("readgroups.defined.in.header.absent.in.body.txt")
        Array[String] rgs_undef_in_header = read_lines("readgroups.defined.in.body.absent.in.header.txt")
    }

    command <<<
        set -euxo pipefail

        time gcloud storage cp ~{bam} ~{local_bam}

        samtools view -H ~{local_bam} > original.header.txt
        grep -v "^@SQ" original.header.txt

        # problems we've seen so far are all about readgroup lines
        grep "^@RG" original.header.txt > readgroups.txt
        cat readgroups.txt

        #####################################################
        # verify all readgroup lines have [ID, SM, PU] fields
        rm -f readgroups.without.ID.txt \
              readgroups.without.SM.txt \
              readgroups.without.PU.txt
        touch readgroups.without.ID.txt \
              readgroups.without.SM.txt \
              readgroups.without.PU.txt
        while IFS= read -r line; do
            if [[ ! $(echo "$line" | tr '\t' '\n' | grep -q 'ID') ]]; then
                echo "$line" >> readgroups.without.ID.txt
            fi
            if [[ ! $(echo "$line" | tr '\t' '\n' | grep -q 'SM') ]]; then
                echo "$line" >> readgroups.without.SM.txt
            fi
            if [[ ! $(echo "$line" | tr '\t' '\n' | grep -q 'PU') ]]; then
                echo "$line" >> readgroups.without.PU.txt
            fi
        done < readgroups.txt

        head readgroups.without.*.txt

        #####################################################
        # verify all reads have RG tag
        samtools view -@1 ~{local_bam} \
        | grep -v "RG:Z:" \
        | awk -F '\t' 'BEGIN{OFS="\t"} {print $1}' \
        > readnames.without.readgroup.tag.txt &

        #####################################################
        # verify all reads have RG tag
        samtools view -@1 ~{local_bam} \
        | awk -F '\t' 'BEGIN{OFS="\t"} {for (i=12;i<=NF;i++){if ($i ~/RG:Z:/) {print $i}}}' \
        | awk -F ':' '{print $3}' \
        | sort | uniq \
        > uniq.readgroups.in.reads.txt &

        #####################################################
        wait

        tr '\t' '\n' < readgroups.txt | grep "^ID:" | awk -F ':' | sort > readgroups.defined.in.header.txt

        comm -23 \
            readgroups.defined.in.header.txt
            uniq.readgroups.in.reads.txt \
        > readgroups.defined.in.header.absent.in.body.txt

        comm -13 \
            readgroups.defined.in.header.txt
            uniq.readgroups.in.reads.txt \
        > readgroups.defined.in.body.absent.in.header.txt
    >>>

    #########################
    runtime {
        cpu:            6
        memory:         "24 GiB"
        disks:          "local-disk ~{disk_size} HDD"
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.1"
    }
}
