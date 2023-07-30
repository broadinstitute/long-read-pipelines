version 1.0

import "../../structs/Structs.wdl"

task GetReadGroupInfo {
    meta {
        desciption:
        "Get some read group information Given a single-readgroup BAM. Will fail if the information isn't present."
    }

    parameter_meta {
        uBAM: "The input BAM file."
        keys: "A list of requested fields in the RG line, e.g. ID, SM, LB."
    }

    input {
        String uBAM  # not using file as call-caching brings not much benefit

        Array[String] keys
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{uBAM} | grep "^@RG" | tr '\t' '\n' > rh_header.txt

        for attribute in ~{sep=' ' keys}; do
            value=$(grep "^${attribute}" rh_header.txt | awk -F ':' '{print $2}')
            echo -e "${attribute}\t${value}" >> "result.txt"
        done
    >>>

    output {
        Map[String, String] read_group_info = read_map("result.txt")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task GetReadGroupLines {
    meta {
        desciption: "Get the @RG lines in a BAM's header. Will error if there's no read group defined in the header."
    }

    input {
        String bam
    }

    output {
        Array[String] read_group_ids = read_lines("rgids.txt")
        Array[String] read_group_lines = read_lines("read_groups.txt")
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        samtools view -H ~{bam} | grep "^@RG" > read_groups.txt

        awk -F '\t' '{print $2}' read_groups.txt | awk -F ':' '{print $2}' > rgids.txt
    >>>

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 10 HDD"
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.2"
    }
}

task GetSortOrder {
    parameter_meta {
        res: "According to the SAM spec: any one of the following values is legit [unknown, unsorted, queryname, coordinate]."
    }
    input {
        String bam
    }

    output {
        String res = read_string("result.txt")
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam}  > header.txt

        head header.txt  # cheap debug line

        grep "^@HD" header.txt | tr '\t' '\n' > hd_line.txt
        if grep -q "SO:" hd_line.txt; then
            grep "^SO:" hd_line.txt | awk -F ':' '{print $2}' > result.txt
        else
            echo "unknown" > result.txt
        fi
    >>>

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task SplitNameSortedUbam {
    meta {
        desciption: "Split a read-name-sorted unaligned BAM into chunks."
    }
    parameter_meta {
        read_cnt: "number of reads in the uBAM; providing this will reduce run time."
        n_reads: "desired number of reads per split; mutually exclusive with n_files"
        n_files: "desired number of split files; mutually exclusive with n_reads"
        uBAM: {
            localization_optional: true
        }
    }
    input {
        File uBAM
        Int? read_cnt
        Int? n_reads
        Int? n_files

        RuntimeAttr? runtime_attr_override
    }
    output {
        Array[File] split = glob("split_outputs/*.bam")
    }

    Boolean fail = defined(n_reads) == defined(n_files)  # mutually exclusive

    Int X = select_first([n_reads, n_files])
    String split_arg = if defined(n_reads) then "--SPLIT_TO_N_READS ~{X}" else "--SPLIT_TO_N_FILES ~{X}"
    String helper_arg = if (defined(read_cnt)) then "--TOTAL_READS_IN_INPUT ~{read_cnt}" else " "

    Int disk_size = 20 + ceil(11 * size(uBAM, "GiB"))

    String base = basename(uBAM, ".bam")
    String local_bam = "/cromwell_root/~{base}.bam"

    command <<<
        set -eux

        if ~{fail}; then echo "one and only one of [n_reads, n_files] must be specified" && exit 1; fi

        # prep
        time gcloud storage cp ~{uBAM} ~{local_bam}
        mkdir -p split_outputs

        # higher memory, also lower # of reads in memory given ~100 longer reads (1.5E4 bp vs 1.5E2 bp)
        gatk SplitSamByNumberOfReads \
            --java-options "-Xmx28G -Xms24G" \
            -use_jdk_deflater -use_jdk_inflater \
            --MAX_RECORDS_IN_RAM 5000 \
            -I ~{local_bam} \
            -O split_outputs \
            ~{split_arg} \
            ~{helper_arg}
    >>>
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          6,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.4.0.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
