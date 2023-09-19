version 1.0

import "tasks/Structs.wdl"

workflow CallxTea {
    input {
        File bam
        File bai
        File ref_fa
        String ref_name
        String sample_name
        String prefix
    }

    call xTea {
        input:
            bam = bam,
            bai = bai,
            ref_fa = ref_fa,
            ref_name = ref_name,
            sample_name = sample_name,
            prefix = prefix
    }

    output {
        File xTea_calls = xTea.outfile
    }
}

task xTea {
    input {
        File bam
        File bai
        File ref_fa
        String ref_name
        String sample_name
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:              "aligned bam of long reads"
        bai:              "index accompanying the BAM"
        ref_fa:           "fasta of reference used"
        ref_name:         "either hg38 or chm13"
        sample_name:      "sample name"
        prefix:           "output file prefix"
    }

    Int disk_size = ceil(size(bam, "GB") + size(ref_fa, "GB")) * 3

    command <<<
        set -x

        dir=$(pwd)
        mem=$(grep '^MemTotal' /proc/meminfo | awk '{ print int($2/1000000) }')
        cpus=$(grep -c '^processor' /proc/cpuinfo | awk '{ print $1 }')


        echo ~{sample_name} > sample_name_file
        echo ~{sample_name}$'\t'~{bam} > bam_file

        # generate running script
        python /xTea/xtea_long/gnrt_pipeline_local_long_read_v38.py -i sample_name_file \
                -b bam_file \
                -p ${dir}/ \
                -o submit_jobs.sh \
                --xtea /xTea/xtea_long/ \
                -n $cpus \
                -m $mem \
                -t 240.00 \
                -r ~{ref_fa} \
                --cns /rep_lib_annotation/consensus/LINE1.fa \
                --rep /rep_lib_annotation/ \
                -- min 4000 \
                -f 31 \
                -y 15 \
                --clean --rmsk /xTea/rep_lib_annotation/LINE/~{ref_name}/~{ref_name}_L1_larger_500_with_all_L1HS.out

        # run pipeline
        cd ~{sample_name}
        chmod +x run_xTEA_pipeline.sh
        ./run_xTEA_pipeline.sh

        # save outputs
        cat classified_results.txt.merged_HERV.txt \
          classified_results.txt.merged_SVA.txt \
          classified_results.txt.merged_ALU.txt \
          classified_results.txt.merged_LINE1.txt > ${dir}/~{prefix}_xTea_output.txt
    >>>

    output {
        File outfile = "~{prefix}_xTea_output.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "quay.io/ymostovoy/lr-xtea:0.1.9"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
