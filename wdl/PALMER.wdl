version 1.0

import "tasks/Structs.wdl"

workflow CallPALMER {
    input {
        File bam
        File bai
        File ref_fa
        String prefix
        String mode
        String MEI_type
        Array[String] contigs
    }

    scatter (contig_name in contigs) {
      call PALMER {
          input:
              bam = bam,
              bai = bai,
              ref_fa = ref_fa,
              prefix = prefix,
              mode = mode,
              MEI_type = MEI_type,
              chrom = contig_name
      }
    }

    call Cat as merge_calls { input: files = PALMER.calls, out = "~{prefix}_calls.txt" }
    call Cat as merge_TSD_reads { input: files = PALMER.TSD_reads, out = "~{prefix}_TSD_reads.txt" }

    output {
        File PALMER_calls = merge_calls.combined
        File PALMER_TSD_reads = merge_TSD_reads.combined
    }
}

task PALMER {
    input {
        File bam
        File bai
        File ref_fa
        String prefix
        String mode
        String MEI_type
        String chrom

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam:              "aligned bam, can contain long reads or assembled contigs"
        bai:              "index accompanying the BAM"
        ref_fa:           "fasta of reference used"
        prefix:           "output file prefix"
        mode:             "either raw or asm"
        MEI_type:         "type of MEI to call, options are LINE, ALU, SVA, HERVK, all"
        chrom:            "chromosome to analyze"
    }

    Int disk_size = ceil(size(bam, "GB") + size(ref_fa, "GB")) * 3

    command <<<
        set -x

        dir=$(pwd)

        if [ "~{MEI_type}" != "all" ]; then 
          /PALMER/PALMER --input ~{bam} \
                 --ref_fa ~{ref_fa} \
                 --ref_ver GRCh38 \
                 --type ~{MEI_type} \
                 --mode ~{mode} \
                 --output ~{prefix} \
                 --chr ~{chrom} \
                 --workdir ${dir}/
        else
          echo $'LINE\nSVA\nALU\nHERVK' > types
          while read MEI_type; do
            mkdir $MEI_type
            /PALMER/PALMER --input ~{bam} \
                    --ref_fa ~{ref_fa} \
                    --ref_ver GRCh38 \
                    --type $MEI_type \
                    --mode ~{mode} \
                   --output $MEI_type \
                    --chr ~{chrom} \
                    --workdir {dir}/${MEI_type}/

            touch ${MEI_type}/${MEI_type}_calls.txt
            touch ${MEI_type}/${MEI_type}_TSD_reads.txt

            sed "s/$/\t${MEI_type}/" ${MEI_type}/${MEI_type}_calls.txt >> calls.txt
            sed "s/$/\t${MEI_type}/" ${MEI_type}/${MEI_type}_TSD_reads.txt >> TSD_reads.txt
          done < types
          
          head -n1 calls.txt > ~{prefix}_calls.txt
          grep -v '^cluster_id' calls.txt >> ~{prefix}_calls.txt

          head -n1 TSD_reads.txt > ~{prefix}_TSD_reads.txt
          grep -v '^cluster_id' TSD_reads.txt >> ~{prefix}_TSD_reads.txt
        fi
    >>>

    output {
        File calls = "~{prefix}_calls.txt"
        File TSD_reads = "~{prefix}_TSD_reads.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "quay.io/ymostovoy/lr-palmer:2.0.0"
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


task Cat {
    input {
      Array[File] files
      String out

      RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(files, "GB"))

    command <<<
        set -euxo pipefail

        head -n1 ~{files[0]} > ~{out}

        echo ~{sep='\n' files} > filenames
        while read filename; do tail -n +2 $filename >> ~{out}; done < filenames
    >>>

    output {
        File combined = "~{out}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "ubuntu:16.04"
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
