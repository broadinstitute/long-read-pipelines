version 1.0

import "tasks/Structs.wdl"
import "tasks/Utils.wdl"

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

    call Utils.Cat as merge_calls { input: files = PALMER.calls, out = "~{prefix}_calls.txt" , has_header = true }
    call Utils.Cat as merge_TSD_reads { input: files = PALMER.TSD_reads, out = "~{prefix}_TSD_reads.txt" , has_header = true }

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
          mkdir LINE
          mkdir SVA
          mkdir ALU
          mkdir HERVK
          /PALMER/PALMER --input ~{bam} \
                  --ref_fa ~{ref_fa} \
                  --ref_ver GRCh38 \
                  --type LINE \
                  --mode ~{mode} \
                  --output LINE \
                  --chr ~{chrom} \
                  --workdir ${dir}/LINE/
          /PALMER/PALMER --input ~{bam} \
                  --ref_fa ~{ref_fa} \
                  --ref_ver GRCh38 \
                  --type SVA \
                  --mode ~{mode} \
                  --output SVA \
                  --chr ~{chrom} \
                  --workdir ${dir}/SVA/
          /PALMER/PALMER --input ~{bam} \
                  --ref_fa ~{ref_fa} \
                  --ref_ver GRCh38 \
                  --type ALU \
                  --mode ~{mode} \
                  --output ALU \
                  --chr ~{chrom} \
                  --workdir ${dir}/ALU/
          /PALMER/PALMER --input ~{bam} \
                  --ref_fa ~{ref_fa} \
                  --ref_ver GRCh38 \
                  --type HERVK \
                  --mode ~{mode} \
                  --output HERVK \
                  --chr ~{chrom} \
                  --workdir ${dir}/HERVK/

          sed "s/$/\tSVA/"   SVA/SVA_calls.txt > ~{prefix}_calls.txt
          sed "s/$/\tLINE/"  LINE/LINE_calls.txt   | grep -v "^cluster_id" >> ~{prefix}_calls.txt
          sed "s/$/\tALU/"   ALU/ALU_calls.txt     | grep -v "^cluster_id" >> ~{prefix}_calls.txt
          sed "s/$/\tHERVK/" HERVK/HERVK_calls.txt | grep -v "^cluster_id" >> ~{prefix}_calls.txt

          sed "s/$/\tSVA/"   SVA/SVA_TSD_reads.txt > ~{prefix}_TSD_reads.txt
          sed "s/$/\tLINE/"  LINE/LINE_TSD_reads.txt   | grep -v "^cluster_id" >> ~{prefix}_TSD_reads.txt
          sed "s/$/\tALU/"   ALU/ALU_TSD_reads.txt     | grep -v "^cluster_id" >> ~{prefix}_TSD_reads.txt
          sed "s/$/\tHERVK/" HERVK/HERVK_TSD_reads.txt | grep -v "^cluster_id" >> ~{prefix}_TSD_reads.txt
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
