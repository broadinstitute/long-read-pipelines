version 1.0

import "tasks/Structs.wdl"
import "tasks/RepeatMasker.wdl" as RM

workflow AnnotateVCFIns {
    parameter_meta {
        VCF: "VCF containing INS calls with sequences in the ALT field"
    }

    input {
        File VCF
    }

    call VCF_INS_to_fa { input: VCF = VCF }

    call RM.RepeatMasker { input: fasta = VCF_INS_to_fa.INS_fa }

    output {
        File RMout = RepeatMasker.RMout
    }
}

task VCF_INS_to_fa {
    input {
        File VCF

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(VCF, "GB"))+20

    String prefix = basename(VCF, ".vcf.gz")

    command <<<
        set -euxo pipefail

        bcftools view -i 'INFO/SVLEN>50 && INFO/SVTYPE=="INS"' ~{VCF} | bcftools view -e 'ALT ~ "<"' > ~{prefix}.INS.vcf
        bcftools query -f '>%CHROM:%POS\n%ALT\n' ~{prefix}.INS.vcf > ~{prefix}_INS.fa
    >>>

    output {
        File INS_fa = '~{prefix}_INS.fa'
    }

     #########################
     RuntimeAttr default_attr = object {
         cpu_cores:          1,
         mem_gb:             2,
         disk_gb:            disk_size,
         boot_disk_gb:       10,
         preemptible_tries:  3,
         max_retries:        0,
         docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
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
