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
        File INSfa = VCF_INS_to_fa.INS_fa
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

        bcftools view -i 'INFO/SVLEN>=50 && INFO/SVTYPE=="INS"' ~{VCF} | bcftools view -e 'ALT ~ "<"' > ~{prefix}.INS.vcf
        # the next line filters out "INS" calls in the integrated callset that are actually INV
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ~{prefix}.INS.vcf |awk 'length($3)==1 {print ">"$1":"$2";"$3"\n"$4}'> ~{prefix}_INS.tmp.fa
        
        #rename duplicate fasta IDs
        seqkit rename -N1 ~{prefix}_INS.tmp.fa > ~{prefix}_INS.fa
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
         docker:             "quay.io/ymostovoy/lr-utils-basic:latest"
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
