version 1.0

##########################################################################################
# This pipeline calls small variants (currently only SNVs) on an input ONT BAM using
# known algorithms that are specifically designed to work with long read data.
##########################################################################################

import "Structs.wdl"
import "Utils.wdl" as Utils
import "Longshot.wdl" as Longshot
import "DeepVariant.wdl" as DV
import "Clair.wdl" as Clair

workflow CallSmallVariantsONT {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai
        File ref_dict
    }

    parameter_meta {
        bam: "input BAM from which to call SNVs"
        bai: "index accompanying the BAM"

        ref_fasta:     "reference to which the BAM was aligned to"
        ref_fasta_fai: "index accompanying the reference"
        ref_dict:      "dictionary accompanying the reference"
    }

    call Utils.MakeChrIntervalList { input: ref_dict = ref_dict }

    scatter (chr_info in MakeChrIntervalList.chrs) {
        call Longshot.Longshot {
            input:
                bam           = bam,
                bai           = bai,
                ref_fasta     = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                chr           = chr_info[0]
        }

#        call DV.PEPPER {
#            input:
#                bam           = bam,
#                bai           = bai,
#                ref_fasta     = ref_fasta,
#                ref_fai       = ref_fasta_fai,
#                chr           = chr_info[0]
#        }

        call Clair.Clair {
            input:
                bam           = bam,
                bai           = bai,
                ref_fasta     = ref_fasta,
                ref_fai       = ref_fasta_fai,
                model_class   = "ont",
                chr           = chr_info[0]

        }
    }

    call MergeSNVCalls as MergeLongshotVCFs {
        input:
            vcfs = Longshot.vcf,
            ref_dict = ref_dict,
            prefix = basename(bam, ".bam") + ".longshot"
    }

#    call MergeSNVCalls as MergeDeepVariantVCFs {
#        input:
#            vcfs = PEPPER.vcf,
#            ref_dict = ref_dict,
#            prefix = basename(bam, ".bam") + ".deepvariant"
#    }
#
#    call MergeSNVCalls as MergeDeepVariantGVCFs {
#        input:
#            vcfs = PEPPER.gvcf,
#            ref_dict = ref_dict,
#            prefix = basename(bam, ".bam") + ".deepvariant.g"
#    }

    call MergeSNVCalls as MergeClairVCFs {
        input:
            vcfs = Clair.vcf,
            ref_dict = ref_dict,
            prefix = basename(bam, ".bam") + ".clair"
    }

    output {
        File longshot_vcf = MergeLongshotVCFs.vcf
        File longshot_tbi = MergeLongshotVCFs.tbi

#        File deepvariant_vcf = MergeDeepVariantVCFs.vcf
#        File deepvariant_tbi = MergeDeepVariantVCFs.tbi
#
#        File deepvariant_gvcf = MergeDeepVariantGVCFs.vcf
#        File deepvariant_gtbi = MergeDeepVariantGVCFs.tbi

        File clair_vcf = MergeClairVCFs.vcf
        File clair_tbi = MergeClairVCFs.tbi
    }
}

task MergeSNVCalls {
    input {
        Array[File] vcfs
        File ref_dict
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(vcfs, "GB")) + 1

    command <<<
        set -x

        VCF_WITH_HEADER=~{vcfs[0]}
        GREPCMD="grep"
        if [[ ~{vcfs[0]} =~ \.gz$ ]]; then
            GREPCMD="zgrep"
        fi

        $GREPCMD '^#' $VCF_WITH_HEADER | grep -v -e '^##contig' -e CHROM > header
        grep '^@SQ' ~{ref_dict} | awk '{ print "##contig=<ID=" $2 ",length=" $3 ">" }' | sed 's/[SL]N://g' >> header
        $GREPCMD -m1 CHROM $VCF_WITH_HEADER >> header

        ((cat header) && ($GREPCMD -h -v '^#' ~{sep=' ' vcfs})) | bcftools sort | bgzip > ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File vcf = "~{prefix}.vcf.gz"
        File tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             24,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-longshot:0.1.2"
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
