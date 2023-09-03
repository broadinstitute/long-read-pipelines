version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/VariantUtils.wdl" as VU
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow SummarizeSmallVariants {
    input {
        File vcf
        File tbi
        File ref_map_file

        String prefix
        String gcs_out_root_dir

        Int threads = 2
        Int memory = 8
        Boolean is_gvcf = false
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/SummarizeSmallVariants/~{prefix}"

    String outbase = basename(basename(vcf, ".g.vcf.gz"), ".vcf.gz")
    String ext = if (is_gvcf) then "g.vcf.gz" else "vcf.gz"

    call VU.FilterSmallVariantsOnQual as QualFilter {
        input:
            vcf = vcf, tbi = tbi,
            ref_fasta = ref_map['fasta'], ref_fai = ref_map['fai'], ref_dict = ref_map['dict'],
            qual_filter = 40,
            outbase = outbase, extension = ext
    }

    call VU.FilterToPASSVCF as HardFilter {
        input:
            raw_vcf = QualFilter.ft_vcf, outbase = outbase, extension = ext
    }

    if (! is_gvcf) { # gVCF SNP count doesn't say much
        call VU.CountSNPs {
            input:
                vcf = HardFilter.pass_vcf, tbi = HardFilter.pass_vcf_tbi,
                ref_fasta = ref_map['fasta'], ref_fai = ref_map['fai'], ref_dict = ref_map['dict']
        }
    }

    call VU.ComputeBasicSmallVariantStats {input: vcf = HardFilter.pass_vcf, threads = threads, memory = memory }

    # Finalize data
    call FF.FinalizeToFile as FinalizeBasicSmallVariantStats { input: outdir = outdir, file = ComputeBasicSmallVariantStats.vcf_stats }

    call GetFolder { input: ff = vcf }
    call FF.FinalizeToFile as FinalizeQualFilteredVCF    { input: outdir = GetFolder.folder, file = QualFilter.ft_vcf}
    call FF.FinalizeToFile as FinalizeQualFilteredVCFTbi { input: outdir = GetFolder.folder, file = QualFilter.ft_tbi}

    output {
        File qual_filtered_vcf = FinalizeQualFilteredVCF.gcs_path
        File qual_filtered_vcf_tbi = FinalizeQualFilteredVCFTbi.gcs_path
        Int? ft_pass_snp_cnt = CountSNPs.count

        File small_variant_stats = FinalizeBasicSmallVariantStats.gcs_path
    }
}

task GetFolder {
    input {
        String ff
    }

    command <<<
        set -eux
        echo ~{ff} | awk -F '/' 'BEGIN{OFS="/"} NF{NF-=1};1' > folder.txt
    >>>

    output {
        String folder = read_string("folder.txt")
    }

    runtime {
        cpu: 1
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}
