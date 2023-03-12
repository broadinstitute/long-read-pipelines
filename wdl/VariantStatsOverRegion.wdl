version 1.0

import "tasks/Structs.wdl"

workflow VariantStatsOverRegion {
    input {
        File snf_vcf
        File snf_tbi

        File pav_vcf
        File pav_tbi

        File phased_snp_vcf
        File phased_snp_tbi

        File regions_file
    }

    call SubsetVCF as GetSnfSVs  { input: vcf_gz = snf_vcf,        vcf_tbi = snf_tbi,        regions_file = regions_file }
    call SubsetVCF as GetPavSVs  { input: vcf_gz = pav_vcf,        vcf_tbi = pav_tbi,        regions_file = regions_file }
    call SubsetVCF as GetSNPs    { input: vcf_gz = phased_snp_vcf, vcf_tbi = phased_snp_tbi, regions_file = regions_file }

    call WhatshapStats { input: vcf = GetSNPs.subset_vcf, tbi = GetSNPs.subset_tbi }
    call BCFtoolsStats { input: vcf = GetSNPs.subset_vcf, tbi = GetSNPs.subset_tbi }
    call StatsOverPhasedSNPsInARegion {
        input:
            vcf = GetSNPs.subset_vcf,
            tbi = GetSNPs.subset_tbi,
            bcftools_stats = BCFtoolsStats.stats,
            whatshap_stats = WhatshapStats.stats
    }

    call StatsOverSVsInARegion as PavStats { input: vcf = GetPavSVs.subset_vcf, tbi = GetPavSVs.subset_tbi }
    call StatsOverSVsInARegion as SnfStats { input: vcf = GetSnfSVs.subset_vcf, tbi = GetSnfSVs.subset_tbi }

    Float pav_sv_cnt = PavStats.sv_cnt
    Float snf_sv_cnt = SnfStats.sv_cnt
    Float pav_del_cnt = PavStats.del_cnt
    Float snf_del_cnt = SnfStats.del_cnt

    Float pav_bnd_frac = if (PavStats.sv_cnt==0) then 0 else PavStats.bnd_cnt/pav_sv_cnt
    Float pav_1kp_frac = if (PavStats.sv_cnt==0) then 0 else PavStats.over_1kbp_cnt/pav_sv_cnt
    Float pav_ins_del_ratio = if (PavStats.del_cnt==0) then 0 else PavStats.ins_cnt/pav_del_cnt
    Float snf_bnd_frac = if (SnfStats.sv_cnt==0) then 0 else SnfStats.bnd_cnt/snf_sv_cnt
    Float snf_1kp_frac = if (SnfStats.sv_cnt==0) then 0 else SnfStats.over_1kbp_cnt/snf_sv_cnt
    Float snf_ins_del_ratio = if (SnfStats.del_cnt==0) then 0 else SnfStats.ins_cnt/snf_del_cnt

    output {
        File subset_snp_vcf = GetSNPs.subset_vcf
        File subset_snp_tbi = GetSNPs.subset_tbi

        Float roi_small_var_median_QUAL = StatsOverPhasedSNPsInARegion.median_QUAL
        Map[String, Int] roi_num_small_vars = StatsOverPhasedSNPsInARegion.num_small_vars
        Map[String, Float] roi_phasing_stats  = StatsOverPhasedSNPsInARegion.phasing_stats

        File subset_snf_vcf = GetSnfSVs.subset_vcf
        File subset_snf_tbi = GetSnfSVs.subset_tbi

        File subset_pav_vcf = GetPavSVs.subset_vcf
        File subset_pav_tbi = GetPavSVs.subset_tbi

        Map[String, Float] roi_sv_stats = {"PAV_tot_cnt": PavStats.sv_cnt,
                                           "SNF_tot_cnt": SnfStats.sv_cnt,
                                           "PAV_BND_frac": pav_bnd_frac,
                                           "SNF_BND_frac": snf_bnd_frac,
                                           "PAV_1k+_frac": pav_1kp_frac,
                                           "SNF_1k+_frac": snf_1kp_frac,
                                           "PAV_insdel_ratio": pav_ins_del_ratio,
                                           "SNF_insdel_ratio": snf_ins_del_ratio}
    }
}

task BCFtoolsStats {
    input {
        File vcf
        File tbi

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([vcf, tbi], "GB")) + 10

    command <<<
        set -eux

        bcftools stats ~{vcf} > "bcftools.stats.out.txt"
    >>>
    output {
        File stats = "bcftools.stats.out.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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

task WhatshapStats {
    input {
        File vcf
        File tbi

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([vcf, tbi], "GB")) + 10

    command <<<
        set -eux

        whatshap stats ~{vcf} > "whatshap.stats.out.tsv"
    >>>
    output {
        File stats = "whatshap.stats.out.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/biocontainers/whatshap:1.7--py310h30d9df9_0"
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

task StatsOverPhasedSNPsInARegion {
    input {
        File vcf
        File tbi

        File bcftools_stats
        File whatshap_stats

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([vcf, tbi], "GB")) + 10

    command <<<
        set -eux

        # custom parsing of VCF (GQ, QUAL)
        zgrep -v "^#" ~{vcf} \
            | awk -F '\t' 'BEGIN{OFS="\t"} {if($7=="PASS") print}' \
            > body.txt
        awk -F '\t' '{print $6}' body.txt \
            | sort -n \
            | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }'  \
            > median_QUAL.txt
        cat median_QUAL.txt

        # extract some facts from bcftools stats
        num_snps=$(grep -F 'number of SNPs:' ~{bcftools_stats} | awk -F ':' '{print $2}' | grep -Eo "[0-9]+")
        num_indels=$(grep -F 'number of indels:' ~{bcftools_stats} | awk -F ':' '{print $2}' | grep -Eo "[0-9]+")
        num_multiallelic=$(grep -F 'number of multiallelic sites:' ~{bcftools_stats} | awk -F ':' '{print $2}' | grep -Eo "[0-9]+")
        echo -e "snps\t${num_snps}\nindels\t${num_indels}\nmultiallelic\t${num_multiallelic}" > parsed.bcfstats.tsv

        # extract some facts from whatshap stats
        tail -n21 ~{whatshap_stats} > all.regions.whatshap.txt
        cat all.regions.whatshap.txt
        num_het=$(grep -F 'Heterozygous:' all.regions.whatshap.txt | awk '{print $2}' | grep -Eo "[0-9]+")
        num_unphased=$(grep -F 'Unphased:' all.regions.whatshap.txt | grep -Eo "[0-9]+")
        median_phase_block_sz=$(grep -F 'Median block size:' all.regions.whatshap.txt | grep -Eo "[0-9]+\.[0-9]+")
        block_ng50=$(grep -F 'Block NG50:' all.regions.whatshap.txt | awk -F ':' '{print $2}' | grep -Eo "[0-9]+")
        echo -e "num_het\t${num_het}\nnum_unphased\t${num_unphased}\nmedian_phase_block_sz\t${median_phase_block_sz}\nblock_ng50\t${block_ng50}" > parsed.whatshapstats.tsv
    >>>
    output {
        Float median_QUAL = read_float("median_QUAL.txt")
        Map[String, Int] num_small_vars = read_map("parsed.bcfstats.tsv")
        Map[String, Float] phasing_stats  = read_map("parsed.whatshapstats.tsv")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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

task StatsOverSVsInARegion {
    input {
        File vcf
        File tbi

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([vcf, tbi], "GB")) + 10

    command <<<
        set -eux

        # filter down to passing SVs (50+)
        zgrep -v "^#" ~{vcf} \
            | awk -F '\t' 'BEGIN{OFS="\t"} {if($7=="PASS" || $7==".") print}' \
            | awk -F '\t' 'BEGIN{OFS="\t"} match($8, /SVLEN=(-)?[0-9]+/){s=substr($8, RSTART+6, RLENGTH-6); l=s+0; if (l>49 || l<-49) print}' \
            > body.txt
        # number of variants
        wc -l body.txt | awk '{print $1}' > total.svs.count
        num_sv=$(cat total.svs.count)
        if [[ ${num_sv} -eq 0 ]]; then
            echo "0" > bnds.count
            echo "0" > ins.count
            echo "0" > del.count
            echo "0" > over.1kbp.svs.count
            exit 0
        fi
        if grep -q 'SVTYPE=INS' body.txt; then
            grep -cF 'SVTYPE=INS' body.txt > ins.count
        else
            echo "0" > ins.count
        fi
        if grep -q 'SVTYPE=DEL' body.txt; then
            grep -cF 'SVTYPE=DEL' body.txt > del.count
        else
            echo "0" > del.count
        fi
        # number of BNDs
        if grep -q 'SVTYPE=BND' body.txt; then
            grep -cF 'SVTYPE=BND' body.txt > bnds.count
        else
            echo "0" > bnds.count
        fi
        # number of calls over 1kbp
        grep -Eo "SVLEN=(-)?[0-9]+" body.txt \
            | awk -F '=' '{print $2}' \
            | sed 's#^-##' \
            | awk '{if($1>1000) print}' | wc -l | awk '{print $1}' \
        > over.1kbp.svs.count
    >>>
    output {
        Int sv_cnt  = read_int("total.svs.count")
        Int ins_cnt = read_int("ins.count")
        Int del_cnt = read_int("del.count")
        Int bnd_cnt = read_int("bnds.count")
        Int over_1kbp_cnt = read_int("over.1kbp.svs.count")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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








task SubsetVCF {
    input {
        File vcf_gz
        File vcf_tbi
        File regions_file
        String prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([vcf_gz, vcf_tbi], "GB")) + 1

    command <<<
        set -euxo pipefail

        bcftools view \
            -f "PASS,." \
            ~{vcf_gz} \
            -R ~{regions_file} \
            | bgzip > ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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