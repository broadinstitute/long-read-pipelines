version 1.0

import "Structs.wdl"

task MergePerChrCalls {
    input {
        Array[File] vcfs
        File ref_dict
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(vcfs, "GB")) + 1

    command <<<
        set -euxo pipefail

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

task MergeAndSortVCFs {
    meta {
        description: "Fast merging & sorting VCFs when the default sorting is expected to be slow"
    }

    input {
        Array[File] vcfs

        File ref_fasta_fai
        File? header_definitions_file

        String prefix

        RuntimeAttr? runtime_attr_override
    }
    parameter_meta {
        header_definitions_file: "a union of definition header lines for input VCFs (related to https://github.com/samtools/bcftools/issues/1629)"
    }

    Int sz = ceil(size(vcfs, 'GB'))
    Int disk_sz = if sz > 100 then 5 * sz else 375  # it's rare to see such large gVCFs, for now

    Boolean suspected_incomplete_definitions = defined(header_definitions_file)

    Int cores = 8

    # pending a bug fix (bcftools github issue 1576) in official bcftools release,
    # bcftools sort can be more efficient in using memory
    Int machine_memory = 48 # 96
    Int work_memory = ceil(machine_memory * 0.8)

    command <<<
        set -euxo pipefail

        echo ~{sep=' ' vcfs} | sed 's/ /\n/g' > all_raw_vcfs.txt

        echo "==========================================================="
        echo "starting concatenation" && date
        echo "==========================================================="
        bcftools \
            concat \
            --naive \
            --threads ~{cores-1} \
            -f all_raw_vcfs.txt \
            --output-type v \
            -o concatedated_raw.vcf.gz  # fast, at the expense of disk space
        for vcf in ~{sep=' ' vcfs}; do rm $vcf ; done

        # this is another bug in bcftools that's hot fixed but not in official release yet
        # (see bcftools github issue 1591)
        echo "==========================================================="
        echo "done concatenation, fixing header of naively concatenated VCF" && date
        echo "==========================================================="
        if ~{suspected_incomplete_definitions}; then
            # a bug from bcftools concat --naive https://github.com/samtools/bcftools/issues/1629
            set +e
            zgrep "^##" concatedated_raw.vcf.gz > header.txt
            grep -vF 'fileformat' header.txt \
                | grep -vF 'fileDate=' \
                | grep -vF 'source=' \
                | grep -vF 'contig' \
                | grep -vF 'ALT' \
                | grep -vF 'FILTER' \
                | grep -vF 'INFO' \
                | grep -vF 'FORMAT' \
                > tmp.others.txt
            touch tmp.other.txt
            set -e
            zgrep "^#CHROM" concatedated_raw.vcf.gz > tmp.sampleline.txt
            cat \
                ~{header_definitions_file} \
                tmp.others.txt \
                tmp.sampleline.txt \
                > fixed.header.txt
            rm -f tmp.*.txt && cat fixed.header.txt

            bcftools reheader \
                -h fixed.header.txt \
                -o tmp.wgs.vcf.gz \
                concatedated_raw.vcf.gz
            rm concatedated_raw.vcf.gz
        else
            mv concatedated_raw.vcf.gz tmp.wgs.vcf.gz
        fi
        bcftools reheader \
            --fai ~{ref_fasta_fai} \
            -o wgs_raw.vcf.gz \
            tmp.wgs.vcf.gz
        rm tmp.wgs.vcf.gz

        echo "==========================================================="
        echo "starting sort operation" && date
        echo "==========================================================="
        bcftools \
            sort \
            --temp-dir tm_sort \
            --output-type z \
            -o ~{prefix}.vcf.gz \
            wgs_raw.vcf.gz
        bcftools index --tbi --force ~{prefix}.vcf.gz
        echo "==========================================================="
        echo "done sorting" && date
        echo "==========================================================="
    >>>

    output {
        File vcf = "~{prefix}.vcf.gz"
        File tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cores,
        mem_gb:             "~{machine_memory}",
        disk_gb:            disk_sz,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CollectDefinitions {
    meta {
        description: "Collect (union) various definitions in vcf files, adddressing a bcftols bug: https://github.com/samtools/bcftools/issues/1629"
    }
    input {
        Array[File] vcfs

        RuntimeAttr? runtime_attr_override
    }

    Int sz = ceil(size(vcfs, 'GB'))
    Int disk_sz = if sz > 100 then 5 * sz else 375

    command <<<
        set -euxo pipefail

        zgrep "^##" ~{vcfs[0]} > header.txt
        grep -F '##fileformat' header.txt > tmp.0.txt
        grep -F '##fileDate=' header.txt > tmp.1.txt
        if grep -q -F '##source=' header.txt; then grep -F 'source=' header.txt > tmp.2.txt; fi
        touch tmp.2.txt
        grep -F '##contig=' header.txt > tmp.3.txt

        cat tmp*txt > easy.txt && rm tmp*txt

        touch tmp.alt.txt tmp.ft.txt tmp.info.txt tmp.format.txt
        for vcf in ~{sep=' ' vcfs}; do
            zgrep -F '##ALT=' "${vcf}" >> tmp.alt.txt
            zgrep -F '##FILTER=' "${vcf}" >> tmp.ft.txt
            zgrep -F '##INFO=' "${vcf}" >> tmp.info.txt
            zgrep -F '##FORMAT=' "${vcf}" >> tmp.format.txt
        done
        for txt in tmp*txt; do
            sort "${txt}" | uniq > "${txt}.union"
        done
        cat tmp.alt.txt.union tmp.ft.txt.union tmp.info.txt.union tmp.format.txt.union > hard.txt
        cat easy.txt hard.txt > definitions.union.txt
    >>>

    output {
        File union_definitions = "definitions.union.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_sz,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task GetVCFSampleName {
    meta {
        description: "Currently mostly used for extracting sample name in fingerprinting genotyped VCF"
    }
    input {
        File fingerprint_vcf
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        fingerprint_vcf: "Assumed to be genotyped, and hold only one sample (other samples will be ignored)."
    }

    command <<<
        set -eux

        GREPCMD="grep"
        if [[ ~{fingerprint_vcf} =~ \.gz$ ]]; then
            GREPCMD="zgrep"
        fi
        "${GREPCMD}" \
            "^#CHROM" \
            ~{fingerprint_vcf} \
            | awk '{print $10}' \
            > sample_name.txt
    >>>

    output {
        String sample_name = read_string("sample_name.txt")
    }

    ###################
    runtime {
        cpu: 2
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task SubsetVCF {
    input {
        File vcf_gz
        File vcf_tbi
        String locus
        String prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([vcf_gz, vcf_tbi], "GB")) + 1

    command <<<
        set -euxo pipefail

        bcftools view ~{vcf_gz} --regions ~{locus} | bgzip > ~{prefix}.vcf.gz
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
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

task ZipAndIndexVCF {

    meta {
        description: "gZip plain text VCF and index it."
    }

    input {
        File vcf
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(vcf, ".vcf")
    Int proposed_disk = 3*ceil(size(vcf, "GB")) + 1
    Int disk_size = if (proposed_disk > 100) then proposed_disk else 100

    command <<<
        cp ~{vcf} ~{prefix}.vcf && \
            bgzip -c ~{prefix}.vcf > ~{prefix}.vcf.gz && \
            tabix -p vcf ~{prefix}.vcf.gz && \
            find ./ -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g'
    >>>

    output {
        File vcfgz = "~{prefix}.vcf.gz"
        File tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             3,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
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

task IndexVCF {

    meta {
        description: "Indexing vcf.gz. Note: do NOT use remote index as that's buggy."
    }

    input {
        File vcf
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(vcf, ".vcf.gz")
    Int proposed_disk = 3*ceil(size(vcf, "GB")) + 1
    Int disk_size = if (proposed_disk > 100) then proposed_disk else 100

    command <<<
        cp ~{vcf} ~{prefix}.vcf.gz && \
            tabix -p vcf ~{prefix}.vcf.gz && \
            find ./ -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g'
    >>>

    output {
        File tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             3,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
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

task FixSnifflesVCF {
    input {
        File vcf
        String sample_name
        File? ref_fasta_fai
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        sample_name:    "Sniffles infers sample name from the BAM file name, so we fix it here"
        ref_fasta_fai:  "provide only when the contig section of the input vcf is suspected to be corrupted"
    }

    Boolean fix_contigs = defined(ref_fasta_fai)

    Boolean vcf_is_bgzipped = sub(vcf, ".gz", "") != sub(vcf, ".vcf.gz", "")
    String local_raw = if vcf_is_bgzipped then "to.be.fixed.vcf.gz" else "to.be.fixed.vcf"
    String local_sp_fixed = if vcf_is_bgzipped then "sample.fixed.vcf.gz" else "sample.fixed.vcf"

    String initial_grep_cmd = if vcf_is_bgzipped then "zgrep" else "grep"

    String prefix = if vcf_is_bgzipped then basename(vcf, ".vcf.gz") else basename(vcf, ".vcf")
    Int proposed_disk = 3*ceil(size(vcf, "GB")) + 1
    Int disk_size = if (proposed_disk > 100) then proposed_disk else 100

    command <<<
        set -euxo pipefail

        # 1. fix sample information (Sniffles derives VCF SM information from the path to the BAM ......)
        cp ~{vcf} ~{local_raw}
        echo ~{sample_name} > sample_names.txt
        bcftools reheader --samples sample_names.txt -o ~{local_sp_fixed} ~{local_raw}
        rm ~{vcf} && rm ~{local_raw}

        ####################################################################
        # 2. prep for fixing undefined VCF INFO/FT/FORMAT, also guard against when the VCF is empty
        ~{initial_grep_cmd} "^##" ~{local_sp_fixed} > header.txt
        ~{initial_grep_cmd} -v "^#" ~{local_sp_fixed} > body.txt || true
        if [[ ! -f body.txt ]] || [[ ! -s body.txt ]]; then
            echo "input VCF seem to contain only header, but I'll proceed anyway and give you only header"
            bcftools \
                sort \
                --temp-dir tm_sort \
                --output-type z \
                -o ~{prefix}.vcf.gz \
                ~{local_sp_fixed}
            bcftools index --tbi --force ~{prefix}.vcf.gz
            exit 0;
        fi

        ####################################################################
        # 2.1. more prep for fixing undefined VCF INFO/FT/FORMATs
        # get FORMATs in header
        if grep -q -F '##FORMAT=<' header.txt; then
            grep -F '##FORMAT=<' header.txt | awk -F ',' '{print $1}' | sed 's/##FORMAT=<ID=//' | sort > formats.in_header.txt
        else
            touch formats.in_header.txt
        fi
        # get FILTERs in header
        if grep -q -F '##FILTER=<' header.txt; then
            grep -F '##FILTER=<' header.txt | awk -F ',' '{print $1}' | sed 's/##FILTER=<ID=//' | sort > filters.in_header.txt
        else
            touch filters.in_header.txt
        fi
        # get non-flag INFO in header
        if grep -q -F '##INFO=<' header.txt; then
            grep -F '##INFO=<' header.txt | grep -vF 'Type=Flag' | awk -F ',' '{print $1}' | sed 's/##INFO=<ID=//' | sort > non_flag_info.in_header.txt
        else
            touch non_flag_info.in_header.txt
        fi
        # get     flag INFO in header
        if grep -q -F '##INFO=<' header.txt; then
            grep -F '##INFO=<' header.txt | grep  -F 'Type=Flag' | awk -F ',' '{print $1}' | sed 's/##INFO=<ID=//' | sort >     flag_info.in_header.txt
        else
            touch flag_info.in_header.txt
        fi

        # get FORMATs in practice
        awk '{print $9}' body.txt | sort | uniq | sed 's/:/\n/g' | sort | uniq > formats.in_vcf.txt
        # get FILTERs in practice, guard against no 'PASS'
        awk '{print $7}' body.txt | sort | uniq | grep -v "^PASS$" > filters.in_vcf.txt || touch filters.in_vcf.txt

        awk '{print $8}' body.txt | sed 's/;/\n/g' > tmp.info.entries.txt
        if grep -q -F '=' tmp.info.entries.txt; then
            # get non-flag INFO in practicez
            grep -F '=' tmp.info.entries.txt | awk -F '=' '{print $1}' | sort | uniq > non_flag_info.in_vcf.txt
        fi
        if grep -q -vF '=' tmp.info.entries.txt; then
            # get     flag INFO in practice
            awk '{print $8}' body.txt | sed 's/;/\n/g' | grep -vF '=' | sort | uniq > flag_info.in_vcf.txt
        fi
        touch non_flag_info.in_vcf.txt
        touch     flag_info.in_vcf.txt

        echo "I survived. More to go..."

        ####################################################################
        # 2.2. more prep for fixing undefined VCF INFO/FT/FORMATs
        comm -13 formats.in_header.txt formats.in_vcf.txt > missing.formats.txt
        while IFS= read -r line
        do
        echo "##FORMAT=<ID=${line},Number=.,Type=String,Description=\"CALLER DID NOT DEFINE THIS.\">" >> missing.formats.header
        done < missing.formats.txt

        comm -13 filters.in_header.txt filters.in_vcf.txt > missing.filters.txt
        while IFS= read -r line
        do
        echo "##FILTER=<ID=${line},Description=\"CALLER DID NOT DEFINE THIS.\">" >> missing.filters.header
        done < missing.filters.txt

        comm -13 non_flag_info.in_header.txt non_flag_info.in_vcf.txt > missing.non_flag_info.txt
        while IFS= read -r line
        do
        echo "##INFO=<ID=${line},Number=.,Type=String,Description=\"CALLER DID NOT DEFINE THIS.\">" >> missing.non_flag_info.header
        done < missing.non_flag_info.txt

        comm -13 flag_info.in_header.txt flag_info.in_vcf.txt > missing.flag_info.txt
        while IFS= read -r line
        do
        echo "##INFO=<ID=${line},Number=0,Type=Flag,Description=\"CALLER DID NOT DEFINE THIS.\">" >> missing.flag_info.header
        done < missing.flag_info.txt

        ####################################################################
        # 2. actually fix undefined VCF INFO/FT/FORMATs, if necessary
        if  find . -maxdepth 1 -type f -name "missing.*.header" 2>/dev/null | grep -q .; then
            grep "^##" ~{local_sp_fixed} | grep -v "^##[A-Z]" | grep -vF 'contig=' > first_lines.txt
            grep -F "##contig=<ID=" header.txt > contigs.txt
            grep "^#CHROM" ~{local_sp_fixed} > sample.line.txt
            grep "^##" ~{local_sp_fixed} | grep "^##[A-Z]" | sort > existing_definitions.txt
            cat existing_definitions.txt missing.*.header | sort > everything.defined.txt
            cat first_lines.txt contigs.txt everything.defined.txt sample.line.txt > fixed.header.txt
            # print to stdout for checking
            grep -vF "##contig=<ID=" fixed.header.txt

            cat fixed.header.txt body.txt > fixed.vcf
            rm ~{local_sp_fixed}
        else
            mv ~{local_sp_fixed} fixed.vcf
        fi

        ####################################################################
        # 3. fix contigs undefined (in later stages)
        if ~{fix_contigs}; then
            bcftools reheader \
                --fai ~{ref_fasta_fai} \
                -o fixed.and_contigs.vcf \
                fixed.vcf
            mv fixed.and_contigs.vcf fixed.vcf
        fi

        ####################################################################
        # 4. fix occationally unsorted VCF
        bcftools \
            sort \
            --temp-dir tm_sort \
            --output-type z \
            -o ~{prefix}.vcf.gz \
            fixed.vcf
        bcftools index --tbi --force ~{prefix}.vcf.gz
    >>>

    output {
        File sortedVCF = "~{prefix}.vcf.gz"
        File tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             3,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        2,
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

task FillTags {
    input {
        File vcf
        Array[String] tags

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(vcf, "GB"))

    command <<<
        set -euxo pipefail

        bcftools +fill-tags ~{vcf} -Oz -o ~{prefix}.vcf.bgz -- -t ~{sep="," tags}
        tabix -p vcf ~{prefix}.vcf.bgz
    >>>

    output {
        File filled_vcf = "~{prefix}.vcf.bgz"
        File filled_tbi = "~{prefix}.vcf.bgz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-shapeit5:1.0.0-beta"
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

task CountVariants {
    input {
        File vcf
        File tbi
        String locus

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size([vcf, tbi], "GB"))

    command <<<
        set -euxo pipefail

        bcftools view ~{vcf} ~{locus} | grep -c -v '^#' | wc -l > num_variants.txt
    >>>

    output {
        Int num_variants = read_int("num_variants.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
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
