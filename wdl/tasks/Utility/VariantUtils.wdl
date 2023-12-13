version 1.0

import "../../structs/Structs.wdl"

task MergePerChrCalls {

    meta {
        description: "Merge per-chromosome calls into a single VCF"
    }

    parameter_meta {
        vcfs: "List of per-chromosome VCFs to merge"
        ref_dict: "Reference dictionary"
        prefix: "Prefix for output VCF"
    }

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

    parameter_meta {
        header_definitions_file: "a union of definition header lines for input VCFs (related to https://github.com/samtools/bcftools/issues/1629)"
    }

    input {
        Array[File] vcfs

        File ref_fasta_fai
        File? header_definitions_file

        String prefix
        String? output_bucket

        RuntimeAttr? runtime_attr_override
    }

    Int sz = ceil(size(vcfs, 'GB'))
    Int disk_sz = if sz > 100 then 5 * sz else 375  # it's rare to see such large gVCFs, for now

    Boolean suspected_incomplete_definitions = defined(header_definitions_file)

    Int cores = 8

    # pending a bug fix (bcftools github issue 1576) in official bcftools release,
    # bcftools sort can be more efficient in using memory
    Int machine_memory = 48 # 96

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
        vcfName="~{prefix}.vcf.gz"
        bcftools \
            sort \
            --temp-dir tm_sort \
            --output-type z \
            -o "$vcfName" \
            wgs_raw.vcf.gz
        bcftools index --tbi --force ~{prefix}.vcf.gz
        echo "==========================================================="
        echo "done sorting" && date
        echo "==========================================================="

        tbiName="${vcfName}.tbi"
        if ~{defined(output_bucket)}; then
            outDir="~{sub(select_first([output_bucket]), "/?$", "/")}"
            gcloud storage cp "$vcfName" "$tbiName" "$outDir"
            vcfName="${outDir}$vcfName"
            tbiName="${outDir}$tbiName"
        fi
    >>>

    output {
        File vcf = "$vcfName"
        File tbi = "$tbiName"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cores,
        mem_gb:             "~{machine_memory}",
        disk_gb:            disk_sz,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:         0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
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

    parameter_meta {
        fingerprint_vcf: "Assumed to be genotyped, and hold only one sample (other samples will be ignored)."
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File fingerprint_vcf
        RuntimeAttr? runtime_attr_override
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
        preemptible:     3
        maxRetries:      2
        docker:"gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task SubsetVCF {

    meta {
        description: "Subset a VCF file to a given locus"
    }

    parameter_meta {
        vcf_gz: "VCF file to be subsetted"
        vcf_tbi: "Tabix index for the VCF file"
        locus: "Locus to be subsetted"
        prefix: "Prefix for the output file"
        runtime_attr_override: "Override default runtime attributes"
    }

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

    parameter_meta {
        vcf: "VCF file to be zipped and indexed"
        runtime_attr_override: "Override default runtime attributes"
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
