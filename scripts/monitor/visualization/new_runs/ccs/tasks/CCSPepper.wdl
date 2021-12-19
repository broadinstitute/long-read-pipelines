version 1.0

#######################################################
# This pipeline calls small variants using DeepVariant.
#######################################################

import "Structs.wdl"


task SubsetBam {
    input {
        File bam
        File bai

        File interval_list_file
        String interval_id
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: {
            description: "bam to subset",
            localization_optional: true
        }
        bai:    "index for bam file"
        interval_list_file:  "a Picard-style interval list file to subset reads with"
        interval_id:         "an ID string for representing the intervals in the interval list file"
        prefix: "prefix for output bam and bai file names"
    }

    Array[String] intervals = read_lines(interval_list_file)

    Int disk_size = 2*ceil(size([bam, bai], "GB"))

    String subset_prefix = prefix + "." + interval_id

    command <<<
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        
        set -euxo pipefail
        # see man page for what '-M' means
        samtools view \
            -bhX \
            -M \
            -@ 1 \
            --verbosity=8 \
            --write-index \
            -o "~{subset_prefix}.bam##idx##~{subset_prefix}.bam.bai" \
            ~{bam} ~{bai} \
            ~{sep=" " intervals}
    >>>

    output {
        File subset_bam = "~{subset_prefix}.bam"
        File subset_bai = "~{subset_prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.0"
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

task SortVCFs {
    input {
        Array[File] vcfs
        Array[File] ?tbis

        File ref_dict

        String prefix

        Boolean use_ssd = false

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + round(size(vcfs, "GB"))

    command <<<
        set -euxo pipefail

        gatk \
            SortVcf \
            -I ~{sep=' -I ' vcfs} \
            -O ~{prefix}.vcf.gz \
            --SEQUENCE_DICTIONARY ~{ref_dict}
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
        docker:             "us.gcr.io/broad-gatk/gatk:latest"
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

task SortVCFsFast {
    meta {
        description: "Fast merging VCFs when the default sorting is expected to be slow"
    }

    input {
        Array[File] vcfs
        Array[File] tbis

        String prefix

        Boolean sort = true

        RuntimeAttr? runtime_attr_override
    }

    Int sz = ceil(size(vcfs, 'GB'))
    Int machine_memory = if sort then 620 else 96 # sz * 2
    Int work_memory = ceil(machine_memory * 0.8)
    Int cores = if sort then 96 else 16

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
            -o concatedated_raw.vcf  # fast, at the expense of disk space
        for vcf in ~{sep=' ' vcfs}; do rm $vcf ; done
        if ! ~{sort}; then
            echo "==========================================================="
            echo "done concatenation; instructed to NOT sort, so just compress" && date
            echo "==========================================================="
            bgzip -c concatedated_raw.vcf > "~{prefix}.unsorted.vcf.gz"
            echo "==========================================================="
            echo "done compression; instructed to NOT sort, so no index tbi"
            echo "==========================================================="
            exit 0;
        fi
        echo "==========================================================="
        echo "done concatenation, starting sort operation" && date
        echo "==========================================================="
        mkdir -p tm_sort
        bcftools \
            sort \
            -m ~{work_memory}G \
            --temp-dir tm_sort \
            --output-type z \
            -o ~{prefix}.vcf.gz \
            concatedated_raw.vcf
        echo "==========================================================="
        echo "done sorting" && date
        echo "==========================================================="
    >>>

    output {
        File vcf = "~{prefix}" + if sort then "." else ".unsorted." + "vcf.gz"
        File? tbi = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cores,
        mem_gb:             "~{machine_memory}",
        disk_gb:            375,  # just use one local SSD disk, as the known solutions use non-trvial amount of disk for swap
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

task PEPPER {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        Boolean keep_logs_and_intermediates = false

        Int memory

        RuntimeAttr? runtime_attr_override
    }

    Int bam_sz = ceil(size(bam, "GB"))
	Int disk_size = if bam_sz > 200 then 2*bam_sz else bam_sz + 200

    String prefix = basename(bam, ".bam") + ".pepper"
    String output_root = "/cromwell_root/pepper_output"

    command <<<
        set -euxo pipefail

        touch ~{bai}
        #########################################
        # Q: do we really need this?
        SM=$(samtools view -H "~{bam}" | grep -m1 '^@RG' | sed 's/\t/\n/g' | grep '^SM:' | sed 's/SM://g')

        samtools view --no-PG -H "~{bam}" | sed 's/\tPU:.\+//g' > header.txt
        samtools reheader -P header.txt "~{bam}" > fixed.bam
        samtools index fixed.bam
        #########################################

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        mkdir -p "~{output_root}"

        # no gVCF as it Pepper simply doesn't produce gVCF on CCS data
        run_pepper_margin_deepvariant \
            call_variant \
            -b fixed.bam \
            -f ~{ref_fasta} \
            -t "${num_core}" \
            -s "${SM}" \
            -o "~{output_root}" \
            -p "~{prefix}" \
            --phased_output \
            --ccs

        if [[ -f "~{output_root}/intermediate_files/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam" ]]; then
            mv "~{output_root}/intermediate_files/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam" \
               "~{output_root}/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam"
            mv "~{output_root}/intermediate_files/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai" \
               "~{output_root}/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai"
        fi
        find "~{output_root}/" -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g' \
            > "~{output_root}/dir_structure.txt"

        # # we can turn it back on, if so desired
        # if [[ ! ~{keep_logs_and_intermediates} ]]; then
        #     rm -rf \
        #         "~{output_root}/logs" \
        #         "~{output_root}/intermediate_files"
        # fi
    >>>

    output {

        File output_dir_structure = "~{output_root}/dir_structure.txt"

        File phased_bam = "~{output_root}/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam"
        File phased_bai = "~{output_root}/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam.bai"

        # pepper SNP vcf might be less useful, it at all
        # File VCF        = "~{output_root}/~{prefix}.vcf.gz"
        # File VCF_tbi    = "~{output_root}/~{prefix}.vcf.gz.tbi"

        # File phasedVCF  = "~{output_root}/~{prefix}.phased.vcf.gz"
        # File phasedtbi  = "~{output_root}/~{prefix}.phased.vcf.gz.tbi"
        
        # # maybe less useful
        # File phaseset_bed = "~{output_root}/~{prefix}.phaseset.bed"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             memory, 
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "kishwars/pepper_deepvariant:r0.4.1"
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

task DV {

    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        Int memory

        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, ".bam") + ".deepvariant"
    String output_root = "/cromwell_root/dv_output"

    Int bam_sz = ceil(size(bam, "GB"))
	Int disk_size = if bam_sz > 200 then 4*bam_sz else bam_sz + 1000

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        mkdir -p "~{output_root}"

        /opt/deepvariant/bin/run_deepvariant \
            --model_type=PACBIO \
            --ref=~{ref_fasta} \
            --reads=~{bam} \
            --output_vcf="~{output_root}/~{prefix}.vcf.gz" \
            --output_gvcf="~{output_root}/~{prefix}.g.vcf.gz" \
            --num_shards="${num_core}" \
            --use_hp_information
        
        find "~{output_root}/" -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g' \
            > "~{output_root}/dir_structure.txt"
    >>>

    output {

        File output_dir_structure = "~{output_root}/dir_structure.txt"

        File VCF        = "~{output_root}/~{prefix}.vcf.gz"
        File VCF_tbi    = "~{output_root}/~{prefix}.vcf.gz.tbi"

        File gVCF       = "~{output_root}/~{prefix}.g.vcf.gz"
        File gVCF_tbi   = "~{output_root}/~{prefix}.g.vcf.gz.tbi"

        File visual_report_html = "~{output_root}/~{prefix}.visual_report.html"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "google/deepvariant:1.2.0"
        # docker:             "google/deepvariant:1.2.0-gpu"
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
