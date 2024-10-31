version 1.0

# Given a PED file and a set of methylation bed files from pb-cpg-tools, 
# find differentially methylated regions for each affected sample vs. the pool of unaffected samples
# Filter cpg sites for promoter and enhancer regions
# and shard by chromosome

workflow DifferentialMethylation {
    input {
        File pedFile
        File? regions_filter
        Array[File] methylationBedFiles
        Array[String] contigs
        String inputType = "PB" # "PB" or "ONT"
    }

    scatter (contig in contigs) {
        call ShardByChromosome {
            input:
                methylationBedFiles = methylationBedFiles,
                contig = contig
        }

        if (defined(regions_filter)) {
            call FilterBed {
                input:
                    bedFile_tarball = ShardByChromosome.MethylationBedsByChr,
                    regions_filter = select_first([regions_filter])
            }
        }

        call DMR {
            input:
                pedFile = pedFile,
                methylationBedFiles = select_first([FilterBed.filtered_beds, ShardByChromosome.MethylationBedsByChr]),
                inputType = inputType
        }
    }

    # merge DMR results
    call MergeDMRs {
        input:
            tsvs = DMR.dmr_tsv
    }

    call MergePDFs {
        input:
            pdfs = DMR.dmr_pdfs
    }

    call MergeRDS {
        input:
            rds_files = DMR.smoothed_rds
    }

    output {
        File combined_dmr_tsv = MergeDMRs.combined_dmr_tsv
        File combined_dmr_pdfs = MergePDFs.combined_dmr_pdfs
        File smoothed_rds = MergeRDS.combined_smoothed_rds
    }
}

task MergePDFs {
    input {
        Array[File] pdfs

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 1,
        disk_gb: ceil(size(pdfs, "GB")*10) + 5,
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10,
        docker: "quay.io/ymostovoy/lr-utils-basic"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         runtime_default.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      runtime_default.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       runtime_default.max_retries])
        docker:                 select_first([runtime_attr.docker,            runtime_default.docker])
    }

    command <<<
        for tarred_pdf in ~{sep=" " pdfs}; do
            tar xzvf $tarred_pdf # extract each pdf to the current directory
        done

        tar czvf combined_dmrs_pdfs.tar.gz *.pdf
    >>>

    output {
        File combined_dmr_pdfs = "combined_dmrs_pdfs.tar.gz"
    }
}


task MergeRDS {
    input {
        Array[File] rds_files

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 1,
        disk_gb: ceil(size(rds_files, "GB")*10) + 5,
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10,
        docker: "quay.io/ymostovoy/lr-utils-basic"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         runtime_default.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      runtime_default.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       runtime_default.max_retries])
        docker:                 select_first([runtime_attr.docker,            runtime_default.docker])
    }

    command <<<
        for rds in ~{sep=" " rds_files}; do
            mv $rds ./
        done

        # make tarball of all of the rds files
        tar czvf smoothed.rds.tar.gz *.rds

    >>>

    output {
        File combined_smoothed_rds = "smoothed.rds.tar.gz"
    }
}

task MergeDMRs {
    input {
        Array[File] tsvs

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 1,
        disk_gb: ceil(size(tsvs, "GB")*2) + 5,
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10,
        docker: "quay.io/ymostovoy/lr-utils-basic"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         runtime_default.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      runtime_default.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       runtime_default.max_retries])
        docker:                 select_first([runtime_attr.docker,            runtime_default.docker])
    }

    command <<<
        # get header from first tsv
        head -n 1 ~{tsvs[0]} > dmrs.tsv

        # get all the rows from all the tsvs
        for tsv in ~{sep=" " tsvs}; do
            tail -n +2 $tsv >> tmp
        done
        sort -k1,1 -k2,2n tmp >> dmrs.tsv
    >>>

    output {
        File combined_dmr_tsv = "dmrs.tsv"
    } 
}

task DMR {
    input {
        File pedFile
        File methylationBedFiles
        String inputType

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 32,
        disk_gb: ceil(size(methylationBedFiles, "GB") + size(pedFile, "GB"))*5 + 20,
        cpu_cores: 4,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10,
        docker: "quay.io/ymostovoy/lr-differential-methylation"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         runtime_default.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      runtime_default.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       runtime_default.max_retries])
        docker:                 select_first([runtime_attr.docker,            runtime_default.docker])
    }

    Int num_cores = select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])

    String contig = basename(methylationBedFiles, "_methylation_beds.tar.gz")

    command <<<
        # convert all the methylation bed files to bismark format and put the outputs in the same directory
        mkdir methylation

        tar zxvf ~{methylationBedFiles} # extract the tarball
        # naming format from ShardByChromosome is ${name}.${contig}.bed
        for singlebedfile in *.bed; do
            sample=$(basename ${singlebedfile} .bed | cut -d'.' -f1)
            if [[ ~{inputType} == "PB" ]]; then
                awk 'OFS="\t" {print $1,$2,$3,100*$7/$6,$7,$8}' ${singlebedfile} > methylation/${sample}.bismark.cov
            elif [[ ~{inputType} == "ONT" ]]; then
                awk 'OFS="\t" {print $1,$2,$3,$11,$12,$13}' ${singlebedfile} > methylation/${sample}.bismark.cov
            else
                echo "inputType must be PB or ONT"
                exit 1
            fi
        done
        
        current_dir=$(pwd)

        # find differentially methylated regions
        # purposefully NOT detecting number of cores on the machine because we don't want to use too many (causes memory issues)
        Rscript /callDMRs.R -d methylation/ -p ~{pedFile} -t ~{num_cores} -o $current_dir

        # compress the pdfs from callDMRs.R
        tar czvf ~{contig}_dmr_pdfs.tar.gz *.pdf
        mv dmrs.tsv ~{contig}_dmrs.tsv
        mv smoothed_bsobj.rds ~{contig}_smoothed.rds
    >>>

    output {
        File dmr_tsv = "~{contig}_dmrs.tsv"
        File dmr_pdfs = "~{contig}_dmr_pdfs.tar.gz"
        File smoothed_rds = "~{contig}_smoothed.rds"
    }
}

task FilterBed {
    input {
        File bedFile_tarball
        File regions_filter
    
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 2,
        disk_gb: ceil(size(bedFile_tarball, "GB"))*10 + ceil(size(regions_filter, "GB")) + 20,
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10,
        docker: "quay.io/ymostovoy/lr-utils-basic"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         runtime_default.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      runtime_default.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       runtime_default.max_retries])
        docker:                 select_first([runtime_attr.docker,            runtime_default.docker])
    }

    String base = basename(bedFile_tarball, ".tar.gz")

    command <<<
        set -euxo pipefail
        tar zxvf ~{bedFile_tarball}

        # filter the bed files
        for bedFile in *.bed; do
            name=$(basename ${bedFile} .bed)
            bedtools intersect -u -a $bedFile -b ~{regions_filter} > ${name}.filtered.bed
        done

        tar czvf ~{base}.filtered.tar.gz *.filtered.bed

    >>>

    output {
        File filtered_beds = "~{base}.filtered.tar.gz"
    }
}

task ShardByChromosome {
    input {
        Array[File] methylationBedFiles
        String contig

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 2,
        disk_gb: ceil(size(methylationBedFiles, "GB"))*5 + 20,
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10,
        docker: "quay.io/ymostovoy/lr-utils-basic"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         runtime_default.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           runtime_default.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      runtime_default.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       runtime_default.max_retries])
        docker:                 select_first([runtime_attr.docker,            runtime_default.docker])
    }

    command <<<
        # filter by chromosome
        for methylationBedFile in ~{sep=" " methylationBedFiles}; do
            echo "filtering "$methylationBedFile
            if [[ ${methylationBedFile} == *.gz ]]; then
                name=$(basename ${methylationBedFile} .bed.gz)
                zcat ${methylationBedFile} | awk -v chr=~{contig} '$1 == chr' > ${name}.~{contig}.bed
            else
                name=$(basename ${methylationBedFile} .bed)
                awk -v chr=~{contig} '$1 == chr' ${methylationBedFile} > ${name}.~{contig}.bed
            fi
        done

        echo "### compressing into tarball ###"
        tar czvf ~{contig}_methylation_beds.tar.gz *.~{contig}.bed
    >>>

    output {
        File MethylationBedsByChr = "~{contig}_methylation_beds.tar.gz"
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}