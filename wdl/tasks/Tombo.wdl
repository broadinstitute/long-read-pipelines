version 1.0

import "Structs.wdl"

workflow Tombo {
    input {
        Array[File] manifest_chunks
        File ref_fasta

        RuntimeAttr? runtime_attr_override
    }

    scatter (manifest in manifest_chunks) {
        call Resquiggle { input: manifest = manifest, ref_fasta = ref_fasta }

        String fast5_dirs = sub(Resquiggle.mod_fast5s[0], basename(Resquiggle.mod_fast5s[0]), "")
    }

    call DetectModifications { input: fast5_dirs = fast5_dirs }

    output {
        File cpg_minus_wig = DetectModifications.cpg_minus_wig
        File cpg_plus_wig  = DetectModifications.cpg_plus_wig
        File dam_minus_wig = DetectModifications.dam_minus_wig
        File dam_plus_wig  = DetectModifications.dam_plus_wig
        File dcm_minus_wig = DetectModifications.dcm_minus_wig
        File dcm_plus_wig  = DetectModifications.dcm_plus_wig

        File cpg_stats = DetectModifications.cpg_stats
        File dam_stats = DetectModifications.dam_stats
        File dcm_stats = DetectModifications.dcm_stats

        File cov_minus_bedgraph = DetectModifications.cov_minus_bedgraph
        File cov_plus_bedgraph  = DetectModifications.cov_plus_bedgraph

        File cpg_sites_pdf = DetectModifications.cpg_sites_pdf
        File dam_sites_pdf = DetectModifications.dam_sites_pdf
        File dcm_sites_pdf = DetectModifications.dcm_sites_pdf
    }
}

task Resquiggle {
    input {
        File manifest
        File ref_fasta

        RuntimeAttr? runtime_attr_override
    }

    Array[File] fast5s = read_lines(manifest)

    Int cpus = 4
    Int disk_size = 10*ceil(size(fast5s, "GB"))

    command <<<
        set -euxo pipefail

        FAST5_DIR=`dirname ~{fast5s[0]}`

        mkdir -p data/single_reads
        multi_to_single_fast5 --input_path ${FAST5_DIR} --save_path data/single_reads --recursive

        tombo resquiggle data/single_reads/ ~{ref_fasta} --processes ~{cpus} --num-most-common-errors 5

        tree -h
    >>>

    output {
        Array[File] mod_fast5s = glob("data/single_reads/*/*.fast5")
        File filename_mapping = "data/single_reads/filename_mapping.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             30,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-ont:0.01.01"
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

task DetectModifications {
    input {
        Array[String] fast5_dirs

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 4
    Int disk_size = 2000

    command <<<
        set -euxo pipefail

        mkdir fast5s
        gsutil -m cp -r ~{sep=' ' fast5_dirs} fast5s/

        tombo detect_modifications alternative_model --fast5-basedirs fast5s/ \
            --statistics-file-basename native.human_sample \
            --alternate-bases CpG dam dcm --processes ~{cpus}

        tombo plot most_significant --fast5-basedirs fast5s/ \
            --statistics-filename native.human_sample.CpG.tombo.stats \
            --plot-standard-model --plot-alternate-model CpG \
            --pdf-filename sample.most_significant_CpG_sites.pdf

        tombo plot most_significant --fast5-basedirs fast5s/ \
            --statistics-filename native.human_sample.dam.tombo.stats \
            --plot-standard-model --plot-alternate-model dam \
            --pdf-filename sample.most_significant_dam_sites.pdf

        tombo plot most_significant --fast5-basedirs fast5s/ \
            --statistics-filename native.human_sample.dcm.tombo.stats \
            --plot-standard-model --plot-alternate-model dcm \
            --pdf-filename sample.most_significant_dcm_sites.pdf

        tombo text_output browser_files --statistics-filename native.human_sample.CpG.tombo.stats \
            --file-types dampened_fraction --browser-file-basename native.human_sample.CpG

        tombo text_output browser_files --statistics-filename native.human_sample.dam.tombo.stats \
            --file-types dampened_fraction --browser-file-basename native.human_sample.dam

        tombo text_output browser_files --statistics-filename native.human_sample.dcm.tombo.stats \
            --file-types dampened_fraction --browser-file-basename native.human_sample.dcm

        tombo text_output browser_files --fast5-basedirs fast5s/ \
            --file-types coverage --browser-file-basename native.human_sample

        tree -h
    >>>

    output {
        File cpg_minus_wig = "native.human_sample.CpG.dampened_fraction_modified_reads.minus.wig"
        File cpg_plus_wig  = "native.human_sample.CpG.dampened_fraction_modified_reads.plus.wig"
        File dam_minus_wig = "native.human_sample.dam.dampened_fraction_modified_reads.minus.wig"
        File dam_plus_wig  = "native.human_sample.dam.dampened_fraction_modified_reads.plus.wig"
        File dcm_minus_wig = "native.human_sample.dcm.dampened_fraction_modified_reads.minus.wig"
        File dcm_plus_wig  = "native.human_sample.dcm.dampened_fraction_modified_reads.plus.wig"

        File cpg_stats = "native.human_sample.CpG.tombo.stats"
        File dam_stats = "native.human_sample.dam.tombo.stats"
        File dcm_stats = "native.human_sample.dcm.tombo.stats"

        File cov_minus_bedgraph = "native.human_sample.coverage.minus.bedgraph"
        File cov_plus_bedgraph  = "native.human_sample.coverage.plus.bedgraph"

        File cpg_sites_pdf = "sample.most_significant_CpG_sites.pdf"
        File dam_sites_pdf = "sample.most_significant_dam_sites.pdf"
        File dcm_sites_pdf = "sample.most_significant_dcm_sites.pdf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/broad-long-read-pipelines/lr-ont:0.01.01"
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

