version 1.0

import "../../structs/Structs.wdl"

task CreateMOBsuiteDB {
    input {
        File additional_plasmids_fasta
        File plasmid_host_tsv

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             16,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us-central1-docker.pkg.dev/broad-dsp-lrma/general/mobsuite:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Int num_cores = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    command <<<
        set -euxo pipefail

        mkdir default_db
        mob_init -d default_db

        cp ~{additional_plasmids_fasta} .
        mob_typer -d default_db -i ~{basename(additional_plasmids_fasta)} --multi -o typing_results.tsv

        mob_cluster -d default_db \
            -f ~{basename(additional_plasmids_fasta)} \
            -t ~{plasmid_host_tsv} \
            -p typing_results.tsv \
            -c default_db/clusters.txt \
            -r default_db/ncbi_plasmid_full_seqs.fas \
            -m update \
            -o updated_db

        # make it look like a standard MOB-suite DB
        cp default_db/host_range_literature_plasmidDB.txt updated_db/
        cp default_db/*.proteins.faa updated_db/
        cp default_db/orit.fas updated_db/
        cp default_db/rep.dna.fas updated_db/
        cp default_db/repetitive* updated_db/
        cp default_db/taxa* updated_db/

        # MOB-suite copies the FASTA from the original database without adding the plasmid contigs from our supplemental
        # plasmid file, so we add them manually and recreate the BLAST DB. To make it look like the default MOB-suite,
        # we also name it "ncbi_plasmid_full_seqs.fas"
        rm updated_db/references_updated.fasta*
        cat default_db/ncbi_plasmid_full_seqs.fas ~{basename(additional_plasmids_fasta)} > updated_db/ncbi_plasmid_full_seqs.fas

        cd updated_db
        mash sketch ncbi_plasmid_full_seqs.fas -o ncbi_plasmid_full_seqs.fas.msh -p ~{num_cores} -k 21 -s 1000
        makeblastdb -in ncbi_plasmid_full_seqs.fas -dbtype nucl
        cd ..

        tar -acf mobsuite_updated_db.tar.gz -C updated_db/ .
    >>>

    output {
        File db = "mobsuite_updated_db.tar.gz"
    }

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

task MOBRecon {
    input {
        File assembly_fasta
        Boolean unicycler_assembly = false
        File? MOBsuite_db

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us-central1-docker.pkg.dev/broad-dsp-lrma/general/mobsuite:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String db_flag = if defined(MOBsuite_db) then "-d mobsuite_db" else ""
    Int num_cores = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    command <<<
        set -euxo pipefail

        if [[ "~{MOBsuite_db}" != "" ]]; then
            mkdir mobsuite_db
            cd mobsuite_db
            tar -xaf ~{MOBsuite_db} .
            cd ..
        fi

        mob_recon -i ~{assembly_fasta} -n ~{num_cores} \
            ~{true="-u" false="" unicycler_assembly} ~{db_flag} \
            -o results
    >>>

    output {
        File chromosome = "results/chromosome.fasta"
        File contig_report = "results/contig_report.txt"
        File mge_report = "results/mge.report.txt"
        File mobtyper_results = "results/mobtyper_results.txt"
        Array[File] plasmids = glob("results/plasmid_*.fasta")
    }

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

task GFFFilterChromosomal {
    input {
        File base_gff3
        File mobsuite_contig_report

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us-central1-docker.pkg.dev/broad-dsp-lrma/general/mobsuite:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euxo pipefail

        fname="~{basename(base_gff3, ".gff3")}.plasmids.gff3"
        gff_filter_chromosome.py ~{base_gff3} ~{mobsuite_contig_report} > $fname
    >>>

    output {
        File plasmids_gff3 = "~{basename(base_gff3, '.gff3')}.plasmids.gff3"
    }

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
