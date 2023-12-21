version 1.0

import "../../structs/Structs.wdl"

task Quast {

    meta {
        description: "A task that runs QUAST to evaluate a given set of assemblies on a species with existing reference assembly. Entire Quast output will be tarballed"
    }
    parameter_meta {
        ref:        "reference assembly of the species"
        assemblies: "list of assemblies to evaluate"
    }

    input {
        File? ref
        Array[File] assemblies
        Boolean is_large = false

        RuntimeAttr? runtime_attr_override
    }

    Int minimal_disk_size = 2*(ceil(size(ref, "GB") + size(assemblies, "GB")))
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100

    String size_optimization = if is_large then "--large" else " "

    command <<<
        set -eux

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        quast --no-icarus \
              "~{size_optimization}" \
              --threads "${num_core}" \
              ~{true='-r' false='' defined(ref)} \
              ~{select_first([ref, ""])} \
              ~{sep=' ' assemblies}

        tree -h quast_results/

        if [[ -d quast_results/contigs_reports ]]; then
            tar -zcvf contigs_reports.tar.gz quast_results/contigs_reports
        fi
    >>>

    output {
        File report_txt = "quast_results/latest/report.txt"
        File report_html = "quast_results/latest/report.html"

        Array[File] report_in_various_formats = glob("quast_results/latest/report.*")

        Array[File] plots = glob("quast_results/latest/basic_stats/*.pdf")

        File? contigs_reports = "contigs_reports.tar.gz"
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:             16,
        mem_gb:                80,
        disk_gb:               disk_size,
        boot_disk_gb:          10,
        preemptible_tries:     0,
        max_retries:           0,
        docker:                "us.gcr.io/broad-dsp-lrma/lr-quast:5.2.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                   select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:        select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:           select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:            select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                select_first([runtime_attr.docker, default_attr.docker])
    }
}

task SummarizeQuastReport {

    meta {
        description: "A task that summarizes the QUAST report into a single tab-delimited file"
    }
    parameter_meta {
        quast_report_txt: "the QUAST report file"
    }

    input {
        File quast_report_txt
    }

    command <<<
        set -eux
        grep -v -e '^All statistics' -e '^$' ~{quast_report_txt} | \
            sed 's/ /_/g' | \
            sed 's/__\+/\t/g' | \
            sed 's/\s\+$//g' | \
            sed 's/>=/gt/g' | \
            tee quast_summary.txt

        for i in $(seq 2 $(awk '{print NF}' quast_summary.txt | sort -nu | tail -n 1))
        do
            j=$(( i - 2 ))  # to make sure the primary, assuming it's the 0-th fed in to this task and the left-most value column
            cut -d$'\t' -f1,${i} < quast_summary.txt > quast_summary_${j}.txt
        done
    >>>

    output {
        File quast_metrics_together = "quast_summary.txt"
        Array[File] quast_metrics = glob("quast_summary_*.txt")
    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task SummarizeAssemblyMetrics {
    meta {
        desciption:
        "Parse and convert most critical metrics from the QUAST report on the [primary, H1 and H2] assemblies"
    }

    parameter_meta {
        quast_summary: {
            desciption: "produced by task SummarizeQuastReport",
            note: "It should be a Map[String, Int] conceptually, but given that WDL/Cromwell/Terra enforces a Int32 max value on Int types, human genomes sizes overflow that value, hence we use String"
        }
    }

    input {
        File quast_summary
    }

    output {
        Map[String, String] summary = read_map("result.tsv")
    }

    command <<<
        set -eux

        python <<CODE
        import csv

        asm_length = list()
        auN = list()
        num_tigs = list()
        with open("~{quast_summary}") as inf:
            rd = csv.reader(inf, delimiter="\t", quotechar='"')
            for row in rd:
                s = row[0]
                if 'Total_length' == s:
                    asm_length = [int(e) for e in row[1:4]]
                elif 'auN' == s:
                    auN = [round(float(e)) for e in row[1:4]]
                elif '#_contigs' == s:
                    num_tigs = [int(e) for e in row[1:4]]

        with open("result.tsv", 'w') as outf:

            outf.write(f"primary_length\t{asm_length[0]}\n")
            outf.write(f"h1_length\t{asm_length[1]}\n")
            outf.write(f"h2_length\t{asm_length[2]}\n")

            outf.write(f"primary_auN\t{auN[0]}\n")
            outf.write(f"h1_auN\t{auN[1]}\n")
            outf.write(f"h2_auN\t{auN[2]}\n")

            outf.write(f"primary_nTigs\t{num_tigs[0]}\n")
            outf.write(f"h1_nTigs\t{num_tigs[1]}\n")
            outf.write(f"h2_nTigs\t{num_tigs[2]}\n")
        CODE
    >>>

    runtime {
        disks: "local-disk 10 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/python:3.9.18-slim-bullseye"
    }
}
