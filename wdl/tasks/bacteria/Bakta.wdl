version 1.0

import "../Structs.wdl"

task BaktaDBDownload {
    input {
        RuntimeAttr? runtime_attr_override
        String filename
    }

    command <<<
        set -euxo pipefail

        bakta_db download --output .

        tar -caf ~{filename} -C db .
    >>>

    output {
        File bakta_db = "~{filename}"
    }


    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            200,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        2,
        docker:             "quay.io/biocontainers/bakta:1.6.1--pyhdfd78af_0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


task BaktaAnnotate {
    input {
        File bakta_db_tar
        File genome_fasta
        String? fname_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int num_cores = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    String prefix = select_first([fname_prefix, basename(genome_fasta)])

    command <<<
        set -euxo pipefail

        mkdir bakta_db
        tar -xaf ~{bakta_db_tar} -C bakta_db

        mkdir output
        bakta --db bakta_db --output output --complete --threads ~{num_cores} \
            --keep-contig-headers --prefix ~{prefix} --verbose \
            ~{genome_fasta}
    >>>

    output {
        File tsv = "output/~{prefix}.tsv"
        File json = "output/~{prefix}.json"
        File gff = "output/~{prefix}.gff3"
        File genbank = "output/~{prefix}.gbff"
        File embl = "output/~{prefix}.embl"
        File ffn = "output/~{prefix}.ffn"
        File faa = "output/~{prefix}.faa"
        File hypotheticals_tsv = "output/~{prefix}.hypotheticals.tsv"
        File hypotheticals_faa = "output/~{prefix}.hypotheticals.faa"

        File summary = "output/~{prefix}.txt"
        File log = "output/~{prefix}.log"
        File plot_png = "output/~{prefix}.png"
        File plot_svg = "output/~{prefix}.svg"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "quay.io/biocontainers/bakta:1.6.1--pyhdfd78af_0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task BaktaAnnotateBatch {
    input {
        File bakta_db_tar
        String output_dir
        File manifest_tsv

        Int worker
        Int batch_size

        RuntimeAttr? runtime_attr_override
    }

    Int num_cores = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    String gcs_output_dir = sub(output_dir, "/+$", "")

    command <<<
        set -euxo pipefail

        mkdir bakta_db
        >&2 echo $(date --rfc-3339=seconds)
        tar -xaf ~{bakta_db_tar} -C bakta_db
        >&2 echo $(date --rfc-3339=seconds)

        lines_start=$(( ~{worker} * ~{batch_size} + 2 ))  # Skip TSV header and sed line numbers start at 1
        lines_end=$(( lines_start + ~{batch_size} - 1 ))
        lines_quit=$(( lines_end + 1 ))
        2>&1 echo "Processing lines ${lines_start}-${lines_end}"

        # Extract plasmid IDs and FASTA filenames to process
        sed -n "${lines_start},${lines_end}p;${lines_quit}q" ~{manifest_tsv} > to_process.tsv

        # Download plasmid FASTA files
        awk -F$'\t' '{print $1}' to_process.tsv > batch_ids.txt
        awk -F$'\t' '{print $2}' to_process.tsv > batch_fasta_gcs.txt
        mkdir -p fasta_files
        gsutil -m cp -I ./fasta_files < batch_fasta_gcs.txt
        for gcs_fname in $(< batch_fasta_gcs.txt); do echo "./fasta_files/${gcs_fname##*/}"; done > batch_fasta.txt

        # List all expected output files for this batch
        awk -F$'\t' '{print "~{gcs_output_dir}/" $1 "/" $1 ".tsv"}' to_process.tsv > batch_tsv.txt
        awk -F$'\t' '{print "~{gcs_output_dir}/" $1 "/" $1 ".json"}' to_process.tsv > batch_json.txt
        awk -F$'\t' '{print "~{gcs_output_dir}/" $1 "/" $1 ".gff3"}' to_process.tsv > batch_gff.txt
        awk -F$'\t' '{print "~{gcs_output_dir}/" $1 "/" $1 ".gbff"}' to_process.tsv > batch_genbank.txt
        awk -F$'\t' '{print "~{gcs_output_dir}/" $1 "/" $1 ".embl"}' to_process.tsv > batch_embl.txt
        awk -F$'\t' '{print "~{gcs_output_dir}/" $1 "/" $1 ".ffn"}' to_process.tsv > batch_ffn.txt
        awk -F$'\t' '{print "~{gcs_output_dir}/" $1 "/" $1 ".faa"}' to_process.tsv > batch_faa.txt
        awk -F$'\t' '{print "~{gcs_output_dir}/" $1 "/" $1 ".hypotheticals.tsv"}' to_process.tsv > batch_hypotheticals_tsv.txt
        awk -F$'\t' '{print "~{gcs_output_dir}/" $1 "/" $1 ".hypotheticals.faa"}' to_process.tsv > batch_hypotheticals_faa.txt
        awk -F$'\t' '{print "~{gcs_output_dir}/" $1 "/" $1 ".txt"}' to_process.tsv > batch_summaries.txt
        awk -F$'\t' '{print "~{gcs_output_dir}/" $1 "/" $1 ".log"}' to_process.tsv > batch_log.txt
        awk -F$'\t' '{print "~{gcs_output_dir}/" $1 "/" $1 ".png"}' to_process.tsv > batch_png.txt
        awk -F$'\t' '{print "~{gcs_output_dir}/" $1 "/" $1 ".svg"}' to_process.tsv > batch_svg.txt

        while IFS=$'\t' read -r plasmid_id fasta tsv json gff genbank embl ffn faa hypotheticals_tsv hypotheticals_faa summary log plot_png plot_svg; do
            >&2 echo $(date --rfc-3339=seconds)

            # Check each output file to see if this genome has already been processed
            # This ensures we can continue from our last genome when the VM gets pre-empted by Google
            if gsutil -q stat "${tsv}" "${json}" "${gff}" "${genbank}" "${embl}" "${ffn}" "${faa}" "${hypotheticals_tsv}" \
                    "${hypotheticals_faa}" "${summary}" "${log}" "${plot_png}" "${plot_svg}"; then
                >&2 echo "Already processed ${plasmid_id}"
            else
                mkdir -p "output/${plasmid_id}"
                bakta --db bakta_db --output "output/${plasmid_id}" --complete --threads ~{num_cores} \
                    --keep-contig-headers --prefix ${plasmid_id} --verbose ${fasta} || true 2>&1 | tee bakta_out.txt

                if grep -q "Exception: diamond error! error code: -9" bakta_out.txt; then
                    >&2 echo "Diamond was Killed. Likely OutOfMemory."  # Trigger Terra retry with more memory
                    exit 1
                elif grep -q "Exception: PILER-CR error! error code: -11" bakta_out.txt; then
                    >&2 echo "PILER-CR crashed. Rerunning bakta without CRISPR predictions."
                    bakta --db bakta_db --output "output/${plasmid_id}" --complete --threads ~{num_cores} --skip-crispr \
                        --keep-contig-headers --prefix ${plasmid_id} --verbose ${fasta} 2>&1 | tee bakta_out.txt
                fi

                gsutil -m cp "output/${plasmid_id}/"* "~{gcs_output_dir}/${plasmid_id}"
                rm -rf "output/${plasmid_id}"
            fi
        done < <(paste batch_ids.txt batch_fasta.txt batch_tsv.txt batch_json.txt batch_gff.txt batch_genbank.txt \
            batch_embl.txt batch_ffn.txt batch_faa.txt batch_hypotheticals_tsv.txt batch_hypotheticals_faa.txt batch_summaries.txt \
            batch_log.txt batch_png.txt batch_svg.txt)
    >>>

    output {
        Array[String] plasmid_ids = read_lines("batch_ids.txt")
        Array[String] tsv = read_lines("batch_tsv.txt")
        Array[String] json = read_lines("batch_json.txt")
        Array[String] gff = read_lines("batch_gff.txt")
        Array[String] genbank = read_lines("batch_genbank.txt")
        Array[String] embl = read_lines("batch_embl.txt")
        Array[String] ffn = read_lines("batch_ffn.txt")
        Array[String] faa = read_lines("batch_faa.txt")
        Array[String] hypotheticals_tsv = read_lines("batch_hypotheticals_tsv.txt")
        Array[String] hypotheticals_faa = read_lines("batch_hypotheticals_faa.txt")

        Array[String] summary = read_lines("batch_summaries.txt")
        Array[String] log = read_lines("batch_log.txt")
        Array[String] plot_png = read_lines("batch_png.txt")
        Array[String] plot_svg = read_lines("batch_svg.txt")
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            150,
        boot_disk_gb:       10,
        preemptible_tries:  50,
        max_retries:        2,
        docker:             "us-central1-docker.pkg.dev/broad-dsp-lrma/fusilli/bakta:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CreateTerraDataTSV {
    input {
        Array[String] plasmid_ids

        Array[String] tsv
        Array[String] json
        Array[String] gff
        Array[String] genbank
        Array[String] embl
        Array[String] ffn
        Array[String] faa
        Array[String] hypotheticals_tsv
        Array[String] hypotheticals_faa

        Array[String] summary
        Array[String] log
        Array[String] plot_png
        Array[String] plot_svg

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        echo $'entity:plasmid_id\tannot_tsv\tannot_json\tannot_gff3\tannot_genbank\tannot_embl\tcds_ffn\tprotein_faa\thypotheticals_tsv\thypotheticals_faa\tannot_summary\tannot_log\tplot_png\tplot_svg' > plasmid_annot.tsv

        paste ~{write_lines(plasmid_ids)} ~{write_lines(tsv)} ~{write_lines(json)} ~{write_lines(gff)} ~{write_lines(genbank)} \
            ~{write_lines(embl)} ~{write_lines(ffn)} ~{write_lines(faa)} ~{write_lines(hypotheticals_tsv)} ~{write_lines(hypotheticals_faa)} \
            ~{write_lines(summary)} ~{write_lines(log)} ~{write_lines(plot_png)} ~{write_lines(plot_svg)} >> plasmid_annot.tsv
    >>>

    output {
        File terra_tsv = "plasmid_annot.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us-central1-docker.pkg.dev/broad-dsp-lrma/fusilli/bakta:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
