version 1.0

import "../../structs/Structs.wdl"

task Panaroo {
    input {
        File input_manifest
        String clean_mode = "sensitive"
        Boolean clean_edges = false

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            200,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        2,
        docker:             "us-central1-docker.pkg.dev/broad-dsp-lrma/fusilli/panaroo:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    Int num_cpu = select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])

    command <<<
        set -euxo pipefail

        mkdir gffs
        tail -n +2 ~{input_manifest} | awk -F$'\t' '{print $2}' > gcp_gffs.txt
        gsutil -m cp -I gffs < gcp_gffs.txt

        ls -1 gffs/* > local_gffs.txt

        mkdir output
        panaroo -i local_gffs.txt ~{true="" false="--no_clean_edges" clean_edges} --clean-mode ~{clean_mode} \
            -t ~{num_cpu} --remove-invalid-genes -o output/
    >>>

    output {
        File summary_stats = "output/summary_statistics.txt"
        File final_graph = "output/final_graph.gml"
        File pre_filt_graph = "output/pre_filt_graph.gml"
        File pangenome_reference = "output/pan_genome_reference.fa"

        File gene_data = "output/gene_data.csv"
        File gene_presence_absence_rtab = "output/gene_presence_absence.Rtab"
        File gene_presence_absence_csv = "output/gene_presence_absence.csv"
        File gene_presence_absence_roary = "output/gene_presence_absence_roary.csv"

        File combined_DNA_CDS = "output/combined_DNA_CDS.fasta"
        File combined_protein_CDS = "output/combined_protein_CDS.fasta"
        File combined_protein_cdhit_out = "output/combined_protein_cdhit_out.txt"
        File combined_protein_cdhit_out_cluster = "output/combined_protein_cdhit_out.txt.clstr"

        File sv_triplets_presence_absence = "output/struct_presence_absence.Rtab"
    }

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
