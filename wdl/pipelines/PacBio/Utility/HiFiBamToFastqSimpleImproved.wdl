version 1.0

workflow HiFiBamToFastqSimpleImproved {
    input {
        File bam
        File outdir
    }

    call BamToFastq { input: bam = bam, prefix=basename(bam, ".bam")}
    call FinalizeToFile as FinalizeFq { input:  outdir = outdir, file = BamToFastq.reads_fq, }

    output {
        File hifi_fq = FinalizeFq.gcs_path
    }
}

task BamToFastq {
    meta {
        description : "Convert a long reads BAM file to a fastq file."
        warn: "Please do not include 'RG' in tags_to_preserve, as that's automatically saved"
    }

    parameter_meta {
        bam: {localization_optional: true}
        prefix: "Prefix for the output fastq file."

        save_all_tags:
        "if true, saves all SAM tags to the FASTQ output; cannot set this true while also specifying tags_to_preserve "
        tags_to_preserve:
        "custom list of tags to preserve; please do not include 'RG' in tags_to_preserve, as that's automatically preserved"

        disk_type: "type of disk to use"
    }

    input {
        File bam
        String prefix

        Boolean save_all_tags = false
        Array[String] tags_to_preserve = []

        String disk_type = "SSD"
        RuntimeAttr? runtime_attr_override
    }

    output {
        File reads_fq = "~{prefix}.fq.gz"
    }

    Boolean custom_tags_to_preserve = 0<length(tags_to_preserve)

    String base = basename(bam)
    String local_bam = "/cromwell_root/~{base}"
    command <<<
        set -euxo pipefail

        # some checks on inputs
        if ~{custom_tags_to_preserve} && ~{save_all_tags} ; then
            echo "cannot ask to save all tags and yet also ask to save a custom list of tags" && exit 1;
        fi
        for tag in ~{sep=' ' tags_to_preserve}; do
            if [[ ${tag} == "RG" ]]; then
                echo "RG tag is automatically saved for you, no need to save it" && exit 1
            fi
        done

        # localize
        time \
        gcloud storage cp ~{bam} ~{local_bam}


        # when saving all tags, the list can be empty as instructed by samtools doc
        time \
        samtools fastq \
            -@1 \
            ~{true='-t ' false =' ' save_all_tags} \
            ~{true='-T ' false =' ' custom_tags_to_preserve} ~{sep=',' tags_to_preserve} \
            -0 ~{prefix}.fq.gz \
            ~{local_bam}
    >>>

    #########################
    Int disk_size = 10 + 3 * ceil(size(bam, "GiB"))

    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FinalizeToFile {

    meta{
        description: "Copies the given file to the specified bucket."
    }

    parameter_meta {
        file: {
            description: "file to finalize",
            localization_optional: true
        }
        keyfile : "[optional] File used to key this finaliation.  Finalization will not take place until the KeyFile exists.  This can be used to force the finaliation to wait until a certain point in a workflow.  NOTE: The latest WDL development spec includes the `after` keyword which will obviate this."
        outdir: "directory to which files should be uploaded"
        name:   "name to set for uploaded file"
    }

    input {
        File file
        String outdir
        String? name

        File? keyfile

        RuntimeAttr? runtime_attr_override
    }

    String gcs_output_dir = sub(outdir, "/+$", "")
    String gcs_output_file = gcs_output_dir + "/" + select_first([name, basename(file)])

    command <<<
        set -euxo pipefail

        gcloud storage cp "~{file}" "~{gcs_output_file}"
    >>>

    output {
        String gcs_path = gcs_output_file
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            10,
        preemptible_tries:  2,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
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
