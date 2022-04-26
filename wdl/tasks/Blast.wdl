version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    String? disk_type
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

workflow BlastN {
    meta {
        description : "Run BLASTN on a fasta file."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }
    input {
        File input_fasta
        String db_name = "refseq_rna"

        File blast_db = "gs://broad-dsde-methods-long-reads-public/resources/BLAST_DB_SENTINEL"
        Int max_target_seqs = 10
        Int out_format = 6

        RuntimeAttr? runtime_attr_override
    }

    # Call our alignment task:
    call BlastN {
        input:
            input_fasta = input_fasta,
            db_name = db_name,
            max_target_seqs = max_target_seqs,
            out_format = out_format,
            blast_db = blast_db,
            runtime_attr_override = runtime_attr_override
    }

    output {
        File blastn_out = BlastN.out
    }
}

task BlastN {

    meta {
        description : "Run BLASTN on a fasta file."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File input_fasta

        String db_name = "nt"

        File blast_db = "gs://broad-dsde-methods-long-reads-public/resources/BLAST_DB_SENTINEL"
        Int max_target_seqs = 10
        Int out_format = 6

        String prefix = "blastn"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 20 + 5*ceil(size(input_fasta, "GB"))

    command <<<
        set -euxo pipefail

        NUM_CPUS=$(cat /proc/cpuinfo | grep processor | tail -n1 | awk '{print $NF+1}')

        # If we have a reference disk for this blast db, we have a sym link.
        # We'll need to process the given input to make sure we can actually point at
        # the DB properly.
        if [[ -L ~{blast_db} ]] ; then
            # Get the mount directory of the file:
            # Assume that the mount dir is of the form /<MOUNT_POINT>/<DIR>
            export BLASTDB=$(readlink ~{blast_db} | awk 'BEGIN{FS="/"}{printf("/%s/%s", $2,$3)}')
        else
            export BLASTDB=~{blast_db}
        fi

        time blastn \
            -num_threads ${NUM_CPUS} \
            -max_target_seqs ~{max_target_seqs} \
            -query ~{input_fasta} \
            -db "~{db_name}" \
            -outfmt ~{out_format} \
            -out ~{prefix}.out
    >>>

    output {
        File out = "~{prefix}.out"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             32,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-blast:0.0.1"
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