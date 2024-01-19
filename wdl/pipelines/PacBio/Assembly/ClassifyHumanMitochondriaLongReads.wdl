version 1.0

workflow ClassifyHumanMitochondriaLongReads {
    meta {
        description:
        "A preprocessing workfow for classifying chrM-mapping reads in a human sample's long reads BAM."
    }
    parameter_meta {
        classifier_docker_tag:
        "tag of docker 'us.gcr.io/broad-dsp-lrma/mito_ccs_read_preprocessor' to use"

        for_assembly_bam :
        "mitochondria reads ready for assembly, in BAM format."
        for_assembly_fasta :
        "mitochondria reads ready for assembly, in FASTA(.gz) format."
        reads_multimap :
        "BAM holding alignment records of sequences that has multiple mapping records."
        reads_too_short :
        "BAM holding alignment records of sequences that are too short (currently hardcoded length of 7kbp)."
        reads_selfoverlap_on_ref :
        "BAM holding alignment records of sequences whose alignment overlap 'significantly' on rCRS reference."
        reads_break_cocircularity :
        "BAM holding alignment records of sequences whose alignments break co-circularity."
        read_length_histogram :
        "A flatfile holding the length distribution of input bam."
        output_notebook :
        "A jupyternotebook display information for debugging how reads are classied."
        read_count_by_class:
        "Count of reads in each category. If they are "
    }
    input {
        File bam
        File bai

        String classifier_docker_tag # 0.0.1 is the latest as of 2024-01-11
    }
    output {
        File for_assembly_bam           = select_first([ClassifyMitoReads.for_assembly_bam])
        File for_assembly_bai           = select_first([ClassifyMitoReads.for_assembly_bai])
        File for_assembly_fasta         = select_first([ClassifyMitoReads.for_assembly_fasta])
        Map[String, Int]? read_count_by_class = CountReads.count_by_case
        File? reads_multimap             = ClassifyMitoReads.reads_multimap
        File? reads_too_short            = ClassifyMitoReads.reads_too_short
        File? reads_selfoverlap_on_ref   = ClassifyMitoReads.reads_selfoverlap_on_ref
        File? reads_break_cocircularity  = ClassifyMitoReads.reads_break_cocircularity
        File? read_length_histogram      = ClassifyMitoReads.read_length_histogram
        File output_notebook            = select_first([ClassifyMitoReads.output_notebook])
    }

    call ExtractMitochondriaReads { input: bam = bam, bai = bai }

    if (ExtractMitochondriaReads.is_samtools_failed) {
        call StopWorkflow { input:
            reason = "Streaming from bucket to subset to chrM failed. This is likely transient. Please rerun without call-caching."
        }
    }

    if (!ExtractMitochondriaReads.is_samtools_failed) {
        call ClassifyMitoReads { input:
            bam = ExtractMitochondriaReads.subset_bam,
            bai = ExtractMitochondriaReads.subset_bai,
            classifier_docker_tag = classifier_docker_tag
        }

        call CountReads { input:
            chrM_bam = ExtractMitochondriaReads.subset_bam,
            classified_bams = select_all([ClassifyMitoReads.for_assembly_bam,
                                          ClassifyMitoReads.reads_multimap, ClassifyMitoReads.reads_too_short,
                                          ClassifyMitoReads.reads_selfoverlap_on_ref, ClassifyMitoReads.reads_break_cocircularity])
        }


    }
}
task ClassifyMitoReads {
    meta {
        description:
        "Classify chrM-mapping reads in a human sample's long reads BAM."
    }
    parameter_meta {
        bam:
        "chrM mapping long reads of a human sample"

        classifier_docker_tag:
        "tag of docker 'us.gcr.io/broad-dsp-lrma/mito_ccs_read_preprocessor' to use"
    }
    input {
        File bam
        File bai

        String classifier_docker_tag
    }
    output {
        File for_assembly_bam           = "local_out_dir/~{prefix}.classified_reads.for_assembly.bam"
        File for_assembly_bai           = "local_out_dir/~{prefix}.classified_reads.for_assembly.bam.bai"
        File for_assembly_fasta         = "local_out_dir/~{prefix}.classified_reads.for_assembly.fastq.gz"

        File? reads_multimap             = "local_out_dir/~{prefix}.classified_reads.multimap.bam"
        File? reads_too_short            = "local_out_dir/~{prefix}.classified_reads.too_short.bam"
        File? reads_selfoverlap_on_ref   = "local_out_dir/~{prefix}.classified_reads.ovp_pair_alns_lr.bam"
        File? reads_break_cocircularity  = "local_out_dir/~{prefix}.classified_reads.break_cocircularity.bam"

        File read_length_histogram      = "local_out_dir/~{prefix}.lengths.histogram.flatfile"
        File output_notebook            = "local_out_dir/~{prefix}.classified_reads.mito_reads_analyzer.papermill_out.ipynb"
    }
    String prefix = basename(bam, ".bam")

    command <<<
    set -eux

        bash classify_reads.sh \
            ~{bam} \
            local_out_dir

        tree local_out_dir
    >>>

    Int disk_size = 20
    runtime {
        cpu:            2
        memory:         "8 GiB"
        disks:          "local-disk ~{disk_size} SSD"
        preemptible:    2
        maxRetries:     1
        docker:         "us.gcr.io/broad-dsp-lrma/mito_ccs_read_preprocessor:~{classifier_docker_tag}"
    }
}

task ExtractMitochondriaReads {
    meta {
        description:
        "For subsetting a high-coverage BAM stored in GCS, without localizing (more resilient to auth. expiration)."
    }

    parameter_meta {
        bam: {
            localization_optional: true
        }

        is_samtools_failed: "if true, the streaming of BAM from the bucket didn't succeed, so consider the result BAM corrupted."
    }

    input {
        File bam
        File bai
    }

    output {
        Boolean is_samtools_failed = read_boolean("samtools.failed.txt")
        File subset_bam = "~{subset_prefix}.bam"
        File subset_bai = "~{subset_prefix}.bam.bai"
    }

    Int disk_size = 20 # 20GiB reads mapping onto chrM is problematic

    String subset_prefix = basename(bam, ".bam") + ".chrM"

    command <<<

        # the way this works is the following:
        # 0) relying on the re-auth.sh script to export the credentials
        # 1) perform the remote sam-view subsetting in the background
        # 2) listen to the PID of the background process, while re-auth every 1200 seconds
        source /opt/re-auth.sh
        set -euxo pipefail

        echo "false" > "samtools.failed.txt"

        # see man page for what '-M' means
        samtools view \
            -bhX \
            -M \
            -@ 1 \
            --verbosity=8 \
            --write-index \
            -o "~{subset_prefix}.bam##idx##~{subset_prefix}.bam.bai" \
            ~{bam} ~{bai} \
            'chrM' && exit 0 || { echo "samtools seem to have failed"; echo "true" > "samtools.failed.txt"; exit 77; } &
        pid=$!

        set +e
        count=0
        while true; do
            sleep 1200 && date && source /opt/re-auth.sh
            count=$(( count+1 ))
            if [[ ${count} -gt 6 ]]; then echo "true" > "samtools.failed.txt" && exit 0; fi  # way too many attempts, get out
            if ! pgrep -x -P $pid; then exit 0; fi
        done
    >>>


    runtime {
        cpu:            2
        memory:         "8 GiB"
        disks:          "local-disk ~{disk_size} SSD"
        preemptible:    2
        maxRetries:     1
        docker:         "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}

task StopWorkflow {

    meta {
        description: "Utility to stop a workflow"
    }

    parameter_meta {
        reason: "reason for stopping"
    }

    input {
        String reason
    }
    command <<<
        echo -e "Workflow explicitly stopped because \n  ~{reason}." && exit 1
    >>>
    runtime {docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"}
}

task CountReads {
    meta {
        desciption: "Count the number of reads"
    }

    input {
        File chrM_bam
        Array[File] classified_bams
    }

    output {
        Map[String, Int] count_by_case = read_map("count_by_case.tsv")
    }

    command <<<
    set -eux

        # not using 'samtools view -c' because we want count of unique molecules, not alignment records
        cnt=$(samtools view ~{chrM_bam} | awk -F '\t' '{print $1}' | sort | uniq | wc -l | awk '{print $1}')
        reads_class='chrM'
        echo -e "${reads_class}\t${cnt}" >> count_by_case.tsv
        rm ~{chrM_bam}

        for bam in ~{sep=' ' classified_bams};
        do
            reads_class=$(echo "${bam}" | awk -F '.' '{print $(NF-1)}')
            cnt=$(samtools view "${bam}" | awk -F '\t' '{print $1}' | sort | uniq | wc -l | awk '{print $1}')
            echo -e "${reads_class}\t${cnt}" >> count_by_case.tsv
        done
        cat count_by_case.tsv

        all_reads_classes=("for_assembly" "multimap" "too_short" "ovp_pair_alns_lr" "break_cocircularity")
        for reads_class in "${all_reads_classes[@]}";
        do
            if ! grep -qF 'for_assembly' count_by_case.tsv; then echo -e "\t0" >> count_by_case.tsv; fi
        done

        sed -i.bak 's/ovp_pair_alns_lr/signif_overlap/' count_by_case.tsv
    >>>

    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk 10 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
}