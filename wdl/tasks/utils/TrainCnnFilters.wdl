version 1.0

workflow TrainCnnFilters {
    meta {
        author: "Jonn Smith"
        description: "A workflow for training the 1D and 2D CNN filtration methods in GATK."
    }

    input {

    }

    output {

    }
}


task Create1DReferenceTensors {

    meta {
        author: "Jonn Smith"
        description: "Task to create 1D reference tensors for the 1D CNN."
    }

    input {
        File vcf_input

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File truth_vcf
        File truth_bed

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(vcf_input, "GB") +
                               size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(ref_dict, "GB") +
                               size(truth_vcf, "GB") +
                               size(truth_bed, "GB")
                          )

    command <<<
        set -euxo pipefail

        export MONITOR_MOUNT_POINT="/cromwell_root"
        curl https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/jts_kvg_sp_malaria/scripts/monitor/legacy/vm_local_monitoring_script.sh > monitoring_script.sh
        chmod +x monitoring_script.sh
        ./monitoring_script.sh &> resources.log &
        monitoring_pid=$!

        gatk CNNVariantWriteTensors \
            -R ~{ref_fasta} \
            -V ~{vcf_input}  \
            -truth-vcf ~{truth_vcf} \
            -truth-bed ~{truth_bed} \
            -tensor-type reference \
            --downsample-snps 1 \
            --downsample-indels 1 \
            --max-tensors 10000000 \
            -output-tensor-dir ~{prefix}_1D_tensor_dir

        # No need to zip - the files are .hd5 formatted:
        tar -cf ~{prefix}_1D_tensor_dir.tar ~{prefix}_1D_tensor_dir

        kill $monitoring_pid
    >>>

    output {
        File monitoring_log = "resources.log"
        File tensor_dir_tar = "~{prefix}_1D_tensor_dir.tar"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.3.0.0"
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

task Create2DReadTensors {

    meta {
        author: "Jonn Smith"
        description: "Task to create 2D read tensors for the 2D CNN."
    }

    input {
        File bam_input
        File vcf_input

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File truth_vcf
        File truth_bed

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(bam_input, "GB") +
                               size(vcf_input, "GB") +
                               size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(ref_dict, "GB") +
                               size(truth_vcf, "GB") +
                               size(truth_bed, "GB")
                          )

    command <<<
        set -euxo pipefail

        export MONITOR_MOUNT_POINT="/cromwell_root"
        curl https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/jts_kvg_sp_malaria/scripts/monitor/legacy/vm_local_monitoring_script.sh > monitoring_script.sh
        chmod +x monitoring_script.sh
        ./monitoring_script.sh &> resources.log &
        monitoring_pid=$!

        gatk CNNVariantWriteTensors \
            -R ~{ref_fasta} \
            -V ~{vcf_input}  \
            -bam-file ~{bam_input} \
            -truth-vcf ~{truth_vcf} \
            -truth-bed ~{truth_bed} \
            -tensor-type read_tensor \
            --downsample-snps 1 \
            --downsample-indels 1 \
            --max-tensors 10000000 \
            -output-tensor-dir ~{prefix}_2D_tensor_dir

        # No need to zip - the files are .hd5 formatted:
        tar -cf ~{prefix}_2D_tensor_dir.tar ~{prefix}_2D_tensor_dir

        kill $monitoring_pid
    >>>

    output {
        File monitoring_log = "resources.log"
        File tensor_dir_tar = "~{prefix}_2D_tensor_dir.tar"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.3.0.0"
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

task TrainCnn1D {

    meta {
        author: "Jonn Smith"
        description: "Task to train the 1D CNN with 1D tensors"
    }

    input {
        Array[File] tensor_tars

        Int epochs = 100
        Int training_steps = 100
        Int validation_steps = 6

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 30*ceil(size(tensor_tars, "GB"))

    command <<<
        set -euxo pipefail

        export MONITOR_MOUNT_POINT="/cromwell_root"
        curl https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/jts_kvg_sp_malaria/scripts/monitor/legacy/vm_local_monitoring_script.sh > monitoring_script.sh
        chmod +x monitoring_script.sh
        ./monitoring_script.sh &> resources.log &
        monitoring_pid=$!

        # Must pre-process the given tensor_tars into a single folder:

        ${GATK} CNNVariantTrain \
            -tensor-type reference \
            --epochs ~{epochs} \
            --training-steps ~{training_steps} \
            --validation-steps ~{validation_steps} \
            -input-tensor-dir p_falciparum_ref_model_tensors_all_variants/ \
            -model-name ~{prefix}_CNN_1D_model \

        kill $monitoring_pid
    >>>

    output {
        File monitoring_log = "resources.log"
        File model_hd5 = "~{prefix}_CNN_1D_model.hd5"
        File model_json = "~{prefix}_CNN_1D_model.json"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.3.0.0"
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

task TrainCnn {

    meta {
        author: "Jonn Smith"
        description: "Task to train the CNN."
    }

    input {
        Array[File] tensor_tars
        String tensor_type  # Can be either "reference" or "read_tensor"

        Int epochs = 100
        Int training_steps = 100
        Int validation_steps = 6

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 30*ceil(size(tensor_tars, "GB"))

    command <<<
        set -euxo pipefail

        export MONITOR_MOUNT_POINT="/cromwell_root"
        curl https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/jts_kvg_sp_malaria/scripts/monitor/legacy/vm_local_monitoring_script.sh > monitoring_script.sh
        chmod +x monitoring_script.sh
        ./monitoring_script.sh &> resources.log &
        monitoring_pid=$!

        # Must pre-process the given tensor_tars into a single folder:
        mkdir tensors
        cd tensors
        while read f ; do
            tar --strip-components 1 -xf $f
        done < ${write_lines(tensor_tars)}
        cd ../

        ${GATK} CNNVariantTrain \
            -tensor-type reference \
            --epochs ~{epochs} \
            --training-steps ~{training_steps} \
            --validation-steps ~{validation_steps} \
            -input-tensor-dir tensors/ \
            -model-name ~{prefix}_CNN_~{tensor_type}_model \

        kill $monitoring_pid
    >>>

    output {
        File monitoring_log = "resources.log"
        File model_hd5 = "~{prefix}_CNN_~{tensor_type}_model.hd5"
        File model_json = "~{prefix}_CNN_~{tensor_type}_model.json"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.3.0.0"
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