version 1.0

import "../../../structs/Structs.wdl" as Structs
import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow TrainCnnFilters {
    meta {
        author: "Jonn Smith"
        description: "A workflow for training the 1D and 2D CNN filtration methods in GATK."
    }

    input {
        Array[File] vcfs
        Array[File] vcf_indices

        Array[File] bams
        Array[File] bais

        Array[File] truth_vcfs
        Array[File] truth_vcf_indices
        Array[File] truth_beds

        File ref_map_file

        String prefix = "out"
    }

    parameter_meta {
        vcfs:               "GCS path to VCF files containing called variants on which to train / test / validate the CNN models."
        vcf_indices:        "GCS path to index files for called variants on which to train / test / validate the CNN models."

        bams:                "GCS path to bam files containing the either the mapped reads from which variants were called, or a bam-out from the variant caller that produced the input VCF files."
        bais:                "GCS path to index files for the bam files containing the either the mapped reads from which variants were called, or a bam-out from the variant caller that produced the input VCF files."

        truth_vcfs:         "GCS path to VCF files containing validated variant calls (\"truth\")  for the corresponding called variants in `vcfs`."
        truth_vcf_indices:  "GCS path to index files for VCF files containing validated variant calls (\"truth\") for the corresponding called variants in `vcfs`."
        truth_beds:         "GCS path to bed files with confident regions for the given `truth_vcfs`"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"
    }

    # Get ref info:
    Map[String, String] ref_map = read_map(ref_map_file)

    # TODO: Validate that lengths of all inputs are the same:
    if ((length(vcfs) != length(vcf_indices)) || (length(vcfs) != length(vcf_indices)) || (length(vcfs) != length(vcf_indices)) || (length(vcfs) != length(vcf_indices)) || (length(vcfs) != length(vcf_indices))) {
        call Utils.StopWorkflow {input: reason="Not all input arrays have the same length!"}
    }

    # First create tensors for the input data:
    scatter (idx_1 in range(length(vcfs))) {
        # 1D CNN:
        call Create1DReferenceTensors {
            input:
                vcf_input = vcfs[idx_1],
                vcf_idx = vcf_indices[idx_1],
                ref_fasta = ref_map['fasta'],
                ref_fasta_fai = ref_map["fai"],
                ref_dict  = ref_map['dict'],
                truth_vcf = truth_vcfs[idx_1],
                truth_vcf_idx = truth_vcf_indices[idx_1],
                truth_bed = truth_beds[idx_1],
                prefix = prefix + "_shard_" + idx_1 + "reference"
        }
        # 2D CNN:
        call Create2DReadTensors {
            input:
                bam_input = bams[idx_1],
                bai_input = bais[idx_1],
                vcf_input = vcfs[idx_1],
                vcf_idx = vcf_indices[idx_1],
                ref_fasta = ref_map['fasta'],
                ref_fasta_fai = ref_map["fai"],
                ref_dict  = ref_map['dict'],
                truth_vcf = truth_vcfs[idx_1],
                truth_vcf_idx = truth_vcf_indices[idx_1],
                truth_bed = truth_beds[idx_1],
                prefix = prefix + "_shard_" + idx_1 + "reference"
        }
    }

    # Train the models with the created Tensors:
    # CNN 1D:
    call TrainCnn as TrainCnn1D {
        input:
            tensor_tars = Create1DReferenceTensors.tensor_dir_tar,
            tensor_type = "reference",
            epochs = 100,
            training_steps = 100,
            validation_steps = 6,
            prefix = prefix + "_CNN_1D_Model"
    }

    # CNN 2D:
    call TrainCnn as TrainCnn2D {
        input:
            tensor_tars = Create2DReadTensors.tensor_dir_tar,
            tensor_type = "read_tensor",
            epochs = 100,
            training_steps = 100,
            validation_steps = 6,
            optimizer_learning_rate = 0.000001,
            prefix = prefix + "_CNN_2D_Model"
    }

    output {
        Array[File] cnn_1d_tensors = Create1DReferenceTensors.tensor_dir_tar
        Array[File] cnn_2d_tensors = Create2DReadTensors.tensor_dir_tar

        File cnn_1d_model_json = TrainCnn1D.model_json
        File cnn_1d_model_hd5 = TrainCnn1D.model_hd5

        File cnn_2d_model_json = TrainCnn2D.model_json
        File cnn_2d_model_hd5 = TrainCnn2D.model_hd5
    }
}

task Create1DReferenceTensors {

    meta {
        author: "Jonn Smith"
        description: "Task to create 1D reference tensors for the 1D CNN."
    }

    input {
        File vcf_input
        File vcf_idx

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File truth_vcf
        File truth_vcf_idx
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
    >>>

    output {
        File tensor_dir_tar = "~{prefix}_1D_tensor_dir.tar"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
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
        File bai_input

        File vcf_input
        File vcf_idx

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File truth_vcf
        File truth_vcf_idx
        File truth_bed

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + 8*ceil(size(bam_input, "GB") +
                                size(vcf_input, "GB") +
                                size(ref_fasta, "GB") + size(ref_fasta_fai, "GB") + size(ref_dict, "GB") +
                                size(truth_vcf, "GB") +
                                size(truth_bed, "GB")
                          )

    command <<<
        set -euxo pipefail

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

# Now check if the tensors contain NaN values:
python << CODE

import os
import h5py

import numpy as np

print()

def find_files(folder_path, ext):
    all_files = []
    for root, directories, files in os.walk(folder_path):
        for file in files:
            if file.endswith(ext):
                all_files.append(os.path.join(root, file))
    return all_files

hd5_files = find_files("~{prefix}_2D_tensor_dir", "hd5")

print(f"Inspecting {len(hd5_files)} hd5 files... ")
for f in hd5_files:
    with h5py.File(f, 'r') as hd5:
        for k in hd5.keys():
            n = np.isnan(np.array(hd5[k]))
            if n.sum() > 0:
                print(f"File: {f}: Found {n.sum()} NaN(s) in key: {k}")
print("Done.")
print()

CODE

        # No need to zip - the files are .hd5 formatted:
        tar -cf ~{prefix}_2D_tensor_dir.tar ~{prefix}_2D_tensor_dir
    >>>

    output {
        File tensor_dir_tar = "~{prefix}_2D_tensor_dir.tar"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
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

        Float optimizer_beta1 = 0.9
        Float optimizer_beta2 = 0.999
        Float optimizer_clipnorm = 1.0
        Float optimizer_epsilon = 0.00000001       # 1.0e-8
        Float optimizer_learning_rate = 0.0001     # 1.0e-4

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    # We need a lot of disk space here for the unpacked tensors and final model:
    Int disk_size = 1 + 4*ceil(size(tensor_tars, "GB"))

    command <<<
        set -euxo pipefail

        # Must pre-process the given tensor_tars into a single folder:
        mkdir tensors
        cd tensors

        # Let's try to do this multi-threaded:
        # NOTE: Yes, I know this is multi-processing, but I'm in a hurry here.

        # Get the max number of threads to use:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        max_threads=$((np-1))
        if [[ $max_threads -le 0 ]] ; then
            max_threads=1
        fi

        # Dispatch some jobs:
        num_active_threads=0
        while read f ; do
            # If we have reached the maximum number of threads, we should wait for a while:
            if [[ $num_active_threads -ge $max_threads ]] ; then
                # Wait for the next background process to finish:
                wait -n

                # Give ourselves some wiggle room:
                sleep 5

                # Refresh the number of active threads:
                num_active_threads=$(jobs | wc -l)
            fi

            # Extract our tensors to the `tensors` folder:
            tar --strip-components 1 -xf $f &

            # Update the number of active threads:
            num_active_threads=$((num_active_threads+1))
        done < ~{write_lines(tensor_tars)}

        # Wait for the rest of our background processes to finish:
        wait

        cd ../

        gatk CNNVariantTrain \
            --verbosity DEBUG \
            -tensor-type reference \
            --epochs ~{epochs} \
            --training-steps ~{training_steps} \
            --validation-steps ~{validation_steps} \
            --optimizer-beta-1 ~{optimizer_beta1} \
            --optimizer-beta-2 ~{optimizer_beta2}  \
            --optimizer-clipnorm ~{optimizer_clipnorm}  \
            --optimizer-epsilon ~{optimizer_epsilon}  \
            --optimizer-learning-rate  ~{optimizer_learning_rate}  \
            -input-tensor-dir tensors/ \
            -model-name ~{prefix}_CNN_~{tensor_type}_model \

        ls -la
    >>>

    output {
        File model_hd5 = "~{prefix}_CNN_~{tensor_type}_model.hd5"
        File model_json = "~{prefix}_CNN_~{tensor_type}_model.json"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "broadinstitute/gatk-nightly:2023-08-18-4.4.0.0-57-g98f63667a-NIGHTLY-SNAPSHOT"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    # NOTE: We NEED GPUs to train the CNNs, so we don't allow for them to be modified by runtime attributes.
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        gpuType:                "nvidia-tesla-t4"
        gpuCount:               4
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
