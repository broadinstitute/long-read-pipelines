version 1.0

import "tasks/Utils.wdl" as Utils

workflow dorado_basecall {

    input {
        String sample_id
        String gcs_fast5_dir
        String basecall_model
        String gcs_out_root_dir
    }

    call Utils.ListFilesOfType { input: gcs_dir = gcs_fast5_dir, suffixes = [ ".fast5" ] }
    call Utils.ChunkManifest { input: manifest = ListFilesOfType.manifest, manifest_lines_per_chunk = 30 }

    scatter (manifest_chunk in ChunkManifest.manifest_chunks) {
        call basecall {
            input:
                sample_id = sample_id,
                fast5_files = read_lines(manifest_chunk),
                basecall_model = basecall_model
        }
    }

    call Utils.MergeBams { input: bams = basecall.unmapped_bam }

    output {
        File unmapped_bam = MergeBams.merged_bam
        File unmapped_bai = MergeBams.merged_bai
    }

    meta {
        author: "Martin Aryee"
        email:"martin.aryee@gmail.com"
    }
}
task basecall  {
    input {
        String sample_id
        Array[File] fast5_files
        String basecall_model
    }

    Int disk_gb = ceil(size(fast5_files, "GB")*3) + 5

    command <<<
        set -euxo pipefail

        fast5_dir=$(dirname ~{fast5_files[0]})
        mkdir pod5

        # convert to pod5
        pod5 convert fast5 ${dir}/*fast5 --output pod5s/ --one-to-one --threads 12

        # make the calls
        dorado basecaller /models/~{basecall_model} pod5s/ --modified-bases 5mCG > ~{sample_id}.bam
    >>>
    runtime {
        gpuType: "nvidia-tesla-v100"
        gpuCount: 1
        cpu: 12
        disks: "local-disk " + disk_gb + " SSD" 
        memory: "32GB"
        nvidiaDriverVersion: "470.161.03"
        zones: ["us-central1-a"] 
        docker: "nanoporetech/dorado:latest"
    }
    output {
        File unmapped_bam = "~{sample_id}.bam"
    }
}
