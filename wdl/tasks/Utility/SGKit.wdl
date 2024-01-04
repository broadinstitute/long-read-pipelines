version 1.0

import "../../structs/Structs.wdl"

task ConvertToZarrStore {
    meta {
        description: "Convert a .vcf.bgz file to a Zarr store and copy it to a final gs:// URL."
    }

    input {
        File gvcf
        File tbi

        String reference = "GRCh38"
        String? ref_fasta
        String? ref_fai
        String prefix = "out"

        Int num_cpus = 4

        String outdir

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 3*ceil(size(gvcf, "GB"))

    command <<<
        set -x

        python3 <<EOF

        from sgkit.io.vcf import vcf_to_zarr

        vcfs = ["~{gvcf}"]
        target = "~{prefix}.zarr"

        vcf_to_zarr(vcfs, target, tempdir="tmp")

        EOF

        gsutil -m rsync -Cr ~{prefix}.zarr ~{outdir}/~{prefix}.zarr
    >>>

    output {
        String gcs_path = "~{outdir}/~{prefix}.zarr"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-sgkit:0.5.0"
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
