version 1.0

import "../../structs/Structs.wdl"

task ConvertToZarrStore {
    meta {
        description: "Convert a .vcf.bgz file to a Zarr store and copy it to a final gs:// URL."
    }

    input {
        File vcf
        File tbi

        String reference = "GRCh38"
        String? ref_fasta
        String? ref_fai
        String prefix = "out"

        Int num_cpus = 4

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 6*ceil(size(vcf, "GB"))

    command <<<
        set -x

        # According to the sgkit documentation, running the conversion in parallel using Dask
        # is significantly faster than running it serially.
        # This can be easily done by importing the Dask library and setting the client info.
        # For more info see: https://sgkit-dev.github.io/sgkit/latest/vcf.html#example-converting-1000-genomes-vcf-to-zarr

        # Set up number of worker processes for Dask:
        num_processors=$(lscpu | grep '^CPU(s):' | awk '{print $NF}')
        num_workers=$((num_processors - 1))

        python3 <<EOF

        from sgkit.io.vcf import vcf_to_zarr

        vcfs = ["~{vcf}"]
        target = "~{prefix}.zarr"

        # Log our task status to a file we can inspect later:
        vcf_to_zarr(vcfs, target, tempdir="tmp")

        EOF

        echo "Done converting to zarr."

        echo "Tarring output..."
        tar -cf ~{prefix}.zarr.tar ~{prefix}.zarr
    >>>

    output {
        File zarr = "~{prefix}.zarr.tar"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
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
