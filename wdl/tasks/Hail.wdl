version 1.0

import "Structs.wdl"

task ConvertToHailMT {
    input {
        File gvcf
        File tbi

        String reference = "GRCh38"
        String? ref_fasta
        String? ref_fai
        String prefix = "out"

        String outdir

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 5*ceil(size(gvcf, "GB"))

    command <<<
        set -x

        python3 <<EOF

        import hail as hl

        print('~{defined(ref_fasta)}')

        if '~{defined(ref_fasta)}' == 'true':
            ref = hl.ReferenceGenome.from_fasta_file('~{reference}', '~{ref_fasta}', '~{ref_fai}')
            print("Hello!")
        else:
            hl.init(default_reference='~{reference}')

        callset = hl.import_vcf('~{gvcf}', array_elements_required=False)
        callset.write('~{prefix}.mt', overwrite=True)

        EOF

        gsutil -m rsync -Cr ~{prefix}.mt ~{outdir}/~{prefix}.mt
    >>>

    output {
        String gcs_path = "~{outdir}/~{prefix}.mt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             64,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "hailgenetics/hail:0.2.105"
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
