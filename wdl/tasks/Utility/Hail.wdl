version 1.0

import "../../structs/Structs.wdl"

task ConvertToHailMT {
    meta {
        description: "Convert a .vcf.bgz file to a Hail MatrixTable and copy it to a final gs:// URL."
    }

    parameter_meta {
        gvcf: "The input .vcf.bgz file."
        tbi: "The input .vcf.bgz.tbi file."
        reference: "The reference genome to use.  Currently only GRCh38 is supported."
        ref_fasta: "The reference genome FASTA file.  If not specified, the reference genome will be downloaded from the Hail website."
        ref_fai: "The reference genome FASTA index file.  If not specified, the reference genome will be downloaded from the Hail website."
        prefix: "The prefix to use for the output MatrixTable."
        runtime_attr_override: "Override the default runtime attributes for this task."
    }

    input {
        File gvcf
        File tbi

        String reference = "GRCh38"
        String? ref_fasta
        String? ref_fai
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 6*ceil(size(gvcf, "GB"))

    command <<<
        set -x

        python3 <<EOF

        import hail as hl
        hl.init(default_reference='GRCh38')

        if '~{defined(ref_fasta)}' == 'true' and '~{defined(ref_fai)}' == 'true':
            ref = hl.ReferenceGenome.from_fasta_file('~{reference}', '~{ref_fasta}', '~{ref_fai}')

        callset = hl.import_vcf(
            '~{gvcf}',
            array_elements_required=False,
            force_bgz=True,
            reference_genome='~{reference}'
        )

        callset.write('~{prefix}.mt')

        EOF

        echo "Created matrix table."
        echo "Tarring file now."
        tar -cf ~{prefix}.mt.tar ~{prefix}.mt

        touch completion_key_file
    >>>

    output {
        File mt_tar = "~{prefix}.mt.tar"
        File completion_file = "completion_key_file"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
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
