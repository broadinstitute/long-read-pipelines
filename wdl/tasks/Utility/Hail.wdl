version 1.0

import "../../structs/Structs.wdl"

task ConvertToHailMT {
    meta {
        description: "Convert a .vcf.bgz file to a Hail MatrixTable."
    }

    parameter_meta {
        gvcf: "The input .vcf.bgz file."
        tbi: "The input .vcf.bgz.tbi file."
        reference: "The reference genome to use."
        ref_fasta: "The reference genome FASTA file.  If not specified, the reference genome will be downloaded from the Hail website."
        ref_fai: "The reference genome FASTA index file.  If not specified, the reference genome will be downloaded from the Hail website."
        prefix: "The prefix to use for the output MatrixTable (a folder)."
        outdir: "The output directory to copy the MatrixTable to, if provided, the MatrixTable will be copied to this location (with the prefix)."
        runtime_attr_override: "Override the default runtime attributes for this task."
    }

    input {
        File gvcf
        File tbi
        String prefix

        String reference = "GRCh38"
        String? ref_fasta
        String? ref_fai

        String? outdir

        RuntimeAttr? runtime_attr_override
    }

    output {
        File mt_tar = "~{prefix}.mt.tar"
        File completion_file = "completion_key_file"

        String gcs_path = if(defined(outdir)) then "~{outdir}/~{prefix}.mt" else "Null"
    }

    Boolean custom_ref = defined(ref_fasta) && defined(ref_fai)

    command <<<
        set -eux

        echo "Creating matrix table......"
        date
        python3 <<EOF

        import hail as hl
        hl.init(idempotent=True)

        if '~{custom_ref}' == 'true':
            if '~{reference}' in {'GRCh37', 'GRCh38', 'GRCm38', 'CanFam3'}:
                ref = hl.default_reference('~{reference}')
            else:
                ref = hl.ReferenceGenome.from_fasta_file('~{reference}', '~{ref_fasta}', '~{ref_fai}')
                hl.default_reference(ref)
        else:
            ref = hl.default_reference('GRCh38')

        callset = hl.import_vcf(
            '~{gvcf}',
            array_elements_required=False,
            force_bgz=True,
            reference_genome=ref
        )

        callset.write('~{prefix}.mt')

        EOF
        date

        ls

        echo "Finished creating matrix table, now tarring......"
        time \
        tar -cf ~{prefix}.mt.tar ~{prefix}.mt &

        ## Copy the matrix table to the output directory if provided
        if ~{defined(outdir)}; then
            echo "Copying the MT (folder) to the specified GS path......"
            gcloud storage cp -r \
                ~{prefix}.mt \
                ~{outdir}/~{prefix}.mt

            echo "Finished copying. Now verify the bucket copy isn't corrupt......"
            # the following code are left-shifted to satisfy the shell
        python3 <<EOF

        import hail as hl
        hl.init(idempotent=True)

        if '~{custom_ref}' == 'true':
            if '~{reference}' in {'GRCh37', 'GRCh38', 'GRCm38', 'CanFam3'}:
                ref = hl.default_reference('~{reference}')
            else:
                ref = hl.ReferenceGenome.from_fasta_file('~{reference}', '~{ref_fasta}', '~{ref_fai}')
                hl.default_reference(ref)
        else:
            ref = hl.default_reference('GRCh38')

        reloaded_mt = hl.read_matrix_table('~{outdir}/~{prefix}.mt')
        hl.summarize_variants(reloaded_mt)

        EOF
            # done with artificial left shift
            echo "If you see me, that means the bucket copy was intact!"
        fi

        wait
        touch completion_key_file
    >>>

    #########################
    Int disk_size = 1 + 6*ceil(size(gvcf, "GB"))

    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             96,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/hail:0.2.130-gcloud"
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
