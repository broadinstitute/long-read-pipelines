version 1.0

import "../../../tasks/Utility/Utils.wdl"

workflow JointCall {

    meta {
        description: "This pipeline joint-calls GVCFs with GLNexus (https://github.com/dnanexus-rnd/GLnexus). It also permits intervals to be specified so that joint calling only takes place on a subset of intervals (this can be useful for finding duplicate samples)."
    }

    parameter_meta {
        gvcfs:    "gVCF files to perform joint calling upon"
        tbis:     "gVCF index files"
        gvcf_and_tbi_manifest: "manifest file containing gVCF and TBI pairs. Provide either this, or gVCFs and TBIs, but not both."
        dict:     "reference sequence dictionary"
        bed:      "intervals to which joint calling should be restricted"

        config:   "configuration preset name or .yml filename"
        more_PL:  "include PL from reference bands and other cases omitted by default"
        squeeze:  "reduce pVCF size by suppressing detail in cells derived from reference bands"
        trim_uncalled_alleles: "remove alleles with no output GT calls in postprocessing"

        num_cpus: "number of CPUs to use"
        max_cpus: "maximum number of CPUs to allow"
        prefix:   "output prefix for joined-called BCF and GVCF files"

        runtime_attr_override: "override default runtime attributes"
    }

    input {
        Array[File]? gvcfs
        Array[File]? tbis

        File? gvcf_and_tbi_manifest

        File dict
        File? bed

        String config = "DeepVariantWGS"
        Boolean more_PL = false
        Boolean squeeze = false
        Boolean trim_uncalled_alleles = false

        Int? num_cpus
        Int max_cpus = 64
        String prefix = "out"
    }

    if (defined(gvcfs) != defined(tbis)) {
        call Utils.StopWorkflow as PairedInputsMissing { input: reason = "gVCF array and tbi must be specified at the same time." }
    }
    if ( defined(gvcf_and_tbi_manifest) == defined(gvcfs) ) {
        call Utils.StopWorkflow as MutexInputsGiven { input: reason = "Provide either gVCF array and TBI array or a manifest file, not both." }
    }

    if (defined(gvcf_and_tbi_manifest)) {
        call SplitManifest { input: manifest = select_first([gvcf_and_tbi_manifest]) }
    }
    Array[Pair[String, String]] gvcfs_and_tbis = if defined(gvcf_and_tbi_manifest)
                                                 then select_first([SplitManifest.g_N_t])
                                                 else zip( select_first([gvcfs]), select_first([tbis]) )


    Int cpus_exp = if defined(num_cpus) then select_first([num_cpus]) else 2*length(gvcfs_and_tbis)
    Int cpus_act = if cpus_exp < max_cpus then cpus_exp else max_cpus

    # List all of the contigs in the reference
    call GetRanges { input: dict = dict, bed = bed }

    # Shard all gVCFs into per-contig shards
    scatter (p in gvcfs_and_tbis) {
        call ShardVCFByRanges { input: gvcf = p.left, tbi = p.right, ranges = GetRanges.ranges }
    }

    Array[Array[File]] per_range_per_sample_gvcfs = transpose(ShardVCFByRanges.sharded_gvcfs)
    # Joint-call in parallel over chromosomes
    scatter (i in range(length(GetRanges.ranges))) {
        Array[File] per_contig_gvcfs = per_range_per_sample_gvcfs[i]

        call Call {
            input:
                gvcfs = per_contig_gvcfs,

                config = config,
                more_PL = more_PL,
                squeeze = squeeze,
                trim_uncalled_alleles = trim_uncalled_alleles,

                num_cpus = cpus_act,
                prefix = prefix
        }
    }

    # Concatenate the contig-sharded joint calls into a single joint callset
    call ConcatBCFs { input: bcfs = Call.joint_bcf, prefix = prefix }

    output {
        File joint_gvcf = ConcatBCFs.joint_gvcf
        File joint_gvcf_tbi = ConcatBCFs.joint_gvcf_tbi
    }
}

task SplitManifest {
    meta {
        description: "Split a 2-col manifest into two arrays of URIs, each holding [gvcf] and [tbi]"
    }
    parameter_meta {
        manifest: "2-col manifest file, holding URI to [(gvcf, tbi)]"
    }
    input {
        File manifest
    }
    output {
        Array[Pair[String, String]] g_N_t = zip(read_lines("gvcfs.txt"), read_lines("tbis.txt"))
    }
    command <<<
    set -euxo pipefail

        cut -f1 ~{manifest} > gvcfs.txt
        cut -f2 ~{manifest} > tbis.txt
    >>>
    runtime {
        disks: "local-disk 10 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}

task GetRanges {
    meta {
        description: "Select loci over which to parallelize downstream operations."
    }

    input {
        File dict
        File? bed

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + ceil(size(dict, "GB"))

    command <<<
        set -euxo pipefail

        if [[ "~{defined(bed)}" == "true" ]]; then
            cat ~{bed} | awk '{ print $1 ":" $2 "-" $3 }' > ranges.txt
        else
            grep '^@SQ' ~{dict} | \
                awk '{ print $2, $3 }' | \
                sed 's/[SL]N://g' | \
                awk '{ print $1 ":0-" $2 }' \
                > ranges.txt
        fi
    >>>

    output {
        Array[String] ranges = read_lines("ranges.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "ghcr.io/dnanexus-rnd/glnexus:v1.4.3"
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

task ShardVCFByRanges {
    meta {
        description: "Split VCF into smaller ranges for parallelization."
    }

    input {
        File gvcf
        File tbi
        Array[String] ranges

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(gvcf, "GB"))

    command <<<
        set -euxo pipefail

        mkdir per_contig

        INDEX=0
        for RANGE in ~{sep=' ' ranges}
        do
            PINDEX=$(printf "%06d" $INDEX)
            FRANGE=$(echo $RANGE | sed 's/[:-]/___/g')
            OUTFILE="per_contig/$PINDEX.~{basename(gvcf, ".g.vcf.gz")}.locus_$FRANGE.g.vcf.gz"

            bcftools view ~{gvcf} $RANGE | bgzip > $OUTFILE

            INDEX=$(($INDEX+1))
        done
    >>>

    output {
        Array[File] sharded_gvcfs = glob("per_contig/*")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "ghcr.io/dnanexus-rnd/glnexus:v1.4.3"
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

task Call {
    meta {
        description: "Joint-call gVCFs with GLNexus."
    }

    parameter_meta {
        gvcfs: {
            localization_optional: true
        }
    }

    input {
        Array[File] gvcfs

        String config = "DeepVariantWGS"
        Boolean more_PL = false
        Boolean squeeze = false
        Boolean trim_uncalled_alleles = false

        Int num_cpus = 96
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    # this will break call-caching
    Array[String] a = gvcfs
    File gs_manifest = write_lines(a)
    File local_manifest = write_lines(gvcfs)

    command <<<
        set -x

        head ~{gs_manifest}
        head ~{local_manifest}

        time \
        cat ~{gs_manifest} \
        | gcloud storage cp \
            -I \
            .

        # For guidance on performance settings, see https://github.com/dnanexus-rnd/GLnexus/wiki/Performance
        ulimit -Sn 65536

        echo ~{gvcfs[0]} | sed 's/.*locus_//' | sed 's/.g.vcf.bgz//' | sed 's/___/\t/g' > range.bed

        glnexus_cli \
            --config ~{config} \
            --bed range.bed \
            ~{if more_PL then "--more-PL" else ""} \
            ~{if squeeze then "--squeeze" else ""} \
            ~{if trim_uncalled_alleles then "--trim-uncalled-alleles" else ""} \
            --list ~{local_manifest} \
            > ~{prefix}.bcf
    >>>

    output {
        File joint_bcf = "~{prefix}.bcf"
    }

    #########################
    Int disk_size = 1 + 5*ceil(size(gvcfs, "GiB"))
    Int mem = 4*num_cpus
    String prefix_padded_disk_type = if (10000<length(gvcfs)) then " LOCAL" else " SSD"

    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             mem,
        disk_gb:            disk_size,
        boot_disk_gb:       20,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "ghcr.io/dnanexus-rnd/glnexus:v1.4.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + prefix_padded_disk_type
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CompressAndIndex {
    meta {
        description: "Convert a BCF file to a vcf.bgz file and index it."
    }

    input {
        File joint_bcf

        Int num_cpus = 8
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 3*ceil(size(joint_bcf, "GB"))

    command <<<
        set -x

        bcftools view ~{joint_bcf} | bgzip -@ ~{num_cpus} -c > ~{prefix}.g.vcf.bgz
        tabix -p vcf ~{prefix}.g.vcf.bgz
    >>>

    output {
        File joint_gvcf = "~{prefix}.g.vcf.bgz"
        File joint_gvcf_tbi = "~{prefix}.g.vcf.bgz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             4*num_cpus,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "ghcr.io/dnanexus-rnd/glnexus:v1.4.3"
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

task ConcatBCFs {
    meta {
        description: "Concatenate BCFs into a single .vcf.bgz file and index it."
    }

    input {
        Array[File] bcfs

        Int num_cpus = 10
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(bcfs, "GB"))

    command <<<
        set -euxo pipefail

        bcftools concat --threads 4 -n ~{sep=' ' bcfs} | bcftools view | bgzip -@ 4 -c > ~{prefix}.g.vcf.bgz
        tabix -p vcf ~{prefix}.g.vcf.bgz
    >>>

    output {
        File joint_gvcf = "~{prefix}.g.vcf.bgz"
        File joint_gvcf_tbi = "~{prefix}.g.vcf.bgz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             60,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "ghcr.io/dnanexus-rnd/glnexus:v1.4.3"
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
