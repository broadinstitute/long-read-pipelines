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
        chromosomes: "Reference contigs to shard over, in output order"

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

        Array[String] chromosomes = [
            "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
            "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
            "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
            "chr22", "chrX", "chrY", "chrM"
        ]
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

    call GetRanges { input: dict = dict, chromosomes = chromosomes }

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

task GetRanges {
    meta {
        description: "Select loci over which to parallelize downstream operations."
    }

    input {
        File dict
        Array[String] chromosomes

        RuntimeAttr? runtime_attr_override
    }

    File chromosome_list = write_lines(chromosomes)
    Int disk_size = 1 + ceil(size(dict, "GB"))

    command <<<
        set -euo pipefail

        awk -v chroms="~{chromosome_list}" '
            BEGIN { FS="\t"; }
            /^@SQ/ {
                name=""; len="";
                for (i=1; i<=NF; i++) {
                    if ($i ~ /^SN:/) name=substr($i,4);
                    else if ($i ~ /^LN:/) len=substr($i,4);
                }
                if (name!="") lengths[name]=len;
            }
            END {
                while ((getline chr < chroms) > 0) {
                    if (!(chr in lengths)) {
                        printf("ERROR: chromosome %s was not found in the sequence dictionary\n", chr) > "/dev/stderr";
                        exit 1;
                    }
                    print chr ":0-" lengths[chr];
                }
            }
        ' ~{dict} > ranges.txt
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
        set -euo pipefail

        mkdir staged_gvcf
        mkdir per_contig
        rm -f sharded_gvcfs.txt
        ln -sf ~{gvcf} staged_gvcf/input.g.vcf.gz
        INDEX_BASENAME=$(basename "~{tbi}")
        if [[ "$INDEX_BASENAME" == *.csi ]]; then
            ln -sf ~{tbi} staged_gvcf/input.g.vcf.gz.csi
        else
            ln -sf ~{tbi} staged_gvcf/input.g.vcf.gz.tbi
        fi

        INDEX=0
        for RANGE in ~{sep=' ' ranges}
        do
            PINDEX=$(printf "%06d" $INDEX)
            FRANGE=$(echo $RANGE | sed 's/[:-]/___/g')
            OUTFILE="per_contig/$PINDEX.~{basename(gvcf, ".g.vcf.gz")}.locus_$FRANGE.g.vcf.gz"

            echo "[$PINDEX] $RANGE -> $OUTFILE" 1>&2
            bcftools view staged_gvcf/input.g.vcf.gz $RANGE | bgzip > $OUTFILE
            echo "$OUTFILE" >> sharded_gvcfs.txt

            INDEX=$(($INDEX+1))
        done

        N_SHARDS=$(wc -l < sharded_gvcfs.txt)
        echo "ShardVCFByRanges wrote ${N_SHARDS} shards for ~{basename(gvcf)}" 1>&2
    >>>

    output {
        Array[File] sharded_gvcfs = read_lines("sharded_gvcfs.txt")
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

    command <<<
        set -euxo pipefail

        # For guidance on performance settings, see https://github.com/dnanexus-rnd/GLnexus/wiki/Performance
        ulimit -Sn 65536

        echo ~{gvcfs[0]} | sed 's/.*locus_//' | sed -E 's/\.g\.vcf\.(bgz|gz)$//' | sed 's/___/\t/g' > range.bed

        glnexus_cli \
            --config ~{config} \
            --bed range.bed \
            ~{if more_PL then "--more-PL" else ""} \
            ~{if squeeze then "--squeeze" else ""} \
            ~{if trim_uncalled_alleles then "--trim-uncalled-alleles" else ""} \
            --list ~{write_lines(gvcfs)} \
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
