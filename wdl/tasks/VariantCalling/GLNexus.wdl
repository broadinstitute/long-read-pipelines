version 1.0

import "../Utility/Utils.wdl"
import "../Utility/VariantUtils.wdl" as VarUtils

workflow JointCall {

    meta {
        description: "This pipeline joint-calls GVCFs with GLNexus (https://github.com/dnanexus-rnd/GLnexus). It also permits intervals to be specified so that joint calling only takes place on a subset of intervals (this can be useful for finding duplicate samples)."
    }

    parameter_meta {
        gvcfs:    "gVCF files to perform joint calling upon"
        tbis:     "gVCF index files"

        background_sample_gvcfs: "Array of GVCFs to use as background samples for joint calling."
        background_sample_gvcf_indices: "Array of GVCF index files for `background_sample_gvcfs`.  Order should correspond to that in `background_sample_gvcfs`."

        dict:     "reference sequence dictionary"
        bed:      "intervals to which joint calling should be restricted"

        config:   "configuration preset name or .yml filename"
        more_PL:  "include PL from reference bands and other cases omitted by default"
        squeeze:  "reduce pVCF size by suppressing detail in cells derived from reference bands"
        trim_uncalled_alleles: "remove alleles with no output GT calls in postprocessing"

        force_add_missing_dp: "force adding missing DP field to gVCFs"

        num_cpus: "number of CPUs to use"
        max_cpus: "maximum number of CPUs to allow"
        prefix:   "output prefix for joined-called BCF and GVCF files"

        runtime_attr_override: "override default runtime attributes"
    }

    input {
        Array[File] gvcfs
        Array[File] tbis

        Array[Array[File]]? background_sample_gvcfs
        Array[Array[File]]? background_sample_gvcf_indices

        File dict
        File? bed

        String config = "DeepVariantWGS"
        File? config_file
        
        Boolean more_PL = false
        Boolean squeeze = false
        Boolean trim_uncalled_alleles = false

        Boolean force_add_missing_dp = false

        Int? num_cpus
        Int max_cpus = 64
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int cpus_exp = if defined(num_cpus) then select_first([num_cpus]) else 2*length(gvcfs)
    Int cpus_act = if cpus_exp < max_cpus then cpus_exp else max_cpus

    # List all of the contigs in the reference
    call GetRanges { input: dict = dict, bed = bed }

    # Check for background samples and combine them with the input gVCFs if they exist:
    if (defined(background_sample_gvcfs)) {
        Array[File] flattened_background_sample_gvcfs = flatten(select_first([background_sample_gvcfs]))
        Array[File] flattened_background_sample_gvcf_indices = flatten(select_first([background_sample_gvcf_indices]))

        Array[File] gvcfs_and_background_samples = flatten([gvcfs, flattened_background_sample_gvcfs])
        Array[File] tbis_and_background_samples = flatten([tbis, flattened_background_sample_gvcf_indices])
    }
    Array[File] final_gvcfs = select_first([gvcfs_and_background_samples, gvcfs])
    Array[File] final_tbis = select_first([tbis_and_background_samples, tbis])

    # Shard all gVCFs into per-contig shards
    scatter (p in zip(final_gvcfs, final_tbis)) {
        call ShardVCFByRanges { input: gvcf = p.left, tbi = p.right, ranges = GetRanges.ranges }
    }

    # Joint-call in parallel over chromosomes
    scatter (i in range(length(ShardVCFByRanges.sharded_gvcfs[0]))) {
        Array[File] per_contig_gvcfs = transpose(ShardVCFByRanges.sharded_gvcfs)[i]

        call Call {
            input:
                gvcfs = per_contig_gvcfs,

                config = config,
                config_file = config_file,

                more_PL = more_PL,
                squeeze = squeeze,
                trim_uncalled_alleles = trim_uncalled_alleles,

                force_add_missing_dp = force_add_missing_dp,

                num_cpus = cpus_act,
                prefix = prefix
        }
    }

    # Concatenate the contig-sharded joint calls into a single joint callset
    call VarUtils.ConcatVariants as ConcatBCFs { input: variant_files = Call.joint_bcf, prefix = prefix }

    output {
        File joint_vcf = ConcatBCFs.combined_vcf
        File joint_vcf_tbi = ConcatBCFs.combined_vcf_tbi
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
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
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
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
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
        File? config_file

        Boolean more_PL = false
        Boolean squeeze = false
        Boolean trim_uncalled_alleles = false

        Boolean force_add_missing_dp = false

        Int num_cpus = 96
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 5*ceil(size(gvcfs, "GB")) * if (force_add_missing_dp) then 2 else 1
    Int mem = 4*num_cpus

    command <<<
        set -x

        # For guidance on performance settings, see https://github.com/dnanexus-rnd/GLnexus/wiki/Performance
        ulimit -Sn 65536

        echo ~{gvcfs[0]} | sed 's/.*locus_//' | sed 's/.g.vcf.bgz//' | sed 's/___/\t/g' > range.bed

        gvcf_file_list=~{write_lines(gvcfs)}

        if [[ "~{force_add_missing_dp}" == "true" ]]; then
            echo "Force adding missing DP field to gVCFs"
            fixed_gvcf_file_list="fixed_gvcf_file_list.txt"
            while read gvcf ; do
                echo "Processing ${gvcf}"
                bn=$(basename ${gvcf} | sed -e 's/\.gz$//' -e 's/\.bgz$//' -e 's/\.vcf$//' -e 's@\.g$@@')
                nn=${bn}.missing_dp_added.g.vcf.gz
                SM=$(bcftools query -l ${gvcf})
                bcftools view ${gvcf} | grep -v ':DP:' | grep -v '^#' | awk -F$'\t' 'BEGIN{OFS="\t"}{print $1,$2,"0"}' | bgzip -c > annot.txt.gz
                tabix -s1 -b2 -e2 annot.txt.gz
                bcftools annotate -Oz2 -o ${nn} -s "${SM}" -a annot.txt.gz -c CHROM,POS,FORMAT/DP ${gvcf}
                bcftools index -t ${nn}
                echo ${nn} >> ${fixed_gvcf_file_list}
            done < ${gvcf_file_list}

            gvcf_file_list=${fixed_gvcf_file_list}
        fi


        glnexus_cli \
            --config ~{if (defined(config_file)) then "~{config_file}" else "~{config}"} \
            --bed range.bed \
            ~{if more_PL then "--more-PL" else ""} \
            ~{if squeeze then "--squeeze" else ""} \
            ~{if trim_uncalled_alleles then "--trim-uncalled-alleles" else ""} \
            --list ${gvcf_file_list} \
            > ~{prefix}.bcf
    >>>

    output {
        File joint_bcf = "~{prefix}.bcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             mem,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
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

