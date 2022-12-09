version 1.0

##########################################################################################
# This pipeline joint-calls GVCFs with GLNexus (https://github.com/dnanexus-rnd/GLnexus).
# It also permits intervals to be specified so that joint calling only takes place on a
# subset of intervals (this can be useful for finding duplicate samples).
##########################################################################################

import "Utils.wdl"
import "VariantUtils.wdl"

workflow JointCall {
    input {
        Array[File] gvcfs
        File? bed

        String config = "DeepVariantWGS"
        Boolean more_PL = false
        Boolean squeeze = false
        Boolean trim_uncalled_alleles = false

        Int num_cpus = 96
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        gvcfs:    "gVCF files to perform joint calling upon"
        bed:      "three-column BED file with ranges to analyze (if not specified, use full length of all contigs)"

        config:   "configuration preset name or .yml filename"
        more_PL:  "include PL from reference bands and other cases omitted by default"
        squeeze:  "reduce pVCF size by suppressing detail in cells derived from reference bands"
        trim_uncalled_alleles: "remove alleles with no output GT calls in postprocessing"

        num_cpus: "number of CPUs to use"
        prefix:   "output prefix for joined-called BCF and GVCF files"
    }

    call Utils.ComputeAllowedLocalSSD as Guess { input: intended_gb = 2*ceil(size(gvcfs, "GB")) }

    call Call {
        input:
            gvcfs = gvcfs,
            bed = bed,

            config = config,
            more_PL = more_PL,
            squeeze = squeeze,
            trim_uncalled_alleles = trim_uncalled_alleles,

            num_cpus = num_cpus,
            prefix = prefix,

            num_ssds = Guess.numb_of_local_ssd
    }

    output {
        File joint_gvcf = Call.joint_gvcf
        File joint_gvcf_tbi = Call.joint_gvcf_tbi
    }
}

task Call {
    input {
        Array[File] gvcfs
        File? bed

        String config = "DeepVariantWGS"
        Boolean more_PL = false
        Boolean squeeze = false
        Boolean trim_uncalled_alleles = false

        Int num_ssds = 1
        Int num_cpus = 32
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 5*ceil(size(gvcfs, "GB"))

    command <<<
        set -x

        df -h

        # For guidance on performance settings, see https://github.com/dnanexus-rnd/GLnexus/wiki/Performance
        ulimit -Sn 65536

        glnexus_cli \
            --dir /mnt/tmp/GLnexus.DB \
            --config ~{config} \
            ~{if more_PL then "--more-PL" else ""} \
            ~{if squeeze then "--squeeze" else ""} \
            ~{if trim_uncalled_alleles then "--trim-uncalled-alleles" else ""} \
            ~{if defined(bed) then "--bed ~{select_first([bed])}" else ""} \
            --list ~{write_lines(gvcfs)} \
            | bcftools view | bgzip -@ ~{num_cpus} -c > ~{prefix}.g.vcf.bgz

        tabix -p vcf ~{prefix}.g.vcf.bgz

        df -h
    >>>

    output {
        File joint_gvcf = "~{prefix}.g.vcf.bgz"
        File joint_gvcf_tbi = "~{prefix}.g.vcf.bgz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             4*num_cpus,
        disk_gb:            0, # ignored
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "ghcr.io/dnanexus-rnd/glnexus:v1.4.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
#        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        disks: "local-disk ~{375*num_ssds} LOCAL, /mnt/tmp ~{disk_size} SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
