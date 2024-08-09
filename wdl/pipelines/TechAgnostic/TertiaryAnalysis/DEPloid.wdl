version 1.0

import "../../../structs/Structs.wdl"

workflow DEPloid {

    meta {
        description: "Run DEPloid on a VCF file to determine complexity of infection and perform variant deconvolution."
    }
    parameter_meta {
        sample_name: "Sample name"
        vcf: "VCF file to process"
        allele_freq_tsv: "TSV file with allele frequencies for each alt allele at each locus."
        ref_map_file:  "Reference map file indicating reference sequence and auxillary file locations"
    }

    input {
        String sample_name
        File vcf
        File ref_map_file
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    # 1 - normalize VCF file
    call NormalizeVcfForAlleleSubsetting as t01_NormalizeVcfForAlleleSubsetting {
        input:
            prefix = sample_name,
            vcf = vcf,
            ref_fasta = ref_map['fasta']
    }

    # 2 - create subsetted PLAF and VCF

    # 3 - run DEPloid

    output {
        String hrp2 = HRP2Status.status
        String hrp3 = HRP3Status.status
    }
}

task NormalizeVcfForAlleleSubsetting {
    input {
        String prefix
        File vcf
        File ref_fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 15*ceil(size(vcf, "GB"))

    command <<<
        input_file=~{vcf}

        echo "${input_file}" | grep -q '.gz$'
        r=$?
        if [ $r -eq 0 ] ; then
          gunzip ${input_file}
          input_file=$(echo "${1}" | sed 's@.gz$@@')
        fi

        set -euxo pipefail

        REF=~{ref_fasta}
        bn=~{prefix}

        # First fudge VQSLOD because deploid needs it.
        # We can swap in the `SCORE` field from VQSR.
        # This is not great, but it will work...
        chrom_line_num=$(grep -n '^#CHROM' ${input_file} | sed 's@:.*@@')
        sed -e "${chrom_line_num}i##INFO=<ID=VQSLOD,Number=1,Type=Float,Description=\"FUDGED FIELD.  THIS IS THE SAME AS SCORE.  ORIGINAL DESCRIPTION: Log odds of being a true variant versus being false under the trained gaussian mixture model\">" \
          -e 's@;SCORE=\([-0-9\.]*\);@;SCORE=\1;VQSLOD=\1;@' ${input_file} > ${bn}.VQSLOD_FUDGE.vcf

        # Now remove potentially invalid fields and normalize the vcf.
        # Normalized to one alt allele per line
        bcftools annotate -x"INFO/HAPCOMP,INFO/HAPDOM,INFO/HEC" ${bn}.VQSLOD_FUDGE.vcf | bcftools view -f "PASS,." | bcftools norm -m -any --atom-overlaps . -f ${REF} -Ov -o ${bn}.VQSLOD_FUDGE.norm.vcf

        # Cleanup the intermediate file:
        rm ${bn}.VQSLOD_FUDGE.vcf
    >>>

    output {
        File normalized_vcf = "${prefix}.VQSLOD_FUDGE.norm.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
