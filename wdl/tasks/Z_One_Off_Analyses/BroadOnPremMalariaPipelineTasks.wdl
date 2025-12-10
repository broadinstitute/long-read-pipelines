version 1.0

import "../../structs/Structs.wdl"

task VariantRecalibrator {
    input {
        String prefix

        File input_vcf

        File reference_fasta
        File reference_fai
        File reference_dict

        File resource_vcf_7g8_gb4
        File resource_vcf_hb3_dd2
        File resource_vcf_3d7_hb3

        File resource_vcf_3d7_hb3_index
        File resource_vcf_hb3_dd2_index
        File resource_vcf_7g8_gb4_index

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 10*ceil(size([input_vcf, reference_fasta, resource_vcf_7g8_gb4, resource_vcf_hb3_dd2, resource_vcf_3d7_hb3], "GB"))

    command <<<
        ################################
        # Standard Preamble

        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        if [[ ${np} -gt 2 ]] ; then
            let np=${np}-1
        fi

        tot_mem_mb=$(free -m | grep '^Mem' | awk '{print $2}')

        ################################

        java_memory_size_mb=$((tot_mem_mb-5120))

        java -Xmx${java_memory_size_mb}M -jar /usr/GenomeAnalysisTK.jar \
            -T VariantRecalibrator \
                -R ~{reference_fasta} \
                -input ~{input_vcf} \
                -mode SNP \
                -resource:7g8_gb4,known=false,training=true,truth=true,prior=15.0 ~{resource_vcf_7g8_gb4} \
                -resource:hb3_dd2,known=false,training=true,truth=true,prior=15.0 ~{resource_vcf_hb3_dd2} \
                -resource:3d7_hb3,known=false,training=true,truth=true,prior=15.0 ~{resource_vcf_3d7_hb3} \
                -recalFile ~{prefix}_CombinedGVCFs.snp.recal \
                -tranchesFile ~{prefix}_CombinedGVCFs.snp.tranches \
                -rscriptFile ~{prefix}_CombinedGVCFs.snp.plots.R \
                -an QD \
                -an FS \
                -an SOR \
                -an DP \
                -an MQ \
                --maxGaussians 8 \
                --MQCapForLogitJitterTransform 70

        java -Xmx${java_memory_size_mb}M -jar /usr/GenomeAnalysisTK.jar \
            -T ApplyRecalibration \
                -R ~{reference_fasta} \
                -input ~{input_vcf} \
                --ts_filter_level 99.5 \
                -tranchesFile ~{prefix}_CombinedGVCFs.snp.tranches \
                -recalFile ~{prefix}_CombinedGVCFs.snp.recal \
                -mode SNP \
                -o ~{prefix}.snp.recalibrated.vcf

        java -Xmx${java_memory_size_mb}M -jar /usr/GenomeAnalysisTK.jar \
            -T VariantRecalibrator \
                -R ~{reference_fasta} \
                -input ~{prefix}.snp.recalibrated.vcf \
                -mode indel \
                -resource:7g8_gb4,known=false,training=true,truth=true,prior=12.0 ~{resource_vcf_7g8_gb4} \
                -resource:hb3_dd2,known=false,training=true,truth=true,prior=12.0 ~{resource_vcf_hb3_dd2} \
                -resource:3d7_hb3,known=false,training=true,truth=true,prior=12.0 ~{resource_vcf_3d7_hb3} \
                -recalFile ~{prefix}_CombinedGVCFs.snp.indel.recal \
                -tranchesFile ~{prefix}_CombinedGVCFs.snp.indel.tranches \
                -rscriptFile ~{prefix}_CombinedGVCFs.snp.indel.plots.R \
                -an QD \
                -an FS \
                -an MQ \
                --maxGaussians 4 \
                --MQCapForLogitJitterTransform 70

        java -Xmx${java_memory_size_mb}M -jar /usr/GenomeAnalysisTK.jar \
            -T ApplyRecalibration \
                -R ~{reference_fasta} \
                -input ~{prefix}.snp.recalibrated.vcf \
                --ts_filter_level 99.0 \
                -tranchesFile ~{prefix}_CombinedGVCFs.snp.indel.tranches \
                -recalFile ~{prefix}_CombinedGVCFs.snp.indel.recal \
                -mode indel \
                -o ~{prefix}.snp.indel.recalibrated.vcf

        java -Xmx${java_memory_size_mb}M -jar /usr/GenomeAnalysisTK.jar \
            -T VariantFiltration \
                -R ~{reference_fasta} \
                -V ~{prefix}.snp.indel.recalibrated.vcf \
                --filterExpression "VQSLOD <= 0.0"\
                --filterName "my_variant_filter" \
                -o ~{prefix}.snp.indel.recalibrated.filtered.vcf
    >>>

    output {
        File vcf = "~{prefix}.snp.indel.recalibrated.filtered.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "broadinstitute/gatk3:3.5-0"
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

task SortCompressIndexVcf {

    input {
        File input_vcf

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + 10*ceil(3*size(input_vcf, "GB"))

    String output_vcf = basename(input_vcf) + ".gz"

    command <<<
        ################################
        # Standard Preamble

        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        if [[ ${np} -gt 2 ]] ; then
            let np=${np}-1
        fi

        tot_mem_mb=$(free -m | grep '^Mem' | awk '{print $2}')

        ################################
        
        # Sort first because otherwise we'll end up with integers in the INFO fields again.
        bcftools sort -m$((tot_mem_mb-2048))M -o tmp.vcf ~{input_vcf}

        # Then we need to fix the integer values in the floating point INFO fields.
        # Without this fix / hack, downstream GATK3 tools will fail (specifically GenotypeGVCFs)
        awk -f - "tmp.vcf" > tmp2.vcf << 'AWK_CODE'
            BEGIN {
                FS = "\t"; OFS = "\t"

                # Define the set of ID keys that require floating point enforcement
                targets["BaseQRankSum"] = 1
                targets["ClippingRankSum"] = 1
                targets["ExcessHet"] = 1
                targets["HaplotypeScore"] = 1
                targets["InbreedingCoeff"] = 1
                targets["MLEAF"] = 1
                targets["MQ"] = 1
                targets["MQRankSum"] = 1
                targets["RAW_MQ"] = 1
                targets["ReadPosRankSum"] = 1
            }

            # Pass header lines through unchanged
            /^#/ { print; next }

            {
                # Column 8 is the INFO column
                # Split the INFO string by semicolon into an array
                n = split($8, info_fields, ";")
                
                new_info_str = ""

                for (i = 1; i <= n; i++) {
                    # Split Key=Value pairs
                    # We check if split returns 2 parts to avoid breaking on Boolean Flags
                    if (split(info_fields[i], kv, "=") == 2) {
                        key = kv[1]
                        val = kv[2]

                        # Check if this key is in our target list
                        if (key in targets) {
                            # Handle Number=A (comma-separated lists) like MLEAF
                            m = split(val, subvals, ",")
                            new_val_str = ""

                            for (j = 1; j <= m; j++) {
                                # Regex Check: Match strictly integers (optional - sign, digits only)
                                # This ignores values that are already floats (contain a dot)
                                if (subvals[j] ~ /^-?[0-9]+$/) {
                                    subvals[j] = subvals[j] ".0"
                                }
                                # Reconstruct comma-separated list
                                new_val_str = (j == 1 ? "" : new_val_str ",") subvals[j]
                            }
                            # Update the field with the new value
                            info_fields[i] = key "=" new_val_str
                        }
                    }
                    # Reconstruct the semicolon-separated INFO string
                    new_info_str = (i == 1 ? "" : new_info_str ";") info_fields[i]
                }

                # Replace the INFO column and print the line
                $8 = new_info_str
                print
            }
AWK_CODE
         
        ################################
        
        # Zip it:
        bgzip -c -l2 tmp2.vcf > ~{output_vcf}

        # Index the output:
        bcftools index --threads ${np} --tbi ~{output_vcf}
    >>>

    output {
        File vcf = output_vcf
        File vcf_index = "~{output_vcf}.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.2"
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