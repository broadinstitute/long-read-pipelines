version 1.0

workflow SplitVCFWorkflow {
    input {
        File joint_vcf
        Array[String] regions
        Array[String] sample_names
    }

    scatter (region in regions) {
        scatter (samplename in sample_names) {
            call SplitVCFbySample {
                input:
                    joint_vcf = joint_vcf,
                    region = region,
                    sample_name = samplename
            }
        }
    }

    output {
        Array[File] single_sample_vcfs = flatten(SplitVCFbySample.single_sample_vcf)
        Array[File] single_sample_vcf_tbis = flatten(SplitVCFbySample.single_sample_vcf_tbi)
    }
}

task SplitVCFbySample {
    input {
        File joint_vcf
        String region
        String sample_name
    }

    command <<<
        set -euxo pipefail

        bcftools index ${joint_vcf}

        bcftools view -s ${sample_name} ${joint_vcf} -r ${region} -o ${sample_name}.subset.g.vcf.gz

        tabix -p vcf ${sample_name}.subset.g.vcf.gz
    >>>

    output {
        File single_sample_vcf = "${sample_name}.subset.g.vcf.gz"
        File single_sample_vcf_tbi = "${sample_name}.subset.g.vcf.gz.tbi"
    }

    runtime {
        cpu: 1
        memory: "64 GiB"
        disks: "local-disk ${disk_size} HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }

    Int disk_size = ceil(2 * size(joint_vcf, "GiB")) + 1
}
