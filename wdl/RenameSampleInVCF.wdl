version 1.0

workflow RenameSampleInVCF {
    input {
        File vcf
        String sample
        String prefix
    }

    call RenameSample {
        input:
            vcf = vcf,
            sample = sample,
            prefix = prefix
    }

    output {
        File renamed_vcf = RenameSample.renamed_vcf
        File renamed_tbi = RenameSample.renamed_tbi
    }
}

task RenameSample {
    input {
        File vcf
        String sample
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(vcf, "GB")) + 5

    command <<<
        set -euxo pipefail

        bcftools query -l ~{vcf} > old_sample_name.txt
        N_SAMPLES=$(wc -l < old_sample_name.txt | tr -d ' ')
        if [ ${N_SAMPLES} -ne 1 ]; then
            echo "Expected exactly 1 sample in input VCF, found ${N_SAMPLES}:"
            cat old_sample_name.txt
            exit 1
        fi

        echo ~{sample} > new_sample_name.txt
        bcftools reheader --samples new_sample_name.txt -o ~{prefix}.renamed.vcf.gz ~{vcf}
        tabix -f ~{prefix}.renamed.vcf.gz
    >>>

    output {
        File renamed_vcf = "~{prefix}.renamed.vcf.gz"
        File renamed_tbi = "~{prefix}.renamed.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
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

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}
