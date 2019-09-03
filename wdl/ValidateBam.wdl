version 1.0

# TODO: describe purpose
task ValidateBam {
    input {
        File input_bam
        String docker
    }

    Int cpus = 2
    Int disk_size = ceil(size(input_bam, "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx4g -jar /usr/local/bin/gatk.jar ValidateSamFile -I ~{input_bam} -O bam_validation_report.txt
    >>>

    output {
        File report = "bam_validation_report.txt"
    }

    runtime {
        cpu: "~{cpus}"
        memory: "8G"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
        maxRetries: 0
        docker: "~{docker}"
    }
}