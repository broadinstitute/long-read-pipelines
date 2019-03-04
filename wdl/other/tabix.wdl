
workflow TabixWorkflow {
    File input_vcf

    call Tabix {
        input: input_vcf=input_vcf
    }
}

task Tabix {
    File input_vcf
    Int disk_size = ceil(size(input_vcf, "GB") + 10)

    command {
        set -xe

        tabix ${input_vcf}
    }

    output {
        File output_bai="${input_vcf}.tbi"
    }

    runtime {
        docker: "weisburd/base-image-for-str-tools@sha256:bc02c67c69bbd13165ef2cd83f7b2ec87814d99fc007070cc7cd8865b29998a9"
        disks: "local-disk ${disk_size} HDD"
    }
}

