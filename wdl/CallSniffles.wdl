version 1.0

import "tasks/Sniffles2.wdl" as Sniffles

workflow CallSniffles {
    input {
        File bam
        File bai
        String prefix
        File tandem_repeats_bed
    }

    call Sniffles.Sniffles {
        input:
            bam = bam,
            bai = bai,
            prefix = prefix,
            tandem_repeats_bed = tandem_repeats_bed
    }

    output {
        File sniffles_vcf = Sniffles.vcf
    }
}