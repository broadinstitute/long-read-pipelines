version 1.0

import "../../structs/Structs.wdl"
import "../../tasks/Utility/Finalize.wdl" as FF

workflow ONTPfTypeDrugResistanceMarkers {
    input {
        File vcf

        String dir_prefix
        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ONTPfTypeDrugResistanceMarkers/~{dir_prefix}"

    call AnnotateEffectsOfSelectedVariants { input: vcf = vcf }

    # Finalize data
    String dir = outdir + "/reports"

    call FF.FinalizeToFile as FinalizeDRReport { input: outdir = dir, file = AnnotateEffectsOfSelectedVariants.report }

    output {
        File drug_res_report = FinalizeDRReport.gcs_path
    }
}

task AnnotateEffectsOfSelectedVariants {
    input {
        File vcf

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 2*ceil(size(vcf, "GB"))
    String base = basename(vcf, ".vcf.gz")

    command <<<
        set -x

        zcat ~{vcf} | \
            sed 's/^Pf3D7_0//' | \
            sed 's/^Pf3D7_1/1/' | \
            sed 's/_v3\t/\t/' | \
            awk '{ if ($0 ~ "^#" || (length($4) == 1 && length($5) == 1 && $7 == "PASS")) print $0 }' \
            > reformatted.vcf

        /usr/local/bin/snpEff ann -v Plasmodium_falciparum reformatted.vcf \
            > ann.vcf

        grep PF3D7_0417200 ann.vcf | grep p.Cys50Arg | wc -l | awk '{ if ($1 > 0) print "pfdhfr\tPF3D7_0417200\tp.Cys50Arg\t+"; else print "pfdhfr\tPF3D7_0417200\tp.Cys50Arg\t-" }' > drug_resistance_report.txt
        grep PF3D7_0417200 ann.vcf | grep p.Asn51Ile | wc -l | awk '{ if ($1 > 0) print "pfdhfr\tPF3D7_0417200\tp.Asn51Ile\t+"; else print "pfdhfr\tPF3D7_0417200\tp.Asn51Ile\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0417200 ann.vcf | grep p.Cys59Arg | wc -l | awk '{ if ($1 > 0) print "pfdhfr\tPF3D7_0417200\tp.Cys59Arg\t+"; else print "pfdhfr\tPF3D7_0417200\tp.Cys59Arg\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0417200 ann.vcf | grep p.Ser108Asn | wc -l | awk '{ if ($1 > 0) print "pfdhfr\tPF3D7_0417200\tp.Ser108Asn\t+"; else print "pfdhfr\tPF3D7_0417200\tp.Ser108Asn\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0417200 ann.vcf | grep p.Ile164Lys | wc -l | awk '{ if ($1 > 0) print "pfdhfr\tPF3D7_0417200\tp.Ile164Lys\t+"; else print "pfdhfr\tPF3D7_0417200\tp.Ile164Lys\t-" }' >> drug_resistance_report.txt

        grep PF3D7_0523000 ann.vcf | grep p.Asn86Tyr | wc -l | awk '{ if ($1 > 0) print "pfmdr1\tPF3D7_0523000\tp.Asn86Tyr\t+"; else print "pfmdr1\tPF3D7_0523000\tp.Asn86Tyr\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0523000 ann.vcf | grep p.Tyr184Phe | wc -l | awk '{ if ($1 > 0) print "pfmdr1\tPF3D7_0523000\tp.Tyr184Phe\t+"; else print "pfmdr1\tPF3D7_0523000\tp.Tyr184Phe\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0523000 ann.vcf | grep p.Ser1034Cys | wc -l | awk '{ if ($1 > 0) print "pfmdr1\tPF3D7_0523000\tp.Ser1034Cys\t+"; else print "pfmdr1\tPF3D7_0523000\tp.Ser1034Cys\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0523000 ann.vcf | grep p.Asn1024Asp | wc -l | awk '{ if ($1 > 0) print "pfmdr1\tPF3D7_0523000\tp.Asn1024Asp\t+"; else print "pfmdr1\tPF3D7_0523000\tp.Asn1024Asp\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0523000 ann.vcf | grep p.Asp1246Tyr | wc -l | awk '{ if ($1 > 0) print "pfmdr1\tPF3D7_0523000\tp.Asp1246Tyr\t+"; else print "pfmdr1\tPF3D7_0523000\tp.Asp1246Tyr\t-" }' >> drug_resistance_report.txt

        grep PF3D7_0709000 ann.vcf | grep p.Lys76Thr | wc -l | awk '{ if ($1 > 0) print "pfcrt\tPF3D7_0709000\tp.Lys76Thr\t+"; else print "pfcrt\tPF3D7_0709000\tp.Lys76Thr\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0709000 ann.vcf | grep p.Met74Ile | wc -l | awk '{ if ($1 > 0) print "pfcrt\tPF3D7_0709000\tp.Met74Ile\t+"; else print "pfcrt\tPF3D7_0709000\tp.Met74Ile\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0709000 ann.vcf | grep p.Asn75Glu | wc -l | awk '{ if ($1 > 0) print "pfcrt\tPF3D7_0709000\tp.Asn75Glu\t+"; else print "pfcrt\tPF3D7_0709000\tp.Asn75Glu\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0709000 ann.vcf | grep p.Cys72Ser | wc -l | awk '{ if ($1 > 0) print "pfcrt\tPF3D7_0709000\tp.Lys76Thr\t+"; else print "pfcrt\tPF3D7_0709000\tp.Lys76Thr\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0709000 ann.vcf | grep p.His97Tyr | wc -l | awk '{ if ($1 > 0) print "pfcrt\tPF3D7_0709000\tp.His97Tyr\t+"; else print "pfcrt\tPF3D7_0709000\tp.His97Tyr\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0709000 ann.vcf | grep p.Cys101Phe | wc -l | awk '{ if ($1 > 0) print "pfcrt\tPF3D7_0709000\tp.Cys101Phe\t+"; else print "pfcrt\tPF3D7_0709000\tp.Cys101Phe\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0709000 ann.vcf | grep p.Phe145Ile | wc -l | awk '{ if ($1 > 0) print "pfcrt\tPF3D7_0709000\tp.Phe145Ile\t+"; else print "pfcrt\tPF3D7_0709000\tp.Phe145Ile\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0709000 ann.vcf | grep p.Met343Leu | wc -l | awk '{ if ($1 > 0) print "pfcrt\tPF3D7_0709000\tp.Met343Leu\t+"; else print "pfcrt\tPF3D7_0709000\tp.Met343Leu\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0709000 ann.vcf | grep p.Ser350Arg | wc -l | awk '{ if ($1 > 0) print "pfcrt\tPF3D7_0709000\tp.Ser350Arg\t+"; else print "pfcrt\tPF3D7_0709000\tp.Ser350Arg\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0709000 ann.vcf | grep p.Gly353Val | wc -l | awk '{ if ($1 > 0) print "pfcrt\tPF3D7_0709000\tp.Glu353Val\t+"; else print "pfcrt\tPF3D7_0709000\tp.Glu353Val\t-" }' >> drug_resistance_report.txt

        grep PF3D7_0810800 ann.vcf | grep p.Ser436Ala | wc -l | awk '{ if ($1 > 0) print "pfdhps\tPF3D7_0810800\tp.Ser436Ala\t+"; else print "pfdhps\tPF3D7_0810800\tp.Ser436Ala\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0810800 ann.vcf | grep p.Lys437Gly | wc -l | awk '{ if ($1 > 0) print "pfdhps\tPF3D7_0810800\tp.Lys437Gly\t+"; else print "pfdhps\tPF3D7_0810800\tp.Lys437Gly\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0810800 ann.vcf | grep p.Lys540Glu | wc -l | awk '{ if ($1 > 0) print "pfdhps\tPF3D7_0810800\tp.Lys540Glu\t+"; else print "pfdhps\tPF3D7_0810800\tp.Lys540Glu\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0810800 ann.vcf | grep p.Ala581Gly | wc -l | awk '{ if ($1 > 0) print "pfdhps\tPF3D7_0810800\tp.Ala581Gly\t+"; else print "pfdhps\tPF3D7_0810800\tp.Ala581Gly\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0810800 ann.vcf | grep p.Ala613Thr | wc -l | awk '{ if ($1 > 0) print "pfdhps\tPF3D7_0810800\tp.Ala613Thr\t+"; else print "pfdhps\tPF3D7_0810800\tp.Ala613Thr\t-" }' >> drug_resistance_report.txt
        grep PF3D7_0810800 ann.vcf | grep p.Ala613Ser | wc -l | awk '{ if ($1 > 0) print "pfdhps\tPF3D7_0810800\tp.Ala613Ser\t+"; else print "pfdhps\tPF3D7_0810800\tp.Ala613Ser\t-" }' >> drug_resistance_report.txt

        grep PF3D7_1343700 ann.vcf | grep p.Tyr493His | wc -l | awk '{ if ($1 > 0) print "pfkelch13\tPF3D7_1343700\tp.Tyr493His\t+"; else print "pfkelch13\tPF3D7_1343700\tp.Tyr493His\t-" }' >> drug_resistance_report.txt
        grep PF3D7_1343700 ann.vcf | grep p.Arg539Thr | wc -l | awk '{ if ($1 > 0) print "pfkelch13\tPF3D7_1343700\tp.Arg539Thr\t+"; else print "pfkelch13\tPF3D7_1343700\tp.Arg539Thr\t-" }' >> drug_resistance_report.txt
        grep PF3D7_1343700 ann.vcf | grep p.Ile543Thr | wc -l | awk '{ if ($1 > 0) print "pfkelch13\tPF3D7_1343700\tp.Ile543Thr\t+"; else print "pfkelch13\tPF3D7_1343700\tp.Ile543Thr\t-" }' >> drug_resistance_report.txt
        grep PF3D7_1343700 ann.vcf | grep p.Arg561His | wc -l | awk '{ if ($1 > 0) print "pfkelch13\tPF3D7_1343700\tp.Arg561His\t+"; else print "pfkelch13\tPF3D7_1343700\tp.Arg561His\t-" }' >> drug_resistance_report.txt
        grep PF3D7_1343700 ann.vcf | grep p.Cys580Tyr | wc -l | awk '{ if ($1 > 0) print "pfkelch13\tPF3D7_1343700\tp.Cys580Tyr\t+"; else print "pfkelch13\tPF3D7_1343700\tp.Cys580Tyr\t-" }' >> drug_resistance_report.txt
        grep PF3D7_1343700 ann.vcf | grep p.Ala675Val | wc -l | awk '{ if ($1 > 0) print "pfkelch13\tPF3D7_1343700\tp.Ala675Val\t+"; else print "pfkelch13\tPF3D7_1343700\tp.Ala675Val\t-" }' >> drug_resistance_report.txt
        grep PF3D7_1343700 ann.vcf | grep p.Phe446Ile | wc -l | awk '{ if ($1 > 0) print "pfkelch13\tPF3D7_1343700\tp.Phe446Ile\t+"; else print "pfkelch13\tPF3D7_1343700\tp.Phe446Ile\t-" }' >> drug_resistance_report.txt
        grep PF3D7_1343700 ann.vcf | grep p.Met476Ile | wc -l | awk '{ if ($1 > 0) print "pfkelch13\tPF3D7_1343700\tp.Met476Ile\t+"; else print "pfkelch13\tPF3D7_1343700\tp.Met476Ile\t-" }' >> drug_resistance_report.txt
        grep PF3D7_1343700 ann.vcf | grep p.Asn458Tyr | wc -l | awk '{ if ($1 > 0) print "pfkelch13\tPF3D7_1343700\tp.Asn458Tyr\t+"; else print "pfkelch13\tPF3D7_1343700\tp.Asn458Tyr\t-" }' >> drug_resistance_report.txt
        grep PF3D7_1343700 ann.vcf | grep p.Phe553Leu | wc -l | awk '{ if ($1 > 0) print "pfkelch13\tPF3D7_1343700\tp.Phe553Leu\t+"; else print "pfkelch13\tPF3D7_1343700\tp.Phe553Leu\t-" }' >> drug_resistance_report.txt
        grep PF3D7_1343700 ann.vcf | grep p.Phe574Leu | wc -l | awk '{ if ($1 > 0) print "pfkelch13\tPF3D7_1343700\tp.Phe574Leu\t+"; else print "pfkelch13\tPF3D7_1343700\tp.Phe574Leu\t-" }' >> drug_resistance_report.txt
        grep PF3D7_1343700 ann.vcf | grep p.Arg633Ile | wc -l | awk '{ if ($1 > 0) print "pfkelch13\tPF3D7_1343700\tp.Arg633Ile\t+"; else print "pfkelch13\tPF3D7_1343700\tp.Arg633Ile\t-" }' >> drug_resistance_report.txt
    >>>

    output {
        File report = "drug_resistance_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/biocontainers/snpeff:5.1d--hdfd78af_0"
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
