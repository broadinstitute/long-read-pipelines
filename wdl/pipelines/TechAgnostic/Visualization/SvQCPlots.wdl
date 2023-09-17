version 1.0

import "../../../structs/Structs.wdl"

workflow PlotSVQCMetrics{

    input{
        String gcs_dir_to_vcf
        Array[String] samples
    }
    scatter(sample in samples){
        call bcfQuerySV{
            input:
                sample_name = sample,
                pbsv_vcf = gcs_dir_to_vcf + "/" + sample + ".pbsv.vcf.gz",
                sniffles_vcf = gcs_dir_to_vcf + "/" + sample + ".sniffles.vcf.gz",
#                pav_vcf = input_vcf + ".pav.vcf.gz"
        }
    }
}


task bcfQuerySV{
    input{
        String sample_name
        File pbsv_vcf
        File sniffles_vcf
#        File pav_vcf
        RuntimeAttr? runtime_attr_override
    }

    String pbsv_stat_out = sample_name + ".pbsv.svlen"
    String sniffles_stat_out = sample_name + ".sniffles.svlen"
#    String pav_stat_out = sample_basename+ ".pav.svlen"

    Int minimal_disk_size = (ceil(size(pbsv_vcf, "GB") + size(sniffles_vcf, "GB")  ) + 100 ) # 100GB buffer #+ size(pav_vcf, "GB")
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100


    command{
        cat ~{pbsv_vcf} | bcftools query -i '(INFO/SVLEN>49 || INFO/SVLEN<-49) && FILTER=="PASS"' --format "%SVTYPE\t%SVLEN\n" > ~{pbsv_stat_out}
        cat ~{sniffles_vcf} | bcftools query -i '(INFO/SVLEN>49 || INFO/SVLEN<-49) && FILTER=="PASS"' --format "%SVTYPE\t%SVLEN\n" > ~{sniffles_stat_out}

    }
#    cat ~{pav_vcf} | bcftools query -i '(INFO/SVLEN>49 || INFO/SVLEN<-49) && FILTER=="PASS"' --format "%SVTYPE\t%SVLEN\n" > ~{pav_stat_out}
    output{
        File pbsv_stat_out = pbsv_stat_out
        File sniffles_stat_out = sniffles_stat_out
#        File pav_stat_out = pav_stat_out
    }
        #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             24,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
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

