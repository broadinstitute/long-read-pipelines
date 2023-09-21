version 1.0

import "../../../structs/Structs.wdl"

workflow PlotSVQCMetrics{

    input{
        String gcs_vcf_dir
        Array[String] samples
    }
    scatter(sample in samples){
        call bcfQuerySV{
            input:
                sample_name = sample,
                pbsv_vcf = gcs_vcf_dir + "/" + sample + ".pbsv.vcf.gz",
                sniffles_vcf = gcs_vcf_dir + "/" + sample + ".sniffles.vcf.gz",
#                pav_vcf = input_vcf + ".pav.vcf.gz"
        }
    }

    call concatSVstats{
        input:
            pbsv_stats = bcfQuerySV.pbsv_stat_out,
            sniffles_stats = bcfQuerySV.sniffles_stat_out,
    }

    call compileSVstats{
        input:
            samples = samples,
            pbsv_stats = bcfQuerySV.pbsv_stat_out,
            sniffles_stats = bcfQuerySV.sniffles_stat_out,
    }
    #call addcoverageinfo
    #
}


task bcfQuerySV{
    input{
        String sample_name
        File pbsv_vcf
        File sniffles_vcf
#        File pav_vcf
        RuntimeAttr? runtime_attr_override
    }

    String pbsv_stat_out_name = sample_name + ".pbsv.svlen"
    String sniffles_stat_out_name = sample_name + ".sniffles.svlen"
#    String pav_stat_out_name = sample_basename+ ".pav.svlen"

    Int minimal_disk_size = (ceil(size(pbsv_vcf, "GB") + size(sniffles_vcf, "GB")  ) + 100 ) # 100GB buffer #+ size(pav_vcf, "GB")
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100


    command{
        cat ~{pbsv_vcf} | bcftools query -i '(INFO/SVLEN>49 || INFO/SVLEN<-49) && FILTER=="PASS"' --format "%SVTYPE\t%SVLEN\n" > ~{pbsv_stat_out_name}
        cat ~{sniffles_vcf} | bcftools query -i '(INFO/SVLEN>49 || INFO/SVLEN<-49) && FILTER=="PASS"' --format "%SVTYPE\t%SVLEN\n" > ~{sniffles_stat_out_name}

    }
#    cat ~{pav_vcf} | bcftools query -i '(INFO/SVLEN>49 || INFO/SVLEN<-49) && FILTER=="PASS"' --format "%SVTYPE\t%SVLEN\n" > ~{pav_stat_out}
    output{
        File pbsv_stat_out = pbsv_stat_out_name
        File sniffles_stat_out = sniffles_stat_out_name
#        File pav_stat_out = pav_stat_out_name
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

task concatSVstats{
    input{
        Array[File] pbsv_stats
        Array[File] sniffles_stats
#        Array[File] pav_stat_out
        RuntimeAttr? runtime_attr_override
    }

    Int minimal_disk_size = (ceil(size(pbsv_stats, "GB") + size(sniffles_stats, "GB")  ) + 100 ) # 100GB buffer #+ size(pav_stat_out, "GB")
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100

    command<<<
        for i in ~{sep=" " pbsv_stats}
        do
            cat ${i} >> pbsv_all_SV_lengths_by_type.svlen
        done

        for i in ~{sep=" " sniffles_stats}
        do
            cat ${i} >> sniffles_all_SV_lengths_by_type.svlen
        done
    >>>
#        for i in ~{sep=" " pav_stat_out}
#        do
#            cat ${i} >> pav_all_SV_lengths_by_type.svlen
#        done
    output{
        File all_pbsv_stats = "pbsv_all_SV_lengths_by_type.svlen"
        File all_sniffles_stats = "sniffles_all_SV_lengths_by_type.svlen"
#        File all_pav_stats = "pav_all_SV_lengths_by_type.svlen"
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

task compileSVstats {
    input {
        Array[File] samples
        Array[File] pbsv_stats
        Array[File] sniffles_stats
#        Array[File] pav_stat_out
        RuntimeAttr? runtime_attr_override
    }

    Array[String] callers = ["pbsv", "sniffles"]
    Array[File] callers_stats = flatten([pbsv_stats, sniffles_stats])
    File sampleFile = write_lines(samples)

    Int disk_size = (ceil(size(callers_stats, "GB") + size(samples, "GB")  ) + 100 ) # 100GB buffer

    # callers variable bellow does not include pav

    command <<<
        mkdir -p stats_by_sample
        for file in ~{sep=" " callers_stats}
        do
            mv file ./stats_by_sample
        done

        python3 <<CODE
import os

def main():
    svtypes=["ALL","DEL","DUP","CNV","INS","INV","OTH"]
    basedir=os.getcwd() + "/stats_by_sample"
    callers=["~{sep='", "' callers}"]

    samplefile = open(~{sampleFile}, 'r')
    samples = []
    for line in samplefile:
        samples.append(line.strip())
    samplefile.close()

    for caller in callers:
        outfile=open("%s/%s_all_sample_stats" % (basedir, caller),'w')
        outfile.write("sample\t%s\n" % '\t'.join(svtypes))

        compile_stats(caller, svtypes, samples, basedir, outfile)

        outfile.close()

def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a)

def compile_stats(caller, svtypes, samples, basedir, outfile):
    for sample in samples:
        ALL=!wc -l {basedir}/{caller}_stats/{sample}|cut -f1 -d' '
        ALL=int(ALL.spstr)
        counts_by_SV=!echo `cut -f1 {basedir}/{caller}_stats/{sample}|sort|uniq -c`
        counts_by_SV=counts_by_SV.spstr.split()

        SVs = {}
        for svtype in svtypes:
            SVs[svtype] = 0
        SVs['ALL'] = ALL

        for num, SV in pairwise(counts_by_SV):
            SV=SV.upper()
            if SV in SVs:
                SVs[SV] = int(num)
            else:
                SVs['OTH'] += int(num)
        outfile.write("%s" % sample)

        for svtype in svtypes:
            outfile.write("\t%d" % SVs[svtype])
        outfile.write("\n")

if __name__ == "__main__":
    main()
CODE
  >>>


    output {
#        Array[File] output_stats = glob("*_all_sample_stats")
#        File pavStatsbysample = "pav_all_sample_stats"
        File pbsvStatsBySample = "pbsv_all_sample_stats"
        File snifflesStatsBySample = "sniffles_all_sample_stats"
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

