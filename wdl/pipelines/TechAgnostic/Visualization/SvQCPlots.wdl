version 1.0

import "../../../structs/Structs.wdl"

workflow PlotSVQCMetrics{

    input{
        String gcs_vcf_dir
        Array[String] samples
        File auxiliary_metrics
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
    #call addCoverageToSVstats{
    #    input:
    #        samples = samples,
    #        pbsv_stats = compileSVstats.pbsv_stat_out,
    #        sniffles_stats = compileSVstats.sniffles_stat_out,
    #        pav_stats = compileSVstats.pav_stat_out,
    #        auxiliary_metrics = auxiliary_metrics
    #}
    #

output{
        File pbsv_stat_out = compileSVstats.pbsv_stat_out
        File sniffles_stat_out = compileSVstats.sniffles_stat_out
#        File pav_stat_out = compileSVstats.pav_stat_out
#        File pbsv_all_stats_with_coverage = addCoverageToSVstats.pbsv_stat_out
#        File sniffles_all_stats_with_coverage = addCoverageToSVstats.sniffles_stat_out
#        File pav_all_stats_with_coverage = addCoverageToSVstats.pav_stat_out
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

    String pbsv_stat_out_name = sample_name + ".pbsv.svlen"
    String sniffles_stat_out_name = sample_name + ".sniffles.svlen"
#    String pav_stat_out_name = sample_basename+ ".pav.svlen"

    Int minimal_disk_size = (ceil(size(pbsv_vcf, "GB") + size(sniffles_vcf, "GB")  ) + 100 ) # 100GB buffer #+ size(pav_vcf, "GB")
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100


    command{
        set -euo pipefail

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
        set -euo pipefail

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
        Array[String] samples
        Array[File] pbsv_stats
        Array[File] sniffles_stats
#        Array[File] pav_stat_out
        RuntimeAttr? runtime_attr_override
    }

    Array[String] callers = ["pbsv", "sniffles"]
    Array[File] callers_stats = flatten([pbsv_stats, sniffles_stats])
    File sampleFile = write_lines(samples)

    Int minimal_disk_size = (ceil(size(callers_stats, "GB") ) + 100 ) # 100GB buffer #+ size(pav_stat_out, "GB")
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100


    command <<<
        set -euo pipefail

        mkdir -p stats_by_sample
        for file in ~{sep=" " callers_stats}
        do
            mv $file ./stats_by_sample
        done

        python <<CODE
import os
import subprocess

def main():
    svtypes=["ALL","DEL","DUP","CNV","INS","INV","OTH"]
    basedir=os.getcwd() + "/stats_by_sample"
    callers=["~{sep='", "' callers}"]

    samplefile = open('~{sampleFile}', 'r')
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
        # Count lines in the file
        with open(os.path.join(basedir, f"{caller}_stats", sample), 'r') as file:
            ALL = sum(1 for _ in file)

        counts_by_SV = subprocess.check_output(f"cut -f1 {os.path.join(basedir, f'{caller}_stats', sample)} | sort | uniq -c", shell=True)
        counts_by_SV = counts_by_SV.decode().split()

        SVs = {}
        for svtype in svtypes:
            SVs[svtype] = 0
        SVs['ALL'] = ALL

        # Process counts_by_SV
        for i in range(0, len(counts_by_SV), 2):
            num = int(counts_by_SV[i])
            SV = counts_by_SV[i + 1].upper()
            if SV in SVs:
                SVs[SV] = num
            else:
                SVs['OTH'] += num

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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:latest"
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


task addCoverageToSVstats{
    input{
        File auxiliary_metrics
        File pbsvStatsBySample
        File snifflesStatsBySample
#        Array[File] pav_stat_out
        RuntimeAttr? runtime_attr_override
    }

    Int minimal_disk_size = (ceil(size([pbsvStatsBySample, snifflesStatsBySample, auxiliary_metrics], "GB")  ) + 100 ) # 100GB buffer #+ size(pav_stat_out, "GB")
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100

    command<<<
        set -euo pipefail
        cut -f1,2 ~{auxiliary_metrics} | tail -n +2 > sample_cov
        echo $'sample\tcov' >> sample_cov

        sort -k1,1 sample_cov -o sample_cov
        sort -k1,1 ~{pbsvStatsBySample} -o pbsv_all_sample_stats
        sort -k1,1 ~{snifflesStatsBySample} -o sniffles_all_sample_stats

        join -1 1 -2 1 -a 1 -e 0 -o 1.1,2.2 sample_cov pbsv_all_sample_stats > pbsv_all_sample_stats_with_cov
        join -1 1 -2 1 -a 1 -e 0 -o 1.1,2.2 sample_cov sniffles_all_sample_stats > sniffles_all_sample_stats_with_cov

        # moving the header up to the top...
        tail -n 1 pbsv_all_sample_stats_with_cov > t1
        sed '$ d' pbsv_all_sample_stats_with_cov >> t1
        mv t1 pbsv_all_sample_stats_with_cov

        tail -n 1 sniffles_all_sample_stats_with_cov > t2
        sed '$ d' sniffles_all_sample_stats_with_cov >> t2
        mv t2 sniffles_all_sample_stats_with_cov

    >>>

    output{
        File pbsv_all_stats_with_cov = "pbsv_all_sample_stats_with_cov"
        File sniffles_all_stats_with_cov = "sniffles_all_sample_stats_with_cov"
#        File pav_all_stats_with_cov = "pav_all_sample_stats_with_cov"
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
