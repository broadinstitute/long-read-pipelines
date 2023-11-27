version 1.0

import "../../../structs/Structs.wdl"

workflow PlotSVQCMetrics{

    meta{
        description: "This workflow takes in a list of SV callers and a list of samples and runs the SV callers on each sample. It then concatenates the SV stats for each sample and plots the stats."
    }
    parameter_meta{
        gcs_vcf_dir: "GCS path to directory containing VCFs. Files should be named <sample_name>.<caller>.vcf.gz"
        samples: "List of sample names"
        coverage_metrics: "List of coverage metrics for each sample"
        callers: "List of SV callers"
        reference_name: "Reference genome name"
    }

    input{
        String gcs_vcf_dir
        Array[String] samples
        Array[Float] coverage_metrics
        Array[String] callers
        String reference_name = "GRCh38"
    }

    scatter(caller in callers){
        scatter(sample in samples){
            call bcfQuerySV{
                input:
                    sample_name = sample,
                    input_vcf = gcs_vcf_dir + "/" + sample + "." + caller + ".vcf.gz",
                    caller = caller,
            }
        }
    }
    Array[File] all_SV_stats = flatten(bcfQuerySV.all_SV_stat_out)

    call concatSVstats{
        input:
            all_stats = all_SV_stats,
            callers = callers,
    }

    call compileSVstats{
        input:
            samples = samples,
            all_stats = all_SV_stats,
            callers = callers,
    }
    call addCoverageToSVstats{
        input:
            coverage_stats = coverage_metrics,
            samples = samples,
            allStatsBySample = compileSVstats.allStatsBySample,
            callers = callers,
    }

    call plotSVQCMetrics{
        input:
            all_stats_with_cov = addCoverageToSVstats.all_stats_with_cov,
            all_stats_by_type = concatSVstats.all_stats_by_type,
            callers = callers,
            reference_name = reference_name,
    }


output{
        Array[File] all_stats_by_type = concatSVstats.all_stats_by_type
        Array[File] all_stats_with_cov = addCoverageToSVstats.all_stats_with_cov
        Array[File] metric_plot_pdfs = plotSVQCMetrics.output_pdfs
        File plot_notebook = plotSVQCMetrics.out_plot_single_sample_stats
    }


}


task bcfQuerySV{

    parameter_meta{
        sample_name: "Sample name"
        input_vcf: "Input vcf file"
        caller: "SV caller used to generate input vcf file"
    }

    input{
        String sample_name
        File input_vcf
        String caller
        RuntimeAttr? runtime_attr_override
    }

    String sample_stat_out = sample_name + "." + caller + ".txt"

    Int minimal_disk_size = (ceil(size(input_vcf, "GB") + 100 )) # 100GB buffer
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100

    command{
        set -euo pipefail

        cat ~{input_vcf} | bcftools query -i '(INFO/SVLEN>49 || INFO/SVLEN<-49) && FILTER=="PASS"' --format "%SVTYPE\t%SVLEN\n" > ~{sample_stat_out}

    }

    output{
        File all_SV_stat_out = sample_stat_out
    }
        #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
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

    parameter_meta{
        all_stats: "List of files containing SV stats for each sample"
        callers: "List of SV callers used to generate input vcf files"
    }

    input{
        Array[File] all_stats
        Array[String] callers
        RuntimeAttr? runtime_attr_override
    }

    Int minimal_disk_size = (ceil(size(all_stats, "GB"))  + 100 ) # 100GB buffer
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100

    command<<<
        set -euo pipefail

        for caller in ~{sep=" " callers}
        do
            for stat_file in ~{sep=" " all_stats}
            do
                if [[ ${stat_file} =~ ${caller} ]]
                then
                    cat ${stat_file} >> ${caller}_all_SV_lengths_by_type.txt
                fi
            done
        done

        # Check if the files created match the number of callers
        if [[ $(find ./*_all_SV_lengths_by_type.txt | wc -l) -ne ~{length(callers)} ]]
        then
            echo "ERROR: Number of callers does not match number of files created"
            ls *_all_SV_lengths_by_type.txt
            exit 1
        fi

    >>>

    output{
        Array[File] all_stats_by_type = glob("*_all_SV_lengths_by_type.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
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

    parameter_meta{
        samples: "List of sample names"
        all_stats: "List of files containing SV stats for each sample"
        callers: "List of SV callers used to generate input vcf files"
    }

    input {
        Array[String] samples
        Array[File] all_stats
        Array[String] callers
        RuntimeAttr? runtime_attr_override
    }

    File sampleFile = write_lines(samples)

    Int minimal_disk_size = (ceil(size(all_stats, "GB") ) + 100 ) # 100GB buffer
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100


    command <<<
        set -euo pipefail

        mkdir -p stats_by_sample
        for file in ~{sep=" " all_stats}
        do
            mv $file ./stats_by_sample
        done

        python3 <<CODE
import os
import subprocess

def main():
    svtypes = ["ALL", "DEL", "DUP", "CNV", "INS", "INV", "OTH"]
    basedir = os.getcwd() + "/stats_by_sample"
    callers = ["~{sep='", "' callers}"]

    samplefile = open('~{sampleFile}', 'r')
    samples = []
    for line in samplefile:
        samples.append(line.strip())
    samplefile.close()

    for caller in callers:
        outfile = open("%s/%s_all_sample_stats" % (os.getcwd(), caller), 'w')
        outfile.write("sample\t%s\n" % '\t'.join(svtypes))

        compile_stats(caller, svtypes, samples, basedir, outfile)

        outfile.close()

def pairwise(iterable):
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return zip(a, a)

def compile_stats(caller, svtypes, samples, basedir, outfile):
    for sample_name in samples:
        # Count lines in the file
        sample_svlen_file = os.path.join(basedir, f"{sample_name}.{caller}.txt")
        with open(sample_svlen_file, 'r') as file:
            ALL = sum(1 for _ in file)

        # Count SVs by type
        cut_sort_uniq_command = "cut -f1 {} | sort | uniq -c".format(sample_svlen_file)
        counts_by_sv = subprocess.check_output(cut_sort_uniq_command, shell=True)
        counts_by_sv_clean = counts_by_sv.decode().split()

        SVs = {}
        for svtype in svtypes:
            SVs[svtype] = 0
        SVs['ALL'] = ALL

        # Process counts_by_sv
        for i in range(0, len(counts_by_sv_clean), 2):
            num = int(counts_by_sv_clean[i])
            SV = counts_by_sv_clean[i + 1].upper()
            if SV in SVs:
                SVs[SV] = num
            else:
                SVs['OTH'] += num

        outfile.write("%s" % sample_name)

        for svtype in svtypes:
            outfile.write("\t%d" % SVs[svtype])
        outfile.write("\n")

if __name__ == "__main__":
    main()
CODE
  >>>


    output {
        Array[File] allStatsBySample = glob("*_all_sample_stats")
    }
        #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
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

    parameter_meta{
        coverage_stats: "List of coverage metrics for each sample"
        samples: "List of sample names"
        allStatsBySample: "List of files containing SV stats by sample"
        callers: "List of SV callers used to generate input vcf files"
    }

    input{
        Array[Float] coverage_stats
        Array[String] samples
        Array[File] allStatsBySample
        Array[String] callers
        RuntimeAttr? runtime_attr_override
    }

    Int minimal_disk_size = (ceil(size(allStatsBySample, "GB")  ) + 100 ) # 100GB buffer
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100

    command<<<
        set -euo pipefail
        # Combine the two arrays and write them to a file
        printf "sample\tCOV\n" > sample_cov
        samples=(~{sep="\t" samples})
        coverage_stats=(~{sep="\t" coverage_stats})
        paste -d $'\t' <(printf "%s\n" "${samples[@]}") <(printf "%s\n" "${coverage_stats[@]}") >> sample_cov

        sort -k1,1 sample_cov -o sample_cov

        cp ~{sep=" " allStatsBySample} .

        for caller in ~{sep=" " callers}
        do
            sort -k1,1 ${caller}_all_sample_stats -o ${caller}_all_sample_stats
            join -1 1 -2 1 -a 1 -e 0 -t $'\t' sample_cov ${caller}_all_sample_stats > ${caller}_all_sample_stats_with_cov
            { grep -m 1 '^sample' ${caller}_all_sample_stats_with_cov; grep -v '^sample' ${caller}_all_sample_stats_with_cov; } > ${caller}_all_sample_stats_with_cov.txt
        done

    >>>

    output{
        Array[File] all_stats_with_cov = glob("*_all_sample_stats_with_cov.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
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

task plotSVQCMetrics{

    parameter_meta{
        all_stats_with_cov: "List of files containing SV summary stats with coverage"
        all_stats_by_type: "List of files containing SV stats by type"
        callers: "List of SV callers used to generate input vcf files"
        reference_name: "Reference genome name"
    }

    input{
        Array[File] all_stats_with_cov
        Array[File] all_stats_by_type
        Array[String] callers
        String reference_name
        RuntimeAttr? runtime_attr_override
    }
    Array[File] input_files = flatten([all_stats_with_cov, all_stats_by_type])

    Int minimal_disk_size = (ceil(size(input_files, "GB")  ) + 100 ) # 100GB buffer
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100

    command{
        set -euo pipefail

        echo "Making directory for input files"
        mkdir ~{reference_name}
        echo "Current Directory"
        ls
        echo "Moving input files to ~{reference_name}"
        mv ~{sep=" " input_files} ~{reference_name}/
        echo "Current Directory After Moving Files"
        ls

        echo "Files in ~{reference_name} directory:"
        ls ~{reference_name}

        echo "Running jupyter notebook"
        papermill /plot_single_sample_stats.ipynb out_plot_single_sample_stats.ipynb \
        -p reference_in ~{reference_name}  \
        -p callers_in "~{sep="," callers}"
    }
    output{
        File out_plot_single_sample_stats = "out_plot_single_sample_stats.ipynb"
        Array[File] output_pdfs = glob("*.pdf")
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             24,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-plot-sv-metrics:beta.0.0.5"
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