version 1.0

import "../../../structs/Structs.wdl"

import "../../../tasks/Utility/Finalize.wdl" as FF

workflow PlotSVQCMetrics{

    meta{
        description: "This workflow generates—for one or more samples—summaries and stats on the provided SV callers and their corresponding VCFs."
    }
    parameter_meta{
        samples: "List of sample names. Order must match the order of coverage_metrics."
        coverage_metrics: "List of coverage metrics for each sample. Order must match the order of samples."
        callers: "List of SV callers"
        reference_name: "Reference genome name"
        output_plot_notebook_name: "Name of the output plot notebook"
    }

    input{
        Array[String] samples # = "this.sample-t2ts.sample-t2t_id"
        Array[Float] coverage_metrics # = "this.sample-t2ts.coverage"
        Array[String] callers # = ["pav", "pbsv", "SNF"]
        Array[Array[File]] matching_vcfs # = ["this.sample-t2ts.pav_vcf", "this.sample-t2ts.pbsv_vcf", "this.sample-t2ts.sniffles_vcf"]

        String cohort_name
        String reference_name
        String output_plot_notebook_name

        String gcs_out_root_dir
    }

    # gather VCFs for each sample and caller
    scatter(i in range(length(callers))) {
        String caller = callers[i]
        scatter(j in range(length(samples))) {
            String sample = samples[j]
            File sv_vcf = matching_vcfs[i][j]
            call bcfQuerySV {
                input:
                    sample_name = sample,
                    input_vcf = sv_vcf,
                    caller = caller,
            }
        }
    }

    Array[File] all_SV_stats = flatten(bcfQuerySV.all_SV_stat_out)
    call concatSVstats {
        input:
            all_stats = all_SV_stats,
            callers = callers,
    }
    call compileSVstats {
        input:
            samples = samples,
            all_stats = all_SV_stats,
            callers = callers,
    }
    call addCoverageToSVstats {
        input:
            coverage_stats = coverage_metrics,
            samples = samples,
            allStatsBySample = compileSVstats.allStatsBySample,
            callers = callers,
    }
    call plotSVQCMetrics {
        input:
            all_stats_with_cov = addCoverageToSVstats.all_stats_with_cov,
            all_stats_by_type = concatSVstats.all_stats_by_type,
            callers = callers,
            output_file_name = output_plot_notebook_name,
            reference_name = reference_name,
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/CohortMetrics/~{cohort_name}/SvQCPlots"
    call FF.FinalizeToDir as SaveStatsBySVtype { input:
        files = concatSVstats.all_stats_by_type,
        outdir = outdir + "/StatsBySVtype"
    }
    call FF.FinalizeToDir as SaveStatsByCoverage { input:
        files = addCoverageToSVstats.all_stats_with_cov,
        outdir = outdir + "/StatsByCoverage"
    }
    call FF.FinalizeToDir as SavePlotPDFs { input:
        files = plotSVQCMetrics.output_pdfs,
        outdir = outdir + "/Plots"
    }
    call FF.FinalizeToFile as SaveNotebook { input:
        file = plotSVQCMetrics.out_plot_single_samples_stats,
        outdir = outdir,
        name = output_plot_notebook_name
    }
    call FF.FinalizeToFile as SaveHTML { input:
        file = plotSVQCMetrics.out_plot_single_samples_stats_html,
        outdir = outdir
    }

    output{
        Map[String, String] SvQCPlotsMisc = {
            "StatsBySVtype"  : SaveStatsBySVtype.gcs_dir,
            "StatsByCoverage": SaveStatsByCoverage.gcs_dir,
            "Plots"          : SavePlotPDFs.gcs_dir,
            "Notebook"       : SaveNotebook.gcs_path
        }
        File SvQCPlotsHTML = SaveHTML.gcs_path
    }
}


task bcfQuerySV{

    meta{
        description: "This task queries the input VCF file for SVs and outputs a file containing the SV stats."
    }
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

    command <<<
    set -euo pipefail
        bcftools query \
            -i '(INFO/SVLEN>49 || INFO/SVLEN<-49) && (FILTER=="PASS" || FILTER==".")' \
            --format "%SVTYPE\t%SVLEN\n" \
            ~{input_vcf} \
        > ~{sample_stat_out}
    >>>

    output{
        File all_SV_stat_out = sample_stat_out
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task concatSVstats{
    meta{
        description: "This task concatenates the SV stats for each sample and outputs a file containing the SV stats by type."
    }
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
        boot_disk_gb:       25,
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

    meta{
        description: "This task compiles the SV stats for each sample and outputs a file containing the SV stats by sample."
    }
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

        python3 /compilesvstats.py \
        --base_dir ./stats_by_sample \
        --sample_file ~{sampleFile} \
        --callers ~{sep=" " callers}  \
        --svtypes ALL DEL DUP CNV INS INV OTH

  >>>


    output {
        Array[File] allStatsBySample = glob("*_all_sample_stats")
    }
        #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-compile-sv-stats:0.0.0"
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

    meta{
        description: "This task adds coverage stats to the SV stats for each sample and outputs a file containing the SV stats with coverage."
    }
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
        boot_disk_gb:       25,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
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

    meta{
        description: "This task generates plots for the SV stats for each sample."
    }
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
        String output_file_name = "out_plot_single_samples_stats"
        RuntimeAttr? runtime_attr_override
    }
    Array[File] input_files = flatten([all_stats_with_cov, all_stats_by_type])
    String output_plot_notebook      = output_file_name + ".ipynb"
    String output_plot_notebook_html = output_file_name + ".html"

    Int minimal_disk_size = (ceil(size(input_files, "GB")  ) + 100 ) # 100GB buffer
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100

    command{
    set -euxo pipefail

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
        papermill /plot_single_sample_stats.ipynb \
            ~{output_plot_notebook} \
            -p reference_in ~{reference_name}  \
            -p callers_in "~{sep="," callers}"

        jupyter nbconvert \
            --to html \
            ~{output_plot_notebook}
    }
    output{
        File out_plot_single_samples_stats      = output_plot_notebook
        File out_plot_single_samples_stats_html = output_plot_notebook_html
        Array[File] output_pdfs = glob("*.pdf")
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       25,
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
