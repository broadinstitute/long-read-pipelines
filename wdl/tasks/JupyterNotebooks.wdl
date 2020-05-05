version 1.0

import "Structs.wdl"

task PB10xMasSeqSingleFlowcellReport {

    meta {
        description : "Create a report for a given MASSeq run which summarizes the results  using a given Jupyter Notebook template."
        author : "Jonn Smith"
        email : "jonn@broadinstitute.org"
    }

    input {
        File notebook_template

        String sample_name

        File subreads_stats
        File ccs_reads_stats
        File array_elements_stats
        File ccs_report_file

        File raw_ccs_bam_file
        File array_element_bam_file
        File array_elements_genome_aligned
        File ccs_rejected_bam_file

        File annotated_bam_file

        File longbow_passed_reads_file
        File longbow_failed_reads_file

        File longbow_passed_ccs_reads
        File longbow_failed_ccs_reads
        File ccs_reclaimable_reads
        File ccs_reclaimed_reads
        File ccs_rejected_longbow_failed_reads
        File raw_array_elements
        File ccs_reclaimed_array_elements

        File zmw_stats_json_gz

        File? zmw_subread_stats_file
        File? polymerase_read_lengths_file
        File? approx_raw_subread_array_lengths

        File? ten_x_metrics_file
        String mas_seq_model

        File workflow_dot_file

        String prefix = ""
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        notebook_template : "Jupyter notebook MASSeq template to run with the given data to produce a MASSeq report."

        sample_name : "Name of the MAS-seq sample being analyzed in this report."

        subreads_stats : "Samtools stats file created from the raw subreads from the PacBio instrument."
        ccs_reads_stats : "Samtools raw stats file created from the aligned CCS corrected reads from the PacBio instrument."
        array_elements_stats : "Samtools raw stats file created from the aligned MASSeq array elements."

        ccs_report_file : "CCS report file from the CCS run for the data from the PacBio instrument."
        raw_ccs_bam_file : "Unaligned reads file in BAM format from the CCS process (pre-array splitting)."

        array_element_bam_file : "Transcriptome aligned reads file in BAM format containing aligned MASSeq array elements as individual reads."
        array_elements_genome_aligned : "Genome aligned reads file in BAM format containing aligned MASSeq array elements as individual reads."
        ccs_rejected_bam_file : "Bam file containing all subreads from zmws that were rejected by CCS."

        annotated_bam_file : "Bam file containing ccs corrected reads with annotated sections in the SG tag."

        longbow_passed_reads_file : "Bam file containing all reads that passed the longbow filter for the model used in this run (both ccs passed and reclaimed)."
        longbow_failed_reads_file : "Bam file containing alll reads that failed the longbow filter for the model used in this run (both ccs passed and reclaimed)."

        longbow_passed_ccs_reads : "Bam file containing ccs corrected reads that passed the longbow filter for the model used in this run (CCS Corrected reads ONLY)."
        longbow_failed_ccs_reads : "Bam file containing ccs corrected reads that failed the longbow filter for the model used in this run (CCS Corrected reads ONLY)."
        ccs_reclaimable_reads : "Bam file containing ccs rejected reads that are deemed to be reclaimable."
        ccs_reclaimed_reads : "Bam file containing ccs rejected reads that have been reclaimed."
        ccs_rejected_longbow_failed_reads : "Bam file containing ccs reclaimable reads that did not pass longbow filtering and were not reclaimed."
        raw_array_elements : "Bam file containing the raw unaligned array elements created from the longbow_passed_reads_file."
        ccs_reclaimed_array_elements : "Bam file containing the unaligned array elements created from reclaimed CCS reads."

        zmw_stats_json_gz : "ZMW stats json.gz file from the PacBio instrument."

        zmw_subread_stats_file : "[optional] File containing statistics about the subreads from each ZMW (created by collect_zmw_subread_stats.py in the PBUtils docker container)."
        polymerase_read_lengths_file : "[optional] File containing the lengths of each polymerase read from the sequencer (as created by collect_polymerase_read_lengths.py)"
        approx_raw_subread_array_lengths : "[optional] File containing the approximate array length information from the raw (pre-ccs) subreads file  (created by get_approx_raw_subread_array_lengths.py in the Cartographer docker container)."

        ten_x_metrics_file : "[optional] Stats file from the 10x tool run for the data in this MASSeq run.  If not supplied stats will not be displayed in the resulting report."
        mas_seq_model : "Built-in mas-seq model to use."

        workflow_dot_file : "DOT file containing the representation of this workflow used to create and analyze the data.  This is included in the QC reports (the DOT file can be generated with womtool)."

        prefix : "[optional] Prefix to prepend to the name of the generated report."
        runtime_attr_override : "[optional] Runtime attributes struct with which to override the docker container runtime.."
    }

    String nb_name = prefix + "report.ipynb"
    String html_out = prefix + "report.html"
    String pdf_out = prefix + "report.pdf"

    Int disk_size = 20 + 8*ceil((
            size(notebook_template, "GB") +
            size(subreads_stats, "GB") +
            size(ccs_reads_stats, "GB") +
            size(ccs_report_file, "GB") +
            size(raw_ccs_bam_file, "GB") +
            size(array_element_bam_file, "GB") +
            size(ccs_rejected_bam_file, "GB") +
            size(annotated_bam_file, "GB") +
            size(raw_ccs_bam_file, "GB") +
            size(zmw_subread_stats_file, "GB") +
            size(polymerase_read_lengths_file, "GB") +
            size(ten_x_metrics_file, "GB") +
            size(workflow_dot_file, "GB")
        ))

    # Handle the optional files:
    String ten_x_metrics_file_flag = if defined(ten_x_metrics_file) then "true" else "false"
    String zmw_subread_stats_file_flag = if defined(zmw_subread_stats_file) then "true" else "false"
    String polymerase_read_lengths_file_flag = if defined(polymerase_read_lengths_file) then "true" else "false"
    String approx_raw_subread_array_lengths_flag = if defined(approx_raw_subread_array_lengths) then "true" else "false"

    command <<<
        set -euxo pipefail

        # Set up memory logging daemon:
        MEM_LOG_INTERVAL_s=5
        DO_MEMORY_LOG=true
        while $DO_MEMORY_LOG ; do
            date
            date +%s
            cat /proc/meminfo
            sleep $MEM_LOG_INTERVAL_s
        done >> memory_log.txt &
        mem_pid=$!

        # Copy the notebook template to our current folder:
        cp ~{notebook_template} ~{nb_name}

        # Create a template to create the html report with collapsed code:
        echo "{%- extends 'full.tpl' -%}" > hidecode.tpl
        echo "" >> hidecode.tpl
        echo "{% block input_group %}" >> hidecode.tpl
        echo "    {%- if cell.metadata.get('nbconvert', {}).get('show_code', False) -%}" >> hidecode.tpl
        echo "        ((( super() )))" >> hidecode.tpl
        echo "    {%- endif -%}" >> hidecode.tpl
        echo "{% endblock input_group %}" >> hidecode.tpl

        # Set some environment variables for the notebook to read in:
        export DATE_RUN="$(date)"
        export WDL_NAME="PB10xMasSeqArraySingleFlowcell.wdl"
        export REPO_INFO="git@github.com:broadinstitute/long-read-pipelines.git"

        # Prepare the config file:
        rm -f mas-seq_qc_inputs.config

        echo "~{sample_name}" >> mas-seq_qc_inputs.config

        echo "~{subreads_stats}" >> mas-seq_qc_inputs.config
        echo "~{ccs_reads_stats}" >> mas-seq_qc_inputs.config
        echo "~{array_elements_stats}" >> mas-seq_qc_inputs.config
        echo "~{ccs_report_file}" >> mas-seq_qc_inputs.config

        echo "~{raw_ccs_bam_file}" >> mas-seq_qc_inputs.config
        echo "~{array_element_bam_file}" >> mas-seq_qc_inputs.config
        echo "~{array_elements_genome_aligned}" >> mas-seq_qc_inputs.config
        echo "~{ccs_rejected_bam_file}" >> mas-seq_qc_inputs.config

        echo "~{annotated_bam_file}" >> mas-seq_qc_inputs.config

        echo "~{longbow_passed_reads_file}" >> mas-seq_qc_inputs.config
        echo "~{longbow_failed_reads_file}" >> mas-seq_qc_inputs.config

        echo "~{longbow_passed_ccs_reads}" >> mas-seq_qc_inputs.config
        echo "~{longbow_failed_ccs_reads}" >> mas-seq_qc_inputs.config
        echo "~{ccs_reclaimable_reads}" >> mas-seq_qc_inputs.config
        echo "~{ccs_reclaimed_reads}" >> mas-seq_qc_inputs.config
        echo "~{ccs_rejected_longbow_failed_reads}" >> mas-seq_qc_inputs.config
        echo "~{raw_array_elements}" >> mas-seq_qc_inputs.config
        echo "~{ccs_reclaimed_array_elements}" >> mas-seq_qc_inputs.config

        echo "~{zmw_stats_json_gz}" >> mas-seq_qc_inputs.config

        if ~{zmw_subread_stats_file_flag} ; then
            echo "~{zmw_subread_stats_file}" >> mas-seq_qc_inputs.config
        else
            echo "NON-EXISTENT-PLACEHOLDER" >> mas-seq_qc_inputs.config
        fi
        if ~{polymerase_read_lengths_file_flag} ; then
            echo "~{polymerase_read_lengths_file}" >> mas-seq_qc_inputs.config
        else
            echo "NON-EXISTENT-PLACEHOLDER" >> mas-seq_qc_inputs.config
        fi
        if ~{approx_raw_subread_array_lengths_flag} ; then
            echo "~{approx_raw_subread_array_lengths}" >> mas-seq_qc_inputs.config
        else
            echo "NON-EXISTENT-PLACEHOLDER" >> mas-seq_qc_inputs.config
        fi

        if ~{ten_x_metrics_file_flag} ; then
            echo "~{ten_x_metrics_file}" >> mas-seq_qc_inputs.config
        else
            echo "NON-EXISTENT-PLACEHOLDER" >> mas-seq_qc_inputs.config
        fi
        echo "~{mas_seq_model}" >> mas-seq_qc_inputs.config

        echo "~{workflow_dot_file}" >> mas-seq_qc_inputs.config

        # Do the conversion:

        # Run the notebook and populate the notebook itself:
        jupyter nbconvert --execute ~{nb_name} --to notebook --inplace --no-prompt --no-input --clear-output --debug --ExecutePreprocessor.timeout=None

        # Convert the notebook output we created just above here to the HTML report:
        jupyter nbconvert ~{nb_name} --to html --no-prompt --no-input --debug --ExecutePreprocessor.timeout=None

        # Create a tar.gz of the figures directory:
        tar -zcf figures.tar.gz figures

        # Create a dummy pickle for process safety:
        touch dummy.pickle

        # Stop the memory daemon softly.  Then stop it hard if it's not cooperating:
        set +e
        DO_MEMORY_LOG=false
        sleep $(($MEM_LOG_INTERVAL_s  * 2))
        kill -0 $mem_pid &> /dev/null
        if [ $? -ne 0 ] ; then
            kill -9 $mem_pid
        fi
    >>>

    output {
        File populated_notebook = nb_name
        File html_report = html_out
        File figures_tar_gz = "figures.tar.gz"
        File generated_config = "mas-seq_qc_inputs.config"

        Array[File] pickles = glob("*.pickle")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             64,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-jupyter_interactive:0.0.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"  # LOCAL here is a local SSD - much faster and even money with normal disk if preemptible
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
