version 1.0

import "../../structs/Structs.wdl"
task RunReportScript {

    meta {
        description: "Use RunReportScript to start the report generation process."
    }

    parameter_meta {
        runtime_attr_override: "Override the default runtime attributes"

        # ------ Summary Page ------ #

        # Sample Info
        sample_name: "name of sequenced sample"
        upload_date: "date sample was sequenced and uploaded"
        species: "species of sample"
        aligned_coverage: "number of reads uniquely mapped to a reference"
        aligned_read_length_n50: "number at which 50%\ of the read lengths are longer than this value" # check
        aligned_read_length_median: "median read length"
        read_qual_median: "median measure of the uncertainty of base calls"

        # Drug Resistance
        drug_resistance_text: "text file used for determining and displaying drug resistances"
        HRP2: "value denoting whether the HRP2 marker is present or not -- true or false"
        HRP3: "value denoting whether the HRP3 marker is present or not -- true or false"

        # Map
        longitude: "longitude value of where the sample was taken"
        latitude: "latitude value of where the sample was taken"
        location: "location where the sample was taken"

        # QC Status
        qc_status: "status to determine whether or not the sequencing run passes quality control standards"

        # ------ Analysis Page ------ #
        # Active Channels
        active_channels: "number of channels active in the sequencing device"

        # Q-Scores Plot
        num_reads_q5: "the number of reads where the probability of a given base call being wrong is approximately 1 in 3"
        num_reads_q7: "the number of reads where the probability of a given base call being wrong is approximately 1 in 5"
        num_reads_q10: "the number of reads where the probability of a given base call being wrong is 1 in 10"
        num_reads_q12: "the number of reads where the probability of a given base call being wrong is approximately 1 in 16"

        # Sequencing Summary
        sample_prep: "type of preparation used for the sample"
        analysis_success: "whether the analysis process completed successfully"
        aligned_bases: "total number of bases aligned to the reference genome"
        aligned_reads: "total number of reads aligned to the reference genome"
        fraction_aligned_bases: "number of bases aligned out of all bases in the sequence"
        # average_identity:

        # Coverage Plot -- incomplete

    }

    input {
        RuntimeAttr? runtime_attr_override

        # ------ Summary Page ------ #

        # Sample Info
        String sample_name
        String upload_date
        String species
        Float aligned_coverage
        Float aligned_read_length_n50
        Float aligned_read_length_median
        Float read_qual_median

        # Drug Resistance
        File? drug_resistance_text
        String? HRP2
        String? HRP3

        # Map
        Float? longitude
        Float? latitude
        String? location
        
        # QC Status
        String qc_status

        # ------ Analysis Page ------ #
        
        # Active Channels
        Int active_channels
        
        # Q-Scores Plot
        Int num_reads_q5
        Int num_reads_q7
        Int num_reads_q10
        Int num_reads_q12
        Int num_reads_q15

        # Sequencing Summary
        String sample_prep
        String analysis_success
        Float aligned_bases
        Int aligned_reads
        Float fraction_aligned_bases
        Float average_identity

        # Coverage Plot -- incomplete  
    }

    Int disk_size_gb = 20 + ceil(size(drug_resistance_text, "GB"))

    # incomplete
    command <<<
        set -euxo
        pwd
        ls
        echo "BEGIN COMMAND"

        python3 /report-files/report_gen.py \
            --sample_name ~{sample_name} \
            --upload_date ~{upload_date} \
            --species ~{species} \
            --aligned_coverage ~{aligned_coverage} \
            --aligned_read_length_n50 ~{aligned_read_length_n50} \
            --aligned_read_length_median ~{aligned_read_length_median} \
            --read_qual_median ~{read_qual_median} \
            --drug_resistance_text ~{drug_resistance_text} \
            --HRP2 ~{HRP2} \
            --HRP3 ~{HRP3} \
            --longitude ~{longitude} \
            --latitude ~{latitude} \
            --location ~{location} \
            --qc_status ~{qc_status} \
            --active_channels ~{active_channels} \
            --num_reads_q5 ~{num_reads_q5} \
            --num_reads_q7 ~{num_reads_q7} \
            --num_reads_q10 ~{num_reads_q10} \
            --num_reads_q12 ~{num_reads_q12} \
            --num_reads_q15 ~{num_reads_q15} \
            --aligned_bases ~{aligned_bases} \
            --sample_prep ~{sample_prep} \
            --analysis_success ~{analysis_success} \
            --aligned_reads ~{aligned_reads} \
            --fraction_aligned_bases ~{fraction_aligned_bases} \
            --average_identity ~{average_identity}
    >>>

    output {
        File report = "~{sample_name}_lrma_report.html"
    }
    

    # ----------------------------------------------------------- #
    # Runtime Config

    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             16,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-malaria-reports:0.0.1"
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
