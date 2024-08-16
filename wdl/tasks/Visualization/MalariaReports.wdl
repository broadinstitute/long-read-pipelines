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
        upload_date: "date sample was uploaded"
        sequencing_date: "date sample was sequenced"
        collection_date: "date sample was collected"
        species: "species of sample"
        aligned_coverage: "number of reads uniquely mapped to a reference"
        aligned_read_length: "number at which 50% of the read lengths are longer than this value" # check
        pct_properly_paired_reads: "median read length"
        read_qual_mean: "mean measure of the uncertainty of base calls"

        # Drug Resistance
        drug_resistance_text: "text file used for determining and displaying drug resistances"
        HRP2: "value denoting whether the HRP2 marker is present or not -- true or false"
        HRP3: "value denoting whether the HRP3 marker is present or not -- true or false"

        # Map
        longitude: "longitude value of where the sample was taken"
        latitude: "latitude value of where the sample was taken"
        location: "location where the sample was taken"
        code: "the three letter code denoting the corresponding health facility/site"
        location_table: "table containing the region, district, site, and facility type for this sample's location"

        # QC Status
        qc_pass: "status to determine whether or not the sequencing run passes quality control standards"

        # ------ Analysis Page ------ #
        # Main Box
        barcode: "the barcode of the sample"

        # Q-Scores Plot
        num_reads_q5: "the number of reads where the probability of a given base call being wrong is approximately 1 in 3"
        num_reads_q7: "the number of reads where the probability of a given base call being wrong is approximately 1 in 5"
        num_reads_q10: "the number of reads where the probability of a given base call being wrong is 1 in 10"
        num_reads_q12: "the number of reads where the probability of a given base call being wrong is approximately 1 in 16"

        # Sequencing Summary
        aligned_bases: "total number of bases aligned to the reference genome"
        aligned_reads: "total number of reads aligned to the reference genome"
        fraction_aligned_bases: "percentage of bases aligned out of all bases in the sequence"

        # Coverage Plot
        # coverage_dir: "directory of BAM files for coverage plot generation"
        fastqc_path: "directory of fastqc_report used for finding BED and BAM files"
        coverage_bin_size: "number to use as size of bins for coverage plot generation; default is 1500"

        ont_qc_report: "ONT QC report file"

        tech_flag: "string denoting if the sample was sequenced using long read or short read techniques"
    }

    input {
        RuntimeAttr? runtime_attr_override

        # ------ Summary Page ------ #

        # Sample Info
        String sample_name
        String upload_date
        String? collection_date
        String? sequencing_date
        String? species
        Float aligned_coverage
        Float? aligned_read_length
        Float? pct_properly_paired_reads
        Float? read_qual_mean

        # Drug Resistance
        File? drug_resistance_text
        String? HRP2
        String? HRP3

        # Map
        Float? longitude
        Float? latitude
        String? location
        String? code
        File location_table
        
        # QC Pass
        String qc_pass

        # ------ Analysis Page ------ #
        
        # Barcode
        String? barcode
        
        # Q-Scores Plot
        Int num_reads_q5
        Int num_reads_q7
        Int num_reads_q10
        Int num_reads_q12
        Int num_reads_q15

        # Sequencing Summary
        Float aligned_bases
        Int aligned_reads
        Float fraction_aligned_bases
        Float average_identity

        # Coverage Plot -- incomplete 
        # String? coverage_dir
        File? fastqc_path
        Int? coverage_bin_size

        # Snapshots
        File regions_bed
        Array[File] snapshots

        # ONT QC Report Page
        File? ont_qc_report

        String tech_flag
    }

    Int disk_size_gb = 20 + ceil(size(drug_resistance_text, "GB"))

    # Wrap location in quotes in case it has spaces
    String wrap_location = "'~{location}'"

    String results_dir = if (defined(fastqc_path)) then sub(select_first([fastqc_path]), "results\/.*", "") else ""
    String coverage_dir = "~{results_dir}results/SRFlowcell/~{sample_name}/metrics/coverage/"
    String coverage_regex = "~{coverage_dir}*?[0-9]_v3.regions.bed.gz"

    String quality_report = select_first([fastqc_path, ont_qc_report, ""])

    command <<<
        set -euxo
        pwd
        ls

        if [[ ~{coverage_regex} == "gs://"* ]]; then
            echo "Retrieving BED files..."
            echo ~{coverage_dir}
            mkdir -p /report-files/data/coverage
            gsutil ls ~{coverage_regex}  > filelist.txt
            echo "COPYING..."
            cat filelist.txt | gsutil -m cp -I /report-files/data/coverage
        else
            echo "No BED files found."
        fi

        if [ ! -z ~{quality_report} ]; then
            echo "Retrieving quality report..."
            mkdir -p /report-files/data/quality_report
            echo ~{quality_report}
            gsutil cp ~{quality_report} /report-files/data/quality_report
        else
            echo "No quality report found."
        fi
        
        #mkdir -p /report-files/data/quality_report
        #echo ~{fastqc_path}
        #gsutil cp ~{fastqc_path} /report-files/data/quality_report
        
        echo "CREATING SUMMARY REPORT..."
        python3 /report-files/report_gen.py \
            --sample_name ~{sample_name} \
            --upload_date ~{upload_date} \
            --collection_date ~{default="Unknown" collection_date} \
            --sequencing_date ~{default="Unknown" sequencing_date} \
            --species ~{default="P. falciparum" species} \
            --aligned_coverage ~{aligned_coverage} \
            --aligned_read_length ~{default=0 aligned_read_length} \
            --pct_properly_paired_reads ~{default=0 pct_properly_paired_reads} \
            --read_qual_mean ~{read_qual_mean} \
            --drug_resistance_text ~{default="None" drug_resistance_text} \
            --HRP2 ~{default="Unknown" HRP2} \
            --HRP3 ~{default="Unknown" HRP3} \
            --longitude ~{default=0 longitude} \
            --latitude ~{default=0 latitude} \
            --location ~{default="Unknown" wrap_location} \
            --barcode ~{default="Unknown" barcode} \
            --num_reads_q5 ~{num_reads_q5} \
            --num_reads_q7 ~{num_reads_q7} \
            --num_reads_q10 ~{num_reads_q10} \
            --num_reads_q12 ~{num_reads_q12} \
            --num_reads_q15 ~{num_reads_q15} \
            --aligned_bases ~{aligned_bases} \
            --aligned_reads ~{aligned_reads} \
            --fraction_aligned_bases ~{fraction_aligned_bases} \
            --average_identity ~{average_identity} \
            --coverage_bin_size ~{default=750 coverage_bin_size} \
            --snapshots ~{sep="," snapshots} \
            --regions_bed ~{regions_bed} \
            --qc_pass ~{default="N/A" qc_pass} \
            --code ~{default="" code} \
            --location_table ~{location_table} \
            --tech_flag ~{tech_flag}
        echo "DONE!"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-malaria-reports:0.0.2"
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
