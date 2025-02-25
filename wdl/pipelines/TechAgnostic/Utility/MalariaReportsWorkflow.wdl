version 1.0
import "../../../tasks/Visualization/MalariaReports.wdl" as MRS
import "../../../tasks/Visualization/MalariaReportsSnapshots.wdl" as Snapshot

workflow GenerateMalariaReports {

    meta {
        author: "Bridget Knight"
        description: "## Report Generation \n This workflow calls the Python script that generates a sequencing report."
    }

    parameter_meta {
        # ------ Summary Page ------ #

        # Sample Info
        sample_name: "name of sequenced sample"
        upload_date: "date sample was uploaded"
        collection_date: "date sample was collected"
        sequencing_date: "date sample was sequenced"
        species: "species of sample"
        aligned_coverage: "number of reads uniquely mapped to a reference"
        aligned_read_length: "number at which 50% of the read lengths are longer than this value" # check
        pct_properly_paired_reads: "percent of reads that are properly paired/aligned"
        read_qual_mean: "mean measure of the uncertainty of base calls"

        # Drug Resistance
        drug_resistance_text: "text file used for determining and displaying drug resistances"
        HRP2: "value denoting whether the HRP2 marker is present or not -- true or false"
        HRP3: "value denoting whether the HRP3 marker is present or not -- true or false"
        predicted_drug_status_chloroquine: "value denoting the sensitivity of this sample to chloroquine"
        predicted_drug_status_pyrimethamine: "value denoting the sensitivity of this sample to pyrimethamine"
        predicted_drug_status_sulfadoxine: "value denoting the sensitivity of this sample to sulfadoxine"
        predicted_drug_status_mefloquine: "value denoting the sensitivity of this sample to mefloquine"
        predicted_drug_status_artemisinin: "value denoting the sensitivity of this sample to artemisinin"
        predicted_drug_status_piperaquine: "value denoting the sensitivity of this sample to piperaquine"

        # Map
        longitude: "longitude value of where the sample was taken"
        latitude: "latitude value of where the sample was taken"
        location: "location where the sample was taken"
        code: "the three letter code denoting the corresponding health facility/site"
        location_table: "table containing the region, district, site, and facility type for this sample's location"

        # QC Status
        qc_pass: "status to determine whether or not the sequencing run passes quality control standards"

        # ------ Analysis Page ------ #
        barcode: "the barcode of the sample"

        # Q-Scores Plot
        num_reads_q5: "the number of reads where the probability of a given base call being wrong is approximately 1 in 3"
        num_reads_q7: "the number of reads where the probability of a given base call being wrong is approximately 1 in 5"
        num_reads_q10: "the number of reads where the probability of a given base call being wrong is 1 in 10"
        num_reads_q12: "the number of reads where the probability of a given base call being wrong is approximately 1 in 16"
        num_reads_q15: "the number of reads where the probability of a given base call being wrong is approximately 1 in 31"

        # Sequencing Summary
        aligned_bases: "total number of bases aligned to the reference genome"
        aligned_reads: "total number of reads aligned to the reference genome"
        fraction_aligned_bases: "number of bases aligned out of all bases in the sequence"

        # Coverage Plot
        fastqc_path: "directory of fastqc_report used for finding BAM files"
        coverage_bin_size: "number to use as size of bins for coverage plot generation; default is 1500"

        # ------ IGV Snapshots ------ #
        regions_bed:        "GCS path to bed file containing drug resistance loci"
        ref_map_file:            "table indicating reference sequence and auxillary file locations"
        aligned_bam:        "GCS path to aligned bam"
        aligned_bai:        "GCS path to aligned bai"
        gcs_out_root_dir:   "GCS bucket to store the results"

        ont_qc_report: "ONT QC report file"

        tech_flag: "string denoting if the sample was sequenced using long read or short read techniques"
    }

    input {
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
        String? predicted_drug_status_chloroquine
        String? predicted_drug_status_pyrimethamine
        String? predicted_drug_status_sulfadoxine
        String? predicted_drug_status_mefloquine
        String? predicted_drug_status_artemisinin
        String? predicted_drug_status_piperaquine

        # Map
        Float? longitude
        Float? latitude
        String? location
        String? code
        File location_table
        
        # QC Pass
        String qc_pass

        # ------ Analysis Page ------ #
        
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
        File? fastqc_path
        Int? coverage_bin_size  

        # ------ IGV Snapshots ------ #
        File aligned_bam
        File aligned_bai
        File regions_bed
        File ref_map_file

        File? ont_qc_report

        String tech_flag
    }
    
    Map[String, String] ref_map = read_map(ref_map_file)
    String fasta = ref_map["fasta"]
    String fai = ref_map["fai"]

    call Snapshot.DrugResIGV as GenerateSnapshots {
        input:
            sample_name = sample_name,
            fasta_path = fasta,
            fasta_index_path = fai,
            regions_bed = regions_bed,
            aligned_bam = aligned_bam,
            aligned_bai = aligned_bai
    }

    call MRS.RunReportScript as RunReportScript { 
        input: 
            sample_name = sample_name,
            upload_date = upload_date,
            collection_date = collection_date,
            sequencing_date = sequencing_date,
            species = species,
            aligned_coverage = aligned_coverage,
            aligned_read_length = aligned_read_length,
            pct_properly_paired_reads = pct_properly_paired_reads,
            read_qual_mean = read_qual_mean,
            drug_resistance_text = drug_resistance_text,
            HRP2 = HRP2,
            HRP3 = HRP3,
            predicted_drug_status_chloroquine=predicted_drug_status_chloroquine,
            predicted_drug_status_pyrimethamine=predicted_drug_status_pyrimethamine,
            predicted_drug_status_sulfadoxine=predicted_drug_status_sulfadoxine,
            predicted_drug_status_mefloquine=predicted_drug_status_mefloquine,
            predicted_drug_status_artemisinin=predicted_drug_status_artemisinin,
            predicted_drug_status_piperaquine=predicted_drug_status_piperaquine,
            longitude = longitude,
            latitude = latitude,
            location = location,
            barcode = barcode,
            num_reads_q5 = num_reads_q5,
            num_reads_q7 = num_reads_q7,
            num_reads_q10 = num_reads_q10,
            num_reads_q12 = num_reads_q12,
            num_reads_q15 = num_reads_q15,
            aligned_bases = aligned_bases,
            aligned_reads = aligned_reads,
            fraction_aligned_bases = fraction_aligned_bases,
            average_identity = average_identity,
            fastqc_path = fastqc_path,
            ont_qc_report = ont_qc_report,
            coverage_bin_size = coverage_bin_size,
            snapshots = GenerateSnapshots.snapshots,
            regions_bed = regions_bed,
            qc_pass = qc_pass,
            code = code,
            location_table = location_table,
            tech_flag = tech_flag
    }

    output {
        File report = RunReportScript.report
    }
}