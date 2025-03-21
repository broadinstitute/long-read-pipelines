version 1.0

import "../../../tasks/QC/QCAssessment.wdl" as QCAssessment
import "../../../tasks/QC/FastQC.wdl" as FastQC
import "../../../tasks/QC/AlignedMetrics.wdl" as AM

workflow QCAssessment {
    meta {
        description: "Assess quality metrics from mosdepth coverage and callable loci data to determine pass/fail status.  To pass, the fraction of callable bases must be greater than min_callable_fraction and at least min_coverage_threshold_regions fraction of the genome must have a coverage depth greater than min_coverage."
    }

    parameter_meta {
        bam: "Bam file to assess"
        bai: "Index file for given bam file"
        coverage_bed_file: "Bed file for regions to assess"
        min_coverage: "Minimum required mean coverage depth (default: 5)"
        min_coverage_threshold_regions: "Minimum required fraction of genome that must have a coverage depth greater than min_coverage (default: 0.2)"
        min_callable_fraction: "Minimum required fraction of genome that is callable (default: 0.50)"
        prefix: "Prefix for output files"
    }

    input {
        File bam
        File bai
        File coverage_bed_file

        String prefix

        Float min_coverage = 5
        Float min_coverage_threshold_regions = 0.2
        Float min_callable_fraction = 0.50

        RuntimeAttr? runtime_attr_override
    }

    call AM.MosDepthOverBed as t_001_MosDepthOverBed {
            input:
                bam = bam,
                bai = bai,
                bed = select_first(coverage_bed_file)
        }

    call QCAssessment.AssessQualityMetrics as t_002_AssessQualityMetrics {
        input:
            callable_loci_summary = t_001_MosDepthOverBed.callable_loci_summary,
            mosdepth_region_bed = t_001_MosDepthOverBed.regions,
            min_coverage = min_coverage,
            min_coverage_threshold_regions = min_coverage_threshold_regions,
            min_callable_fraction = min_callable_fraction,
            prefix = SM
    }

    call FastQC.FastQC as t_003_FastQC { 
        input: 
            bam = bam, 
            bai = bai 
    }

    output {
 
        # Primary outputs:
        String qc_status = t_002_AssessQualityMetrics.qc_status
        String qc_message = t_002_AssessQualityMetrics.qc_message

        # Mosdepth outputs:
        File mosdepth_global_dist = t_001_MosDepthOverBed.global_dist
        File mosdepth_region_dist = t_001_MosDepthOverBed.region_dist
        File mosdepth_regions = t_001_MosDepthOverBed.regions
        File mosdepth_regions_csi = t_001_MosDepthOverBed.regions_csi

        # Callable loci outputs:
        File callable_loci_summary = select_first([t_001_MosDepthOverBed.callable_loci_summary, t_002_AssessQualityMetrics.callable_loci_summary])
        File callable_loci_bed = select_first([t_001_MosDepthOverBed.callable_loci_bed, t_002_AssessQualityMetrics.callable_loci_bed])

        # FastQC outputs:
        File fastqc_zip = t_003_FastQC.fastqc_zip
        File fastqc_html = t_003_FastQC.fastqc_html
    }
}
