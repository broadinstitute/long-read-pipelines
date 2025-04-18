version 1.0

import "../../../tasks/QC/QCAssessment.wdl" as QCAssessmentTasks
import "../../../tasks/QC/FastQC.wdl" as FastQC
import "../../../tasks/QC/AlignedMetrics.wdl" as AM

workflow QCAssessmentDepthBased {
    meta {
        description: "Assess quality metrics from mosdepth coverage and callable loci data to determine pass/fail status.  To pass, the fraction of callable bases must be greater than min_callable_fraction and at least min_coverage_threshold_regions fraction of the genome must have a coverage depth greater than min_coverage."
    }

    parameter_meta {
        bam: "Bam file to assess"
        bai: "Index file for given bam file"
        coverage_bed_file: "Bed file for regions to assess"
        ref_map_file:  "Reference map file for the primary reference sequence and auxillary file locations."
        min_coverage: "Minimum required mean coverage depth (default: 5) "
        min_coverage_threshold_regions: "Minimum required fraction of genome that must have a coverage depth greater than min_coverage (default: 0.2)"
        min_callable_fraction: "Minimum required fraction of genome that is callable (default: 0.50)"
        prefix: "Prefix for output files"
    }

    input {
        File bam
        File bai
        File coverage_bed_file

        File ref_map_file

        String prefix

        Float min_coverage = 5
        Float min_coverage_threshold_regions = 0.2
        Float min_callable_fraction = 0.50
    }

    # Get ref info:
    Map[String, String] ref_map = read_map(ref_map_file)

    call AM.CallableLoci as t_001_CallableLoci {
        input:
            bam_file = bam,
            bam_index = bai,
            ref_fasta = ref_map['fasta'],
            ref_fasta_index = ref_map['fai'],
            ref_dict = ref_map['dict'],
            min_depth = 5,
            prefix = prefix
    }

    call AM.MosDepthOverBed as t_001_MosDepthOverBed {
            input:
                bam = bam,
                bai = bai,
                bed = coverage_bed_file
        }

    call QCAssessmentTasks.AssessQualityMetrics as t_002_AssessQualityMetrics {
        input:
            callable_loci_summary = t_001_CallableLoci.callable_loci_summary,
            mosdepth_region_bed = t_001_MosDepthOverBed.regions,
            min_coverage = min_coverage,
            min_coverage_threshold_regions = min_coverage_threshold_regions,
            min_callable_fraction = min_callable_fraction,
            prefix = prefix
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
        File callable_loci_summary = t_001_CallableLoci.callable_loci_summary
        File callable_loci_bed = t_001_CallableLoci.callable_loci_bed

        # FastQC outputs:
        File fastqc_stats = t_003_FastQC.stats
        File fastqc_report = t_003_FastQC.report
    }
}
