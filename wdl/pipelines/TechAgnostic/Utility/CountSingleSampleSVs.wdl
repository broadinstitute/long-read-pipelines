version 1.0

import "../../../tasks/QC/SVvcfMetrics.wdl"

workflow CountSingleSampleSVs {
    meta {
        desciption:
        "Count the number of SVs on a single sample."
        warn:
        "SV VCFs are such a wild-west that there's no guarantee this will work on all VCFs. Please report if you have issue with a VCF."
    }

    parameter_meta {
        sv_vcfs: "gzipped single sample SV VCFs for a sample (you may have multiple VCF produced with different callers, or even different versions of the sample caller)."
        sv_callers: "1-word name of SV callers that generated the VCFs; must be in the same order as sv_vcfs; used as keys in output sv_counts_by_caller_then_type"

        size_lower_bound: "SVs with SVLEN annotation in the INFO column whose absolute value lower than this threshold will be dropped from reports; don't privde negative values"

        sv_counts_by_caller_then_type: "Counts of SVs, grouped by caller, then by sv type"
    }

    input {
        Array[File] sv_vcfs
        Array[String] sv_callers
        Int size_lower_bound = 50
    }

    output {
        Array[Pair[String, Map[String, Int]]] sv_counts_by_caller_then_type = zip(caller, GatherCounts.counts_by_type)
    }

    scatter(pair in zip(sv_callers, sv_vcfs)) {
        String caller = pair.left
        File svvcf = pair.right
        String filter_pass_expression = if (caller == "pav" || caller == "PAV") then "." else "PASS"
        call SVvcfMetrics.GatherCounts { input:
            vcf = svvcf,
            size_lower_bound = size_lower_bound,
            filter_pass_expression = filter_pass_expression
        }
    }
}
