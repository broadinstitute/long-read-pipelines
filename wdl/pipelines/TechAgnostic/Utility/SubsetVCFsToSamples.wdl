version 1.0

import "../../../structs/Structs.wdl"
import "../../../tasks/Utility/VariantUtils.wdl" as VARUTIL

workflow SubsetVCFsToSamples {
    meta {
        description: "Subset each of a list of (multi-sample) VCFs to a specified set of samples, emitting one multi-sample VCF per input that retains only those samples (no per-sample splitting). Each VCF is streamed from gs:// and processed in its own scatter shard. Exactly one of sample_names or sample_name_list must be provided; every requested sample must exist in each VCF."
        author: "Jonn Smith"
    }

    parameter_meta {
        input_vcfs:        "VCF files to subset. Each is streamed from its gs:// location and subset independently."
        input_vcf_indices: "Index files for input_vcfs, matched by position (same length and order)."
        sample_names:      "Samples to keep, as an inline list of names. Mutually exclusive with sample_name_list; exactly one of the two must be given."
        sample_name_list:  "Samples to keep, as a file with one sample name per line. Mutually exclusive with sample_names; exactly one of the two must be given."
        error_if_sample_missing: "When true (default), a requested sample absent from any VCF fails that shard (all absent samples listed on stderr). When false, absent samples produce a warning and processing continues with the samples that are present."
    }

    input {
        Array[File] input_vcfs
        Array[File] input_vcf_indices

        Array[String]? sample_names
        File? sample_name_list

        Boolean error_if_sample_missing = true
    }

    scatter (idx in range(length(input_vcfs))) {
        call VARUTIL.SubsetVCFToSamples {
            input:
                input_vcf = input_vcfs[idx],
                input_vcf_index = input_vcf_indices[idx],
                sample_names = sample_names,
                sample_name_list = sample_name_list,
                error_if_sample_missing = error_if_sample_missing
        }
    }

    output {
        Array[File] subset_vcfs = SubsetVCFToSamples.subset_vcf
        Array[File] subset_vcf_indices = SubsetVCFToSamples.subset_vcf_index
    }
}
