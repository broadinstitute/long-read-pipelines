version 1.0

import "../PBUtils.wdl" as PB
import "../Utils.wdl"

workflow CollectSMRTCellUnalignedMetrics {
    input {
        File smrtcell_pbi
    }

    call PB.SummarizePBI { input: pbi = smrtcell_pbi, runtime_attr_override = { 'mem_gb': 72 } }
    call Utils.MapToTsv {input: my_map = SummarizePBI.results, name_of_file = "pbi.summary.tsv" }

    output {

        File pbi_summary = MapToTsv.result

        Float polymerase_read_length_mean = SummarizePBI.results['polymerase_mean']
        Float polymerase_read_length_N50 = SummarizePBI.results['polymerase_n50']
        Float subread_read_length_mean = SummarizePBI.results['subread_mean']
        Float subread_read_length_N50 = SummarizePBI.results['subread_n50']
        Float num_reads = SummarizePBI.results['reads']
        Float num_bases = SummarizePBI.results['bases']
        Float read_length_mean = SummarizePBI.results['subread_mean']
        Float read_length_median = SummarizePBI.results['subread_median']
        Float read_length_stdev = SummarizePBI.results['subread_stdev']
        Float read_length_N50 = SummarizePBI.results['subread_n50']
        Float read_qual_mean = SummarizePBI.results['mean_qual']
        Float read_qual_median = SummarizePBI.results['median_qual']
    }
}
