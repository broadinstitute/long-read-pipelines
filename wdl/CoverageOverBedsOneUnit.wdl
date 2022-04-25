version 1.0

import "tasks/Finalize.wdl" as FF
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/SampleLevelAlignedMetrics.wdl" as worker

workflow CoverageOverBedsOneUnit {

    meta {
        description: "Computes coverage over three BED files, for one unit (flowcell, sample)"
    }
    input {
        File bam
        File bai
        String participant_name

        File fiveK_bed
        File threeH_bed
        File exon_bed

        Int mapQ_filter

        String workflow_root_name
        String gcs_out_root_dir
    }

    String dir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_root_name}/~{participant_name}/alignments"

    ##########
    call AM.MosDepthOverBed as fiveK {
        input:
            bam = bam,
            bai = bai,
            bed = fiveK_bed,
            mapQ_filter = mapQ_filter
    }
    call worker.SummarizeDepthOverWholeBed as fiveK_mean {
        input:
            mosdepth_output = fiveK.regions
    }

    call FF.FinalizeToFile as fiveK_final { input: outdir = dir, file = fiveK_mean.cov_summary }

    ##########
    call AM.MosDepthOverBed as threeH {
        input:
            bam = bam,
            bai = bai,
            bed = threeH_bed,
            mapQ_filter = mapQ_filter
    }
    call worker.SummarizeDepthOverWholeBed as threeH_mean {
        input:
            mosdepth_output = threeH.regions
    }

    call FF.FinalizeToFile as threeH_final { input: outdir = dir, file = threeH_mean.cov_summary }

    ##########
    call AM.MosDepthOverBed as exon {
        input:
            bam = bam,
            bai = bai,
            bed = exon_bed,
            mapQ_filter = mapQ_filter
    }
    call worker.SummarizeDepthOverWholeBed as exon_mean {
        input:
            mosdepth_output = exon.regions
    }

    call FF.FinalizeToFile as exon_final { input: outdir = dir, file = exon_mean.cov_summary }

    output {
        File five_K_genes_cov = fiveK_final.gcs_path
        File three_H_genes_cov = threeH_final.gcs_path
        File exons_cov = exon_final.gcs_path
    }
}