version 1.0

import "tasks/Finalize.wdl" as FF
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/SampleLevelAlignedMetrics.wdl" as worker

workflow CoverageOverBed {
    input {
        File bam
        File bai
        String participant_name

        File fiveK_bed
        File threeH_bed
        File exon_bed

        String workflow_root_name
        String gcs_out_root_dir
    }

    String dir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_root_name}/~{participant_name}/alignments"

    ##########
    call AM.MosDepthOverBed as One {
        input:
            bam = bam,
            bai = bai,
            bed = fiveK_bed
    }
    call worker.SummarizeDepthOverWholeBed as OneS {
        input:
            mosdepth_output = One.regions
    }
    call FF.FinalizeToFile as OneF { input: outdir = dir, file = OneS.cov_summary }

    ##########
    call AM.MosDepthOverBed as Two {
        input:
            bam = bam,
            bai = bai,
            bed = threeH_bed
    }
    call worker.SummarizeDepthOverWholeBed as TwoS {
        input:
            mosdepth_output = Two.regions
    }
    call FF.FinalizeToFile as TwoF { input: outdir = dir, file = TwoS.cov_summary }

    ##########
    call AM.MosDepthOverBed as Three {
        input:
            bam = bam,
            bai = bai,
            bed = exon_bed
    }
    call worker.SummarizeDepthOverWholeBed as ThreeS {
        input:
            mosdepth_output = Three.regions
    }
    call FF.FinalizeToFile as ThreeF { input: outdir = dir, file = ThreeS.cov_summary }

    output {
        File all_5kgenes_cov_summary = OneF.gcs_path
        File med_chall_genes_cov_summary = TwoF.gcs_path
        File exons_cov_summary = ThreeF.gcs_path
    }
}