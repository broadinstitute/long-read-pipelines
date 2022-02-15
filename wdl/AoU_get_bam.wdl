version 1.0


##########################################################################################
## A workflow
##########################################################################################



import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils

workflow MitochondriaProcessing{
    input{
        File bam
        File bai
        String locus
        String prefix = "subset"
    }

    parameter_meta{
        bam:      "GCS path to raw subread bam"
        bai:      "index for bam file"
        locus:    "genomic locus to select"
        prefix:   "prefix for output bam and bai file names"
    }
#    String outdir =

    # subset chrM reads from bam
    call Utils.SubsetBam as SubsetBam {input: bam = bam, bai = bai} #task_var = workflow_var

    output{File chrM_bam = SubsetBam.subset_bam
           File chrM_bam_bai = SubsetBam.subset_bai}
}
















