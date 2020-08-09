
#gtf=/ahg/regevdata/projects/612-eqtl/eQTL_testing/LongReadData/Code/ref.flip.short.gtf
gtf=/ahg/regevdata/projects/612-eqtl/eQTL_testing/Code/fixed.ref.flip.gtf
bamfile=$1
#outfile=${bamfile}.GE.bam
outfile=$2

TagReadWithGeneFunction I=$bamfile O=$outfile ANNOTATIONS_FILE=$gtf USE_STRAND_INFO=false
