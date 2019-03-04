// Bpipe pipeline config file
// Paths are relative to the directory the pipeline is running in, so absolute
// paths are recommended.

// Adjust parameters
PLATFORM='illumina'

// Number of threads to use for BWA
threads=2

// For exome pipeline only ***Edit before running the exome pipeline***
EXOME_TARGET=""

// STRETCH installation location
STRETCH="/STRetch"

// Paths to tools used by the pipeline
bpipe="$STRETCH/tools/bin/bpipe"
python="$STRETCH/tools/bin/python"
goleft="$STRETCH/tools/bin/goleft"
bedtools="$STRETCH/tools/bin/bedtools"
bwa="$STRETCH/tools/bin/bwa"
samtools="$STRETCH/tools/bin/samtools"

bazam="$STRETCH/tools/bin/bazam.jar"
picard="$STRETCH/tools/bin/picard.jar"

// Path to reference data
//refdir="$STRETCH/reference-data"
// Decoy reference assumed to have matching .genome file in the same directory
//REF="$refdir/full_genome_and_STRdecoys.sorted.fasta"
//DECOY_BED="$refdir/STRdecoys.sorted.bed"
//STR_BED="$refdir/GRCh37/repeats/hg19.simpleRepeat_period1-6_dedup.sorted.bed"


// By default, uses other samples in the same batch as a control
//CONTROL=""
// Uncomment the line below to use a set of WGS samples as controls, or specify your own
CONTROL="/PCRfreeWGS.controls.tsv"

//R="/usr/local/installed/R/3.3.2/lib64/R/bin/R"
