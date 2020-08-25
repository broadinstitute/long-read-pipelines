#!/usr/bin/env bash
# MIT License
#
# Copyright 2018 Broad Institute
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

tmpdir=
outdir=`pwd`
genomedir=
reference=
echo_prefix=
dropseq_root=$(dirname $0)
star_executable=STAR
keep_intermediates=0
progname=`basename $0`

function usage () {
    cat >&2 <<EOF
USAGE: $progname [options] <unmapped-queryname-sorted.bam>
Perform Drop-seq tagging, trimming and alignment

-g <genomedir>      : Directory of STAR genome directory.  Required.
-r <referencefasta> : Reference fasta of the Drop-seq reference metadata bundle.  Required.
-d <dropseq_root>   : Directory containing Drop-seq executables.  Default: directory containing this script.
-o <outputdir>      : Where to write output bam.  Default: current directory.
-t <tmpdir>         : Where to write temporary files.  Default: a new subdirectory in $TMPDIR.
-s <STAR_path>      : Full path of STAR.  Default: STAR is found via PATH environment variable.
-e                  : Echo commands instead of executing them.
-k                  : Keep intermediate files
EOF
}

function error_exit() {
    echo "ERROR: $1
    " >&2
    usage
    exit 1
}

function check_set() {
    value=$1
    name=$2
    flag=$3

    if [[ -z "$value" ]]
    then error_exit "$name has not been specified.  $flag flag is required"
    fi
}

set -e
# Fail if any of the commands in a pipeline fails
set -o pipefail

while getopts ":d:t:o:g:r:es:k" options; do
  case $options in
    d ) dropseq_root=$OPTARG;;
    t ) tmpdir=$OPTARG;;
    o ) outdir=$OPTARG;;
    g ) genomedir=$OPTARG;;
    r ) reference=$OPTARG;;
    s ) star_executable=$OPTARG;;
    e ) echo_prefix="echo";;
    k ) keep_intermediates=1;;
    h ) usage
          exit 1;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;
  esac
done
shift $(($OPTIND - 1))

check_set "$dropseq_root" "Drop-seq root" "-d"
check_set "$genomedir" "Genome directory" "-g"
check_set "$reference" "Reference fasta"  "-r"

if (( $# != 1 ))
then error_exit "Incorrect number of arguments"
fi

if [[ -z "$tmpdir" ]]
then tmpdir=`mktemp -d`
     echo "Using temporary directory $tmpdir"
fi

if [[ "$star_executable" != "STAR" ]]
then if [[ ! ( -x $star_executable && -f $star_executable ) ]]
     then error_exit "STAR executable $star_executable passed via -s does not exist or is not executable"
     fi
elif which STAR > /dev/null
then echo > /dev/null
else error_exit "STAR executable must be on the path"
fi

reference_basename=$(basename $(basename $reference .gz) .fasta)
gene_intervals=$(dirname $reference)/$reference_basename.genes.intervals
exon_intervals=$(dirname $reference)/$reference_basename.exon.intervals
refflat=$(dirname $reference)/$reference_basename.refFlat
picard_jar=${dropseq_root}/3rdParty/picard/picard.jar

unmapped_bam=$1
tagged_unmapped_bam=${tmpdir}/unaligned_mc_tagged_polyA_filtered.bam
aligned_sam=${tmpdir}/star.Aligned.out.sam
aligned_sorted_bam=${tmpdir}/aligned.sorted.bam
files_to_delete=


# Stage 4
merge_bam=""
tag_with_gene_exon=""

# Stage 1: pre-alignment tag and trim

# cellular tag
$echo_prefix ${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Cellular.bam_summary.txt \
  BASE_RANGE=1-12 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY=1 \
  INPUT=${unmapped_bam} OUTPUT=$tmpdir/unaligned_tagged_Cell.bam
files_to_delete="$files_to_delete $tmpdir/unaligned_tagged_Cell.bam"

# molecular tag
$echo_prefix ${dropseq_root}/TagBamWithReadSequenceExtended SUMMARY=${outdir}/unaligned_tagged_Molecular.bam_summary.txt \
  BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=true TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1 \
  INPUT=$tmpdir/unaligned_tagged_Cell.bam OUTPUT=$tmpdir/unaligned_tagged_CellMolecular.bam
files_to_delete="$files_to_delete $tmpdir/unaligned_tagged_CellMolecular.bam"

# quality filter
$echo_prefix ${dropseq_root}/FilterBam TAG_REJECT=XQ INPUT=$tmpdir/unaligned_tagged_CellMolecular.bam \
  OUTPUT=$tmpdir/unaligned_tagged_filtered.bam
files_to_delete="$files_to_delete $tmpdir/unaligned_tagged_filtered.bam"

# read trimming
$echo_prefix ${dropseq_root}/TrimStartingSequence OUTPUT_SUMMARY=${outdir}/adapter_trimming_report.txt \
  SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5 INPUT=$tmpdir/unaligned_tagged_filtered.bam \
  OUTPUT=$tmpdir/unaligned_tagged_trimmed_smart.bam
files_to_delete="$files_to_delete $tmpdir/unaligned_tagged_trimmed_smart.bam"

$echo_prefix ${dropseq_root}/PolyATrimmer OUTPUT=${tagged_unmapped_bam} OUTPUT_SUMMARY=${outdir}/polyA_trimming_report.txt \
  MISMATCHES=0 NUM_BASES=6 NEW=true INPUT=$tmpdir/unaligned_tagged_trimmed_smart.bam
files_to_delete="$files_to_delete ${tagged_unmapped_bam}"


# Stage 2: alignment
$echo_prefix java -Xmx500m -jar ${picard_jar} SamToFastq INPUT=${tmpdir}/unaligned_mc_tagged_polyA_filtered.bam \
  FASTQ=$tmpdir/unaligned_mc_tagged_polyA_filtered.fastq
files_to_delete="$files_to_delete $tmpdir/unaligned_mc_tagged_polyA_filtered.fastq"

$echo_prefix $star_executable --genomeDir ${genomedir} --outFileNamePrefix ${tmpdir}/star. \
  --readFilesIn $tmpdir/unaligned_mc_tagged_polyA_filtered.fastq
files_to_delete="$files_to_delete ${aligned_sam}"

# Stage 3: sort aligned reads (STAR does not necessarily emit reads in the same order as the input)
$echo_prefix java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx4000m -jar ${picard_jar} \
  SortSam INPUT=${aligned_sam} OUTPUT=${aligned_sorted_bam} SORT_ORDER=queryname TMP_DIR=${tmpdir}
files_to_delete="$files_to_delete ${aligned_sorted_bam}"

# Stage 4: merge and tag aligned reads
$echo_prefix java -Xmx4000m -jar ${picard_jar} MergeBamAlignment REFERENCE_SEQUENCE=${reference} UNMAPPED_BAM=${tagged_unmapped_bam} \
  ALIGNED_BAM=${aligned_sorted_bam} INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false CLIP_ADAPTERS=false \
  TMP_DIR=${tmpdir} OUTPUT=$tmpdir/merged.bam
files_to_delete="$files_to_delete $tmpdir/merged.bam"

$echo_prefix ${dropseq_root}/TagReadWithInterval I=$tmpdir/merged.bam O=$tmpdir/gene_tagged.bam TMP_DIR=${tmpdir} \
  INTERVALS=${gene_intervals} TAG=XG
files_to_delete="$files_to_delete $tmpdir/gene_tagged.bam"

$echo_prefix ${dropseq_root}/TagReadWithGeneFunction O=${tmpdir}/function_tagged.bam ANNOTATIONS_FILE=${refflat} \
  INPUT=$tmpdir/gene_tagged.bam
files_to_delete="$files_to_delete $tmpdir/function_tagged.bam"

# Stage 5: bead repair
$echo_prefix ${dropseq_root}/DetectBeadSubstitutionErrors INPUT=${tmpdir}/function_tagged.bam OUTPUT=${tmpdir}/substitution_repaired.bam \
  TMP_DIR=$tmpdir MIN_UMIS_PER_CELL=20 OUTPUT_REPORT=${outdir}/substitution_error_report.txt
files_to_delete="$files_to_delete ${tmpdir}/substitution_repaired.bam"

$echo_prefix ${dropseq_root}/DetectBeadSynthesisErrors INPUT=${tmpdir}/substitution_repaired.bam MIN_UMIS_PER_CELL=20 \
  OUTPUT_STATS=${outdir}/synthesis_error_stats.txt SUMMARY=${outdir}/synthesis_error_summary.txt \
  REPORT=${outdir}/synthesis_error_report.txt CREATE_INDEX=true TMP_DIR=$tmpdir OUTPUT=$outdir/final.bam

if (( $keep_intermediates == 0 ))
then $echo_prefix rm $files_to_delete
fi

echo "Completed successfully."

