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

progname=`basename $0`
function usage () {
    cat >&2 <<EOF
USAGE: $progname [options]
Create Drop-seq reference metadata bundle

-n <name>           : Name for reference metadata set to be created.  Required.
-r <referencefasta> : Reference fasta of the Drop-seq reference metadata bundle.  Required.
-s <species>        : Species.  Required.
-g <gtf>            : Gene annotation file.  Required.
-f <filtered-gene-biotype>: Annotations with the given gene_biotype will be filtered. Multiple values may be
                      specified by using this argument more than once, and/or by providing a comma-separated list.
                      Use ValidateReference command to see the gene_biotypes in your GTF in order to decide what to
                      exclude.  Default: not gene biotypes are filtered.
-d <dropseq_root>   : Directory containing Drop-seq executables.  Default: directory containing this script.
-o <outputdir>      : Where to write output bam.  Default: current directory.
-t <tmpdir>         : Where to write temporary files.  Default: Value of $$TMPDIR environment variable.
-a <STAR_path>      : Full path of STAR.  Default: STAR is found via PATH environment variable.
-b <bgzip_path>     : Full path of bgzip: Default: bgzip is found via PATH environment variable.
-i <samtools_path>  : Full path of samtools.  Default: samtools is found via PATH environment variable.
-v                  : verbose
-e                  : merely echo commands instead of executing
EOF
}

function error_exit() {
    echo "ERROR: $@
    " >&2
    exit 1
}


unset reference_name
unset reference_fasta
unset species
unset gtf
filtered_gene_biotypes=
unset dropseq_root
outdir=.
tmpdir=$TMPDIR
star_executable=`which STAR 2> /dev/null`
samtools_executable=`which samtools 2> /dev/null`
bgzip_executable=`which bgzip 2> /dev/null`
verbose=0
ECHO=

set -e
# Fail if any of the commands in a pipeline fails
set -o pipefail



while getopts ":n:r:s:g:f:d:o:t:a:i:b:ev" options; do
  case $options in
    n ) reference_name=$OPTARG;;
    r ) reference_fasta=$OPTARG;;
    s ) species=$OPTARG;;
    g ) gtf=$OPTARG;;
    f ) filtered_gene_biotypes="$filtered_gene_biotypes G=$OPTARG";;
    d ) dropseq_root=$OPTARG;;
    o ) outdir=$OPTARG;;
    t ) tmpdir=$OPTARG;;
    a ) star_executable=$OPTARG;;
    i ) samtools_executable=$OPTARG;;
    b ) bgzip_executable=$OPTARG;;
    e ) ECHO=echo;;
    v ) verbose=1;;
    h ) usage
          exit 1;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;
  esac
done
shift $(($OPTIND - 1))

: ${reference_name:?"ERROR: -n is required"}
: ${reference_fasta:?"ERROR: -r is required"}
: ${species:?"ERROR: -s is required"}
: ${gtf:?"ERROR: -g is required"}
: ${dropseq_root:?"ERROR: -d is required"}
: ${star_executable:?"ERROR: -s is required if STAR is not on PATH"}
: ${samtools_executable:?"ERROR: -i is required if samtools is not on PATH"}
: ${bgzip_executable:?"ERROR: -b is required if bgzip is not on PATH"}

function check_invoke() {
    if (( $verbose ))
    then echo $@
    fi
    if $ECHO $@
    then :
    else error_exit "non-zero exit status " $? " executing $@"
    fi
}

function invoke_picard() {
    check_invoke java -jar $dropseq_root/3rdParty/picard/picard.jar $@
}

function invoke_dropseq() {
    dropseq_program=$1
    shift
    check_invoke $dropseq_root/$dropseq_program $@
}

output_fasta=$outdir/$reference_name.fasta
sequence_dictionary=$outdir/$reference_name.dict
output_gtf=$outdir/$reference_name.gtf
reduced_gtf=$outdir/$reference_name.reduced.gtf

invoke_picard NormalizeFasta INPUT=$reference_fasta OUTPUT=$output_fasta
if [ -e $sequence_dictionary ]
then $ECHO rm $sequence_dictionary
fi
invoke_picard CreateSequenceDictionary REFERENCE=$output_fasta OUTPUT=$sequence_dictionary SPECIES=$species
invoke_dropseq FilterGtf GTF=$gtf SEQUENCE_DICTIONARY=$sequence_dictionary OUTPUT=$output_gtf $filtered_gene_biotypes
invoke_dropseq ConvertToRefFlat ANNOTATIONS_FILE=$output_gtf SEQUENCE_DICTIONARY=$sequence_dictionary OUTPUT=$outdir/$reference_name.refFlat
invoke_dropseq ReduceGtf GTF=$output_gtf SEQUENCE_DICTIONARY=$sequence_dictionary OUTPUT=$reduced_gtf
invoke_dropseq CreateIntervalsFiles SEQUENCE_DICTIONARY=$sequence_dictionary REDUCED_GTF=$reduced_gtf PREFIX=$reference_name \
           OUTPUT=$outdir
$ECHO mkdir -p $outdir/STAR
check_invoke $star_executable --runMode genomeGenerate --genomeDir $outdir/STAR --genomeFastaFiles $output_fasta \
           --sjdbGTFfile $output_gtf --sjdbOverhang 100
check_invoke $bgzip_executable $output_fasta
check_invoke $samtools_executable faidx $output_fasta.gz

echo "Reference metadata created sucessfully"

