#!/bin/bash


#$ -cwd
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=90g
#$ -e ErrFiles/err.ann.err
#$ -o ErrFiles/log.ann.out
#$ -l os=RedHat7
#$ -l h_rt=40:00:00
#$ -t 3-3


source /broad/software/scripts/useuse

use UGER
use .htslib-1.8
use GCC-5.2 
use .samtools-1.8
use BEDTools
use UGER
use .java-jdk-1.8.0_181-x86-64
use R-3.5
use .python-3.6.0 


export R_LIBS_USER=/psych/genetics_data/ssimmons/R/x86_64-pc-linux-gnu-library/3.5_v2

export PATH=$PATH:Dropseq/Drop-seq_tools-2.3.0:/ahg/regevdata/users/libo/software/scan_snps

SGE_TASK_ID=3
#SEEDFILE=/ahg/regevdata/projects/612-eqtl/eQTL_testing/LongReadData/Oct1_new/forPipeline.txt

SEEDFILE=/ahg/regevdata/projects/612-eqtl/eQTL_testing/LongReadData/Annotated/forPipeline.txt


bamfile=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $1}')
prefix=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $4}')
ind=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $5}')

prefix=${prefix}_test

echo $bamfile
echo $vcffile
echo $prefix

#codedir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
codedir=$PWD

vcffile=$codedir/ref/SNPs.vcf


echo First get gene information, ignoring strand
outbam=$prefix.GE.bam

gtf=$codedir/ref/genes.chr.gtf


head $gtf  

TagReadWithGeneFunction I=$bamfile O=$outbam ANNOTATIONS_FILE=$gtf USE_STRAND_INFO=false


bamfile=$outbam
echo " "
echo Nice Make Nice bam, 10X like!
outbam=$prefix.cleaned.bam
python $codedir/GetCBC_UMI_Ann.py $bamfile $outbam 

echo " "
echo Next Get SNP information!
samtools view  -F 2304 -h $outbam | grep -v _ANTISENSE | grep -v XF:Z:AMBIGUOUS | grep -v XF:Z:INTERGENIC | samtools view -Sb > ${prefix}.snps.bam

##Format for scan_snps
sed 's/|/\//g' $codedir/ref/SNPs.vcf | gzip > $prefix.snps.vcf.gz 

#/ahg/regevdata/users/libo/software/scan_snps/scan_snp $vcffile ${prefix}.snps.bam ${prefix}.snps
scan_snp $prefix.snps.vcf.gz ${prefix}.snps.bam ${prefix}.snps


echo " "
echo Get gene count matrix!
python $codedir/GetCounts.py ${prefix}.snps.bam ${prefix}.counts

echo Get Allele Count Matrix!
#export CLASSPATH=/ahg/regevdata/projects/612-eqtl/eQTL_testing/Code/Phasing:/ahg/regevdata/projects/612-eqtl/eQTL_testing/Code/Phasing/HTSJDK/htsjdk/build/libs/htsjdk-2.21.1-3-g3df5a35-SNAPSHOT.jar
export CLASSPATH=$codedir/Phasing:$codedir/Phasing/HTSJDK/htsjdk/build/libs/htsjdk-2.21.1-3-g3df5a35-SNAPSHOT.jar
java GetPhasing/CountPhased ${prefix}.snps.bam $vcffile $ind $prefix.phasin.txt 

echo " "
echo "Collect QC!"
echo First for reads!
samtools view -F 2304 $prefix.GE.bam | grep 'XF:Z:' | awk -F 'XF:Z:' '{print $2}' | awk '{print $1}' | sort | uniq -c > ${prefix}.counts.QC.all.Reads.txt
samtools view $prefix.GE.bam | awk '{print $1}' | sort | uniq | wc -l  >> ${prefix}.counts.QC.all.Reads.txt 
samtools view -F 2304 $prefix.cleaned.bam | grep 'XF:Z:' | awk -F 'XF:Z:' '{print $2}' | awk '{print $1}' | sort | uniq -c > ${prefix}.counts.QC.clean.Reads.txt
samtools view $prefix.cleaned.bam | awk '{print $1}' | sort | uniq | wc -l  >> ${prefix}.counts.QC.clean.Reads.txt
samtools view -f 4 $prefix.GE.bam | wc -l >> ${prefix}.counts.QC.all.Reads.txt
samtools view $prefix.snps.bam | wc -l > ${prefix}.useful.txt

echo Gene Expression and SNP info
Rscript $codedir/GetQC.R $prefix
Rscript $codedir/QC.Phasing.R $prefix
#Rscript $codedir/QC.Phasing.R $prefix



