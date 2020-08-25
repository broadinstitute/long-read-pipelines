#!/bin/bash


#$ -cwd
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=90g
#$ -e ErrFiles/err.ann.err
#$ -o ErrFiles/log.ann.out
#$ -l os=RedHat7
#$ -l h_rt=40:00:00


source /broad/software/scripts/useuse


##Load the needed packages/ dotkits
use UGER
use .htslib-1.8
use GCC-5.2 
use .samtools-1.8
use BEDTools
use UGER
use .java-jdk-1.8.0_181-x86-64
use R-3.5
use .python-3.6.0 


#Point to R libraries
export R_LIBS_USER=/psych/genetics_data/ssimmons/R/x86_64-pc-linux-gnu-library/3.5_v2

##Add drop-seq code and scan_snps to PATH
export PATH=$PATH:Dropseq/Drop-seq_tools-2.3.0:/ahg/regevdata/users/libo/software/scan_snps

#Take in inputs
bamfile=$1
vcffile=$2
ind=$3
prefix=$4


echo $bamfile
echo $vcffile
echo $prefix

##sets the current directory as the code directory
codedir=$PWD



echo First get gene information, ignoring strand
outbam=$prefix.GE.bam

gtf=$codedir/ref/genes.chr.gtf


##Uses drop-seq code to tag reads with genes, ignoring strand info for now
TagReadWithGeneFunction I=$bamfile O=$outbam ANNOTATIONS_FILE=$gtf USE_STRAND_INFO=false


bamfile=$outbam
echo " "
echo Nice Make Nice bam, 10X like!
outbam=$prefix.cleaned.bam

##Cleans up bam for later steps
python $codedir/GetCBC_UMI_Ann.py $bamfile $outbam 

echo " "
echo Next Get SNP information!

##Gets reads that are uniquelly mapped, not antisense, not ambiguous, and not intergenic for SNP and phasing calling
samtools view  -F 2304 -h $outbam | grep -v _ANTISENSE | grep -v XF:Z:AMBIGUOUS | grep -v XF:Z:INTERGENIC | samtools view -Sb > ${prefix}.snps.bam

##Format for scan_snps
sed 's/|/\//g' $codedir/ref/SNPs.vcf | gzip > $prefix.snps.vcf.gz 

#/ahg/regevdata/users/libo/software/scan_snps/scan_snp $vcffile ${prefix}.snps.bam ${prefix}.snps
scan_snp $prefix.snps.vcf.gz ${prefix}.snps.bam ${prefix}.snps


echo " "
echo Get gene count matrix!
##Creates a gene by cell UMI count matrix
python $codedir/GetCounts.py ${prefix}.snps.bam ${prefix}.counts

echo Get Allele Count Matrix!
##Puts required code on class path for java
export CLASSPATH=$codedir/Phasing:$codedir/Phasing/HTSJDK/htsjdk/build/libs/htsjdk-2.21.1-3-g3df5a35-SNAPSHOT.jar

#Runs our code for getting phasing info for each gene in each cell
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



