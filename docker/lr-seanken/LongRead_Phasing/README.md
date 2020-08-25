# LongRead_Phasing


This is the current pipeline for our downstream analysis of the annotated bam files (still very messy, etc).


Inparticular, assumes you are running the pipeline from the main directory. The script to run the entire pipeline is CombinedPipeline_annot.sh

# What the pipeline relies on

The pipeline assumes:
1) you have Samtools (used 1.8) in your path
2) you have python in your path
3) You have bedtools in your path
4) you have R (I used v3.5) installed
5) You have java (I used 1.8) installed
6) If you want to call SNPs, need scan_snps on your path (https://github.com/broadinstitute/demuxEMS/tree/scan_snps)

In python, assumes you have the following packages installed (note some of these are now not used, will clean up soon):
1) pysam
2) scipy
3) pandas
4) collections
5) Levenshtein
6) gzip

In R, it assume you have the following packages installed:
1) dplyr
2) tidyr
3) Matrix
4) openxlsx


The steps in the pipeline are as folows:

# Step 0: setup

To set up need to do a few things:
1) unzip ref/genes.chr.gtf.gz if not already done
2) If needed recompile java code, Phasing/GetPhasing/CountPhased.java (note need $codedir/Phasing/HTSJDK/htsjdk/build/libs/htsjdk-2.21.1-3-g3df5a35-SNAPSHOT.jar on the class path)

# Step 1: Read In Data

The pipeline takes:
1) An annotated bam file
2) A vcf file (uncompressed)
3) The individual in the vcf file corresponding to the bam file (0 indexed). In paricular, the vcf file contains the genotype info for numerous inidividuals (starting in column 9, I believe). So column 9 corresponds to individual 0, column 10 corresponds to individual 1, etc.
4) The output prefix. In particular, the pipeline produces lots of files, so the prefix is the common prefix used for all those files. Ussually I set a particular outdir, and set my prefix to $outdir/results



The command to run the pipeline (in bash script form) is

./CombinedPipeline.sh $bamfile $vcffile $ind $prefix

so, for example:

./CombinedPipeline.sh test.bam test.vcf 0 test/results

Note CombinedPipeline.sh is built to work with UGER

Note CombinedPipeline_annot.sh is meant to be used with UGER array jobs and a sample sheet.

# Step 2: Label with Gene Of Origin

The next step is to label each read with the gene of origin (or label as intergenic if overlap no genes). This is done using the program developed to do this in Drop-seq, TagReadWithGeneFunction. Note the labelling is not strand specific, will take care of that later. 

The output is a labelled bam file, $prefix.GE.bam.

# Step 3: Fix up for later steps

To use the tools that come later, the bam file has to be cleaned up. As of now, this is done with a simple python script GetCBC_UMI_Ann.py.


This script reads in the reads from the $prefix.GE.bam file one by one, and does the following:
1) Filters out unmapped reads
2) Filters out unnanotated reads
3) Renames the tag for UMI from ZU to UB
4) Determines if the CBC/UMI are at the beginning of the read or the end (this is important for determining if sense or antisense). Adds a flag FL containing this information.
5) Checks if antisense, and if so modifies the XF to indicate antisense (the XF flag, output by the drop-seq tool, tells oyou if a read is intronic, intergenic, etc)
6) Check if a read is ambiguous (overlaps multiple genes), and if so changes the XF flag to say ambiguous
7) Adds strand and gene information in the flags GN and SG
8) Finally, reformats the CIGAR string to be ok with scan_snps (change X and = to M)

Outputs the resulting bam file into $prefix.cleaned.bam


# Step 4: Gets SNP information (Optional)

The next step is to extract snp information. First, we create the bam used for this (${prefix}.snps.bam), which is formed by rempvong antisense reads, ambiguous reads, intergenic reads, and multimapped reads.

Then scan_snp is run to generate count matrices for SNPs (each has a prefix ${prefix}.counts) in MM format.

Note this step is optional: we use the phased UMI, not the SNPs themselves

# Step 5: Get Gene Counts

Uses the reads in ${prefix}.snps.bam to get counts of number of UMIs per cell mapping to each gene. Outputs a sparse matrix (prefix ${prefix}.counts)


# Step 6: Get Phased Gene Counts

Like the above, except counts the number of UMI mapping to the first and second allele for each gene in each cell. Saved in $prefix.phasin.txt.

# Step 7: Get QC

Creates excel files with some basic QC metrics

# Output files:


The output files of importance are:

$prefix.snps.bam: the alignments used for our downstream analysis

$prefix.phasin.txt: a dataframe with allelic info in it. Column 1 is CBC, column 2 is gene, column 3 is allele (1 or 2), and column 4 is the count (number UMI's from that cell/gene/allele)

$prefix.snps.matrix.alt.mtx: The matrix (in MM format) of UMI counts for each SNP with the number of UMi with the alternative allele. $prefix.snps.barcodes.tsv and $prefix.snps.snps.tsv contain the row and column info

$prefix.snps.matrix.ref.matrix: Same as above, except for the reference allele

$prefix.counts.counts.mtx: The matrix in MM format counting the number of UMI mapping to each gene in each cell, with $prefix.counts.genes.txt being the gene labels, $prefix.counts.cells.txt the cell labels.

$prefix.QC.Phasing.xlsx: QC information about the phasing. For each gene it counts the total number of UMI mapping to each allele, the number that are ambiguous, and the percent ambiguous.

$prefix.QC.combined.xlsx: QC information about the mapping (the first page includes all reads, the second page only annotated reads), counts the number of reads in each category. The 3rd page is QC about the expression of the cells (note: counts empty droplets as of now)


