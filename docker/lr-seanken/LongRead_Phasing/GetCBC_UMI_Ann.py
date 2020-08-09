import pysam
from gzip import GzipFile
import sys

##Reverse complements a string
def revComp(bases):
	bases=[MapBases[b] for b in list(bases)]
	bases.reverse()
	bases=''.join(bases)
	return(bases)


##
##Extracts UMI and CBC information from bam file (bamfile), outputs it is in outfile (outfile), correcting the CBC using the list of cbc in barcodes
##
def GetCBC(bamfile,outfile,isPacBio=False,toCorrectCBC=True):
	poss_chr=["chr"+str(i) for i in range(1,23)]
	poss_chr.extend(["chrX","chrY","chrM"])
	samfile = pysam.AlignmentFile(bamfile, "rb")
	outbam=pysam.AlignmentFile(outfile,"wb",template=samfile)

	print("Make dictionary for fast search!")
	#bars=FastBarcodeSearcher(barcodes) ##loads the barcodes into a datastructure that is fast to search within a certain edit distance

	i=0

	dists=[]
	noseq=0
	short=0
	nomatch=0
	nopolyA=0
	longer=0
	match=0
	noninter=0
	unmap=0
	flip=0
	noflip=0
	anti=0
	antiFlip=0
	for seq in samfile:
		ifFlip="not"
		i=i+1
		#seq_new=seq
		if i%100000 == 0:
			print(i);	
			print(" ")
		bases=seq.seq

		if seq.tid < 0:
			unmap=unmap+1
			continue;
		
		if not seq.has_tag("ZU"):
			continue;
		if not seq.has_tag("CB"):
			continue;

		umi=seq.get_tag("ZU")
		
		
		if not seq.is_reverse and umi in bases:
			ifFlip="flipped"

		if seq.is_reverse and not umi in bases:
			ifFlip="flipped"

		seq.set_tag("UB",umi,"Z")
		seq.set_tag("FL",ifFlip,"Z")
			
		##Update to include antisense information in annotation
		lab=seq.get_tag("XF")
		if lab!="INTERGENIC":
			noninter=noninter+1
			strand_gene=list(set(seq.get_tag("gs").split(","))) ##strand genes
			genes=list(set(seq.get_tag("gn").split(","))) ##strand genes
			

			#strand=seq.get_tag("ts") ##strand read
			strand="+"
			if seq.is_reverse:
				strand="-";
			revStrand={"+":"-","-":"+"}
			if ifFlip=="flipped":
				strand=revStrand[strand]
			if strand!=strand_gene[0]:
				lab=lab+"_ANTISENSE"
				anti=anti+1
				#FLIP if ifFlip=="flipped":
				#FLIP 	antiFlip=antiFlip+1
			if len(strand_gene)>1 or len(genes)>1:
				lab="AMBIGUOUS";
			seq.set_tag("XF",lab,"Z")
			seq.set_tag("SG",strand,"Z")
			seq.set_tag("GN",genes[0],"Z")

		##Fix CIGAR, replave X and = with M
		cigar=seq.cigarstring
		fixed_cigar=cigar.replace("X","M").replace("=","M")
		seq.cigarstring=fixed_cigar
		
		outbam.write(seq)	
	outbam.close()
	samfile.close()



if __name__=="__main__":
	args=sys.argv
	bamfile=args[1]
	outfile=args[2]
	isPacBio=False
	#isPacBio=args[4]
	print("Read in and fix barcodes!")
	

	print(args)

	print("Make it!")
	GetCBC(bamfile,outfile,isPacBio)
	print("Done!")
