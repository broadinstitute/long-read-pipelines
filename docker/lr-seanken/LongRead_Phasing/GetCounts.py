from scipy.io import mmwrite
from scipy.sparse import csr_matrix
from collections import Counter
import sys
import pysam
import pandas
import os

def getCounts(bamfile,outprefix):
	samfile = pysam.AlignmentFile(bamfile, "rb")
	i=0
	j=0
	dct={}
	for seq in samfile:
		i=i+1;
		if i%100000 == 0:
			print(i);
		location=seq.get_tag("XF")
		if location=="INTERGENIC":
			continue;
		j=j+1
		cbc=seq.get_tag("CB")
		umi=seq.get_tag("UB")
		gene=seq.get_tag("gn")
		if len(list(set(gene.split(","))))>1:
			gene="Ambiguous"
		gene=list(set(gene.split(",")))[0]
	
		if location in ["AMBIGUOUS","ANTISENSE"]:
			gene="Ambiguous"
		if "ANTISENSE" in location:
			gene="Ambiguous"
		comb=cbc+"_"+umi+"_"+gene

		#res=dct.get(comb,"")
			
		#if len(res)<2:
		dct[comb]=1
		#elif res!=gene:
		#	dct[comb]="Ambiguous"

	samfile.close()
	print("Have it all loaded!")
	
	print("Get counts!")

	counts=[[comb.split("_")[0],comb.split("_")[2]] for comb in dct.keys()]

	dat=pandas.DataFrame(counts,columns=["CBC","Gene"])

	print(dat.shape)
	print(dat.head())

	dat["Count"]=1

	dat=dat.groupby(["CBC","Gene"]).sum().reset_index()
	print(dat.head())

	genes=list(set(dat["Gene"]))
	cells=list(set(dat["CBC"]))
	
	numGenes=len(genes);
	numCells=len(cells)
	
	GeneToNum={}
	CellToNum={}
	for i in range(0,len(genes)):
		GeneToNum[genes[i]]=(i)

	dat["IndexGene"]=[GeneToNum[t] for t in dat["Gene"]]

	for i in range(0,len(cells)):
		CellToNum[cells[i]]=(i)



	dat["IndexCell"]=[CellToNum[t] for t in dat["CBC"]]

	#dat.to_csv(outprefix+"temp.txt",sep="\t")
	dat=dat[["IndexGene","IndexCell","Count"]]

	dat=csr_matrix(([num for num in dat["Count"]],([gene for gene in dat["IndexGene"]],[cell for cell in dat["IndexCell"]])))

	#dat=dat.pivot(index='Gene',columns='CBC',values='Count').fillna(0)

	#genes=list(dat.index)
	#cells=list(dat.columns)
	#dat=csr_matrix(dat.values)
	
	print("Save!")
	mmwrite(outprefix+".counts.mtx",dat)	
	with open(outprefix+'.genes.txt', 'w') as f:
		for g in genes:
			f.write("%s\n" % g)

	with open(outprefix+'.cells.txt', 'w') as f:
		for g in cells:
			f.write("%s\n" % g)
		
	
	print("Done!")

if __name__=="__main__":
	args=sys.argv
	bamfile=args[1]
	outprefix=args[2]

	getCounts(bamfile,outprefix)
