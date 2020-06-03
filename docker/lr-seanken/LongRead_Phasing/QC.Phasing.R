library(dplyr);library(tidyr);library(Matrix);library(openxlsx)

prefix=commandArgs(trailingOnly=TRUE)[1]
count_allele=paste(prefix,".phasin.txt",sep="")
count_gene=paste(prefix,".counts.counts.mtx",sep="")
gene_file=paste(prefix,".counts.genes.txt",sep="")
savefile=paste(prefix,".QC.Phasing.xlsx",sep="")

print("Load data")
dat=read.table(count_allele)
colnames(dat)=c("CBC","Gene","Allele","Count")

print("Combine results")
tab<-dat %>% spread(Allele,Count,fill=0) %>% group_by(Gene) %>% summarise(All1=sum(All1),All2=sum(All2),Ambig=sum(Ambig),Tot=All1+All2) %>% as.data.frame()

print("Get Gene Expression Info")
dat=readMM(count_gene)
genes<-scan(gene_file,"")

tots=rowSums(dat)
names(tots)=genes
print(head(tab))
print(head(tots))
tab["UMI_Tot"]=tots[as.character(tab[,"Gene"])]

tab=tab[!is.na(tab[,"UMI_Tot"]),]

tab["Percent"]=-1

tab[tab[,"UMI_Tot"]>0,"Percent"]=tab[tab[,"UMI_Tot"]>0,"Tot"]/tab[tab[,"UMI_Tot"]>0,"UMI_Tot"]

tab2=data.frame(c("Number Allele 1","Number Allele 2","Number Phased UMI","Number Ambiguous"))

colnames(tab2)="Label"

tab2["Count"]=c(sum(tab[,"All1"]),sum(tab[,"All2"]),sum(tab[,"Tot"]),sum(tab[,"Ambig"]))

res=list()
res[["Global"]]=tab2
res[["Local"]]=tab

print("Save!")
write.xlsx(res,savefile)
