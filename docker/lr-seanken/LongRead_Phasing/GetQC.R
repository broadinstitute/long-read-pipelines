library(dplyr)
library(tidyr)
library(Matrix)
library(openxlsx)

args = commandArgs(trailingOnly=TRUE)

prefix=args[1]

print("Read in data about mapping rates, etc!")

QC_all=read.table(paste(prefix,".counts.QC.all.Reads.txt",sep=""),fill=T,stringsAsFactors=F)
QC_all[5,2]="All Reads"
QC_all[6,2]="Unmapped Reads"
QC_all[c(1,3,4),2]=sub("$"," (Sense or Antisense)",QC_all[c(1,3,4),2])
colnames(QC_all)=c("Counts_Unique_Reads","Type")
QC_all["Percent"]=QC_all[,1]/QC_all[5,1]
QC_all=QC_all[,c(2,1,3)]

QC_CBC=read.table(paste(prefix,".counts.QC.clean.Reads.txt",sep=""),fill=T,stringsAsFactors=F)

num=dim(QC_CBC)[1]

QC_CBC[num,2]="All Reads with CBC"

colnames(QC_CBC)=c("Counts_Unique_Reads","Type")
QC_CBC["Percent"]=QC_CBC[,1]/QC_CBC[num,1]
QC_CBC=QC_CBC[,c(2,1,3)]



print("Get Expression information!")
dat=readMM(paste(prefix,".counts.counts.mtx",sep=""))

numGenes=dim(dat)[1]
numCells=dim(dat)[2]

res<-c(numCells,numGenes)
nam<-c("Number Cells","Number Genes")


nUMI=colSums(dat)
medUMI=median(nUMI)
meanUMI=mean(nUMI)
sdUMI=sd(nUMI)
numHigh=sum(nUMI>100)

res<-c(res,medUMI,meanUMI,sdUMI,numHigh)
nam<-c(nam,"median nUMI","mean nUMI","SD nUMI","Number cells with UMI>100")


nGene=colSums(dat>0)
medGene=median(nGene)
meanGene=mean(nGene)
sdGene=sd(nGene)
numHigh=sum(nGene>50)



res<-c(res,medGene,meanGene,sdGene,numHigh)
nam<-c(nam,"median nGene","mean nGene","SD nGene","Number cells with Gene>50")


Expres=data.frame(nam)
colnames(Expres)="QC"
Expres["Quant"]=res



print("Get SNP info!")


lst=list()
lst[["QC All Reads"]]=QC_all
lst[["QC CBC Reads"]]=QC_CBC
lst[["QC Expression"]]=Expres

write.xlsx(lst,paste(prefix,".QC.combined.xlsx",sep=""))
