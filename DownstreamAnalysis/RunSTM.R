library(stm)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(Seurat)
library(dplyr)
library(tidyr)
library(Matrix)

##
##Runs STM for a particular CellType
##
runSTM<-function(seur,celltypes,numberPCs=0,form=~sample,numTopics=10,genes=c(),getMeta=T)
{
print("Subset!")
if(length(celltypes)>0)
{
seur=subset(seur,cells=names(seur@active.ident)[seur@active.ident %in% celltypes])
}
print("Number Cells:")
print(length(seur@active.ident))


print("Remove common and rare genes, remove mito and ribo genes")
mn=rowMeans(seur@assays$RNA@counts>0)

dat=seur@assays$RNA@counts[mn>.05 & mn<.9,names(seur@active.ident)]


dat=dat[grep("^mt-",rownames(dat),invert=T),]
dat=dat[grep("^Rp[s,l]",rownames(dat),invert=T),]




if(length(genes)>0)
{
print("Use supplied gene list!")
genes=intersect(genes,rownames(dat))
dat=dat[genes,]
}



dat=t(dat)

print(dim(dat))


meta=seur@meta.data

print(form)

print("Run STM!")
print(dim(dat))
print(dim(meta))
res<-stm(dat,K=numTopics,prevalence=form,data=data.frame(meta),LDAbeta=T,interactions=F)
if(getMeta)
{
lst=list()
lst[["res"]]=res
lst[["meta"]]=meta
return(lst)
}
return(res)

}


##draws the Topic heatmap
drawRes<-function(res)
{
colnames(res)[2]="Log_Odds_Ratio"

mx=max(abs(res[,2]))+.01
p=ggplot(res,aes(x=factor(Topic),y=perturbation,fill=Log_Odds_Ratio))+geom_tile()+scale_fill_gradient2(low="blue",high="red",mid="white",lim=c(-mx,mx))+xlab("")+xlab("Topic")
return(p)
}





if(!interactive())
{
args=commandArgs(trailingOnly=TRUE)
seur_file=args[1]
outdir=args[2]
print("Read in seurat object!")
load(seur_file)
seur=excit
print("Clean up Perts")

print("Iterate through cell types with enough cells!")

seur@meta.data["nGene"]=scale(seur@meta.data[,"nGene"])


out<-runSTM(seur,c(),numberPCs=0,forForm=c("ko_vs_wt"),numTopics=10,genes=c(),minPert=10,getMeta=T)

save(out,file=paste(outdir,"/stm.Robj",sep=""))
print("next!")
}





GetTopGenes<-function(res,n)
{
lab=labelTopics(res,n=n)


frex=data.frame(lab$frex)
colnames(frex)=sub("^","Gene",as.character(1:dim(frex)[2]))
rownames(frex)=sub("^","Topic",as.character(1:dim(frex)[1]))


n=dim(frex)[2]

top=data.frame(lab$prob)
colnames(top)=sub("^","Gene",as.character(1:n))
rownames(top)=sub("^","Topic",as.character(1:dim(top)[1]))

lift=data.frame(lab$lift)
colnames(lift)=sub("^","Gene",as.character(1:n))
rownames(lift)=sub("^","Topic",as.character(1:dim(lift)[1]))



score=data.frame(lab$score)
colnames(score)=sub("^","Gene",as.character(1:n))
rownames(score)=sub("^","Topic",as.character(1:dim(score)[1]))

res=list()

res[["top"]]=top
res[["lift"]]=lift
res[["score"]]=score
res[["frex"]]=frex

for(i in names(res)){dat=res[[i]];res[[i]]["Topic"]=rownames(dat);n=dim(res[[i]])[2];res[[i]]=res[[i]][,c(n,1:(n-1))]}

return(res)


}


