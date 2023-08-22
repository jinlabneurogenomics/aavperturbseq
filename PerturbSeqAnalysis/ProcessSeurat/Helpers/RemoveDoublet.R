library(Seurat)
library(class)
library(scater)
library(scds)

DoubletScores<-function(seur,group_by="orig.ident")
{
samps=unique(seur@meta.data[,group_by])
#seur@meta.data["scds_group"]=seur@meta.data[,group_by]
seur@meta.data["scds"]=-1
seur@meta.data["cxds"]=-1
seur@meta.data["bcds"]=-1

for(samp in samps)
{
print(samp)
seur1=subset(seur,cells=names(seur@active.ident)[seur@meta.data[,group_by]==samp])
sce=as.SingleCellExperiment(seur1)
logcounts(sce) = log1p(counts(sce))
sce = cxds(sce,retRes = TRUE)
sce = bcds(sce,retRes = TRUE,verb=F)
sce = cxds_bcds_hybrid(sce)
seur@meta.data[names(sce$hybrid_score),"scds"]=as.numeric(sce$hybrid_score)
seur@meta.data[names(sce$hybrid_score),"cxds"]=as.numeric(sce$cxds_score)
seur@meta.data[names(sce$hybrid_score),"bcds"]=as.numeric(sce$bcds_score)
}

return(seur)

}






DoubletScores_knn<-function(seur,group_by="orig.ident",k=50,npc=20)
{
samps=unique(seur@meta.data[,group_by])
#seur@meta.data["scds_group"]=seur@meta.data[,group_by]
seur@meta.data["knn"]=-1

for(samp in samps)
{
print(samp)
seur1=subset(seur,cells=names(seur@active.ident)[seur@meta.data[,group_by]==samp])
pc=seur1@reductions$pca@cell.embeddings[,1:npc]
seur1@meta.data["type"]=seur1@meta.data[,"demux_type"] 
seur1@meta.data[seur1@meta.data[,"demux_type"] =="unknown","type"]="doublet"
out=knn(pc,pc,seur1@meta.data[,"type"],prob=T,k=k)
val=attributes(out)$prob
val[out=="singlet"]=1-val[out=="singlet"]

names(val)=names(seur1@active.ident)
seur@meta.data[names(val),"knn"]=as.numeric(val)
}

return(seur)

}


CallDoublets<-function(seur,samp_loc="orig.ident")
{
samps=unique(seur@meta.data[,samp_loc])
seur@meta.data[,"num_score"]=-1
for(samp in samps)
{
print(samp)
meta=seur@meta.data
meta=meta[meta[,samp_loc]==samp,]
numDoublets=dim(meta)[1]*dim(meta)[1]/1000*.01
print(numDoublets)
meta=meta[order(meta$scrub,decreasing=T),]
meta["scrub_doublet"]=1
meta[1:numDoublets,"scrub_doublet"]=0

meta=meta[order(meta$scds,decreasing=T),]
meta["scds_doublet"]=1
meta[1:numDoublets,"scds_doublet"]=0

meta=meta[order(meta$knn,decreasing=T),]
meta["knn_doublet"]=1
meta[1:numDoublets,"knn_doublet"]=0

meta["num_score"]=meta[,"knn_doublet"]+meta[,"scds_doublet"]+meta[,"scrub_doublet"]
seur@meta.data[rownames(meta),"num_score"]=meta[,"num_score"]
}
seur@meta.data["Singlet"]=seur@meta.data[,"num_score"]>1
return(seur)
}


library(here)

getHere<-function(x)
{
print(getwd())
print(here())
}
