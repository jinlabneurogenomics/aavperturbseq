library(Seurat)
library(scGBM)


Run_scGBM=function(seur,numDims=15,nCells=NULL,downstream=T,assay="RNA")
{
    dat=seur[[assay]]@counts[seur[[assay]]@var.features,]
    #dat=seur@assays$RNA@counts[seur@assays$RNA@var.features,]
    dat=as.matrix(dat)
    out.proj <- gbm.sc(dat,M=numDims,subset=nCells) 
    seur[["gbm"]] <- CreateDimReducObject(embeddings=out.proj$V,key="GBM_")
    rownames(seur@reductions$gbm@cell.embeddings)=names(seur@active.ident)
    if(downstream)
    {
        seur=RunUMAP(seur,dims=1:numDims,reduction="gbm")
        seur=FindNeighbors(seur,dims=1:numDims,reduction="gbm");seur=FindClusters(seur)
    }
    return(seur)
}
