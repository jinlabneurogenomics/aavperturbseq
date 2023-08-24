source("Helpers/RemoveDoublet.R")  
source("Helpers/Run_scGBM.R")  
source("Helpers/load_Seurat.R")
source("Helpers/Azimuth.R")



##fils is array consisting of filtered count matrices from CelLRanger, with names Ch1 to Ch5
##ref is a reference for Azimuth
##qcTable is table of QC results from the cell level QC code
##Not included in thsi command: Labelling clusters as cell types, additional filtering, adding in guide labels
MakeSeurat<-function(fils,ref,qcTable)
{
    print("Load CellRanger output")
    data=Read10X(fils)
    dat=data[[1]] ##RNA-seq counts
    mat=data[[2]] ##Guide RNA counts

    print("Process")
    seur=dir10X(dat=dat) ##From Helpers/load_Seurat.R, does basic Seurat object construction
    seur[["gRNA"]]=CreateAssayObject(counts = mat)
    print("Doublet ID")
    seur=DoubletScores(seur) ##Runs scdds to ID doublets

    print("Run scGBM and cluster/UMAP")
    seur=Run_scGBM(seur,numDims=20,nCells=30000) ##does dimentionality reduction with scGBM, followed by UMAP and clustering

    print("Run Azimuth")
    meta=RunAzimuth(seur,ref) ##uses Azimuth to ID clusters with help of Allen Brain dataset
    seur@meta.data=meta


    print("Add additional QC")
    for(i in colnames(qcTable)){seur@meta.data[i]=qcTable[names(seur@active.ident),i]}
    mito=grep("mt-",rownames(seur@assays$RNA@counts))
    seur@meta.data["Mito"]=100*colSums(seur@assays$RNA@counts[mito,])/colSums(seur@assays$RNA@counts)

    return(seur)

}