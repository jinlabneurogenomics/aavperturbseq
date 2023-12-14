library(Seurat)
library(glmnet)
library(dplyr);library(tidyr)
library(ggplot2)

HiddenGLMNET<-function(seur,condition,sample="orig.ident",reduction="pca",numPC=50,usePCA=T,KOGene="Gpr34",type.measure=c("auc","deviance"))
{
    print("Prepare data")
    genes=c()
    if(usePCA)
    {
        X=seur[[reduction]]@cell.embeddings
        X=data.frame(X)[,1:numPC]
    }
    else{
        X=as.data.frame(t(seur@assays$RNA@scale.data))
        #X=X[,colnames(X)!=KOGene]
        genes=colnames(X)
        colnames(X)=sub("^","PC_",1:dim(X)[2])
    }
    X["sample"]=seur@meta.data[,sample]
    X["Cond"]=seur@meta.data[,condition]
    print(dim(X))

    #form=as.formula(paste("Cond~",paste(grep("PC",colnames(X),value=T),collapse="+"),sep=""))

    print("Prepare CV folds")
    #tab<-X %>% group_by(Cond,sample) %>% summarise(NumCells=length(Cond)) %>%  as.data.frame()
    conds=unique(X$Cond)
    if(length(conds)<3)
    {
    print("Has to have more than 2 conditions")
    return()
    }
    #tab["Num"]=-1
    #for(cond in conds)
    #{
    #    tab[tab$Cond==cond,"Num"]=1:sum(tab$Cond==cond)
    #}
   

    #minNum=min(sum(tab$Cond==conds[1]),sum(tab$Cond==conds[2]))
    #tab[tab$Num>minNum,"Num"]=sample(1:minNum,sum(tab$Num>minNum),replace=T)

    #print(tab)

    #rownames(tab)=tab[,"sample"]
    foldIDs=as.numeric(X[,"sample"])
    #weights=max(tab[as.character(X[,"sample"]),"NumCells"])*1/tab[as.character(X[,"sample"]),"NumCells"]
    print(table(foldIDs))
    conds=unique(X[,"Cond"])
    print("Fit cv")
    #out=cv.glmnet(as.matrix(X[,1:(dim(X)[2]-2)]),as.numeric(X[,"Cond"]==conds[1]),weights=weights,foldid=foldIDs,family = "binomial",type.measure=type.measure)
    out=cv.glmnet(as.matrix(X[,1:(dim(X)[2]-2)]),X[,"Cond"],foldid=foldIDs,family="multinomial")
    #out=cv.glmnet(as.matrix(X[,grep("PC",colnames(X))]),as.numeric(X[,"Cond"]==conds[1]),weights=weights,foldid=foldIDs,family = "binomial",type.measure=type.measure)
    #pred=predict(out,as.matrix(X[,grep("PC",colnames(X))]),s="lambda.min")
    pred=predict(out,as.matrix(X[,1:(dim(X)[2]-2)]),s="lambda.min")
    ret=list()
    ret[["model"]]=out
    rows=rownames(pred)
    pred=apply(pred,2,c)
    rownames(pred)=rows
    ret[["Prediction"]]=pred
    ret[["gene"]]=genes
    ret[["data"]]=as.matrix(X[,1:(dim(X)[2]-2)])
    return(ret)

}



plot_feature2 <- function(obj, features, title="", nc=0, cols=NULL, filename="test.pdf",
        size=3, dev='pdf', res=300, alpha=NA, limits=NULL, center_cols=F, ...)
{
geneAndMeta = c(rownames(obj), colnames(obj@meta.data))
if (!all(features %in% geneAndMeta)) {
print("Some genes not found!")
print(setdiff(features, geneAndMeta))
features = intersect(features, geneAndMeta)
}
print(paste0("alpha=", alpha))
if(is.null(cols)){cols <- viridis(50)}
p <- FeaturePlot(obj, features=features, combine=F, ...)
for (i in 1:length(p))
{
if (!is.null(limits)) {
p[[i]] <- p[[i]] + NoAxes() + scale_colour_gradient2(low="#004b88", mid="gray90", high="red3", limits=limits)
} else if (center_cols) {
p[[i]] <- p[[i]] + NoAxes() + scale_colour_gradient2(low="#004b88", mid="gray90", high="red3")
} else {
p[[i]] <- p[[i]] + NoAxes() + scale_colour_gradientn(colours=cols)
}
if (!is.na(alpha)) {
print(paste0("alpha=", alpha))
p[[i]]$layers[[1]]$aes_params$alpha = alpha
print(p[[i]]$layers[[1]]$aes_params$alpha)
}
}
if (nc < 1) {
nc=ceiling(sqrt(length(features)))
}
print(paste0("n_columns=", nc))
plot.width = nc*size
plot.height = ceiling(length(features)/nc)*size
print(paste0("plot_dimensions=", plot.width, " x ", plot.height))
if (dev=='tiff') {
tiff(filename, width=plot.width, height=plot.height, units='in', res=res)
} else if (dev=='pdf') {
pdf(filename, width=plot.width, height=plot.height)
} else if (dev=='png') {
png(filename, width=plot.width, height=plot.height, units='in', res=res)
}
print(ggpubr::annotate_figure(
p = cowplot::plot_grid(plotlist=p, ncol=nc), top=ggpubr::text_grob(label=title, face='bold', size=20)))
dev.off()
}


