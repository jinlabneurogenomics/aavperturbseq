library(Seurat)
library(glmnet)
library(dplyr);library(tidyr)

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
