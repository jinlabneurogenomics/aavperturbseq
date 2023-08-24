####
###This code starts with a matrix with expression info and creates a basic Seurat object with it
####
library(Seurat)
library(stringr)
library(Matrix)

#loads all 10X lanes from a given directort
dir10X<-function(seur=NULL,dir="",dat=NULL,lst=c(),makeSeurat=T,minGenes=500,regress=c("nCount_RNA"),species="human",num=c())
{
if(is.null(seur))
{
if(length(lst)==0 & is.null(dat))
{
print(paste("ls ",dir,"/*/*/outs/filt*/* | grep : | sed 's/://g'",sep=""))
lst=system(paste("ls ",dir,"/*/*/outs/filt*/* | grep : | sed 's/://g'",sep=""),intern=T)

}
if(length(num)>0)
{
lst=lst[num]
}
print(lst)
if(is.null(dat))
{
print("Read in!")
dat=Read10X(lst)
print("Fix colnames!")
cols=colnames(dat)
cols_new=c()
for(col in cols)
{
start=str_sub(col,1,1)
cur=col
if(start %in% c("A","T","G","C")){cur=paste("1_",cur,sep="")}
cols_new<-c(cols_new,cur)
}

colnames(dat)=cols_new
#print("Return!")
}
print(paste("Dims: ",toString(dim(dat))))
if(!makeSeurat){return(dat)}

print("Make object!")
seur<-CreateSeuratObject(dat,"Seurat",min.features=minGenes)#,normalization.method="LogNormalize",scale.factor=1000000)
}

seur<-NormalizeData(seur,normalization.method="LogNormalize",scale.factor=10000)

print("Get variable genes!")
seur<-FindVariableFeatures(seur)


print("Regress out!")
if(length(regress)>0)
{
seur<-ScaleData(seur,features=seur@assays$RNA@var.features,vars.to.regress=regress)

}
else{
seur<-ScaleData(seur,features=seur@assays$RNA@var.features)
}
print("Run PCA!")
seur<-RunPCA(seur,npcs=60)


if(species=="mouse")
{
print("Get ribosomal and mito")
bot=Matrix::colSums(seur@assays$RNA@counts)
mito=grep("^mt-",rownames(seur@assays$RNA@counts))
if(length(mito)>0)
{
top=Matrix::colSums(seur@assays$RNA@counts[mito,])
mn=top/bot
seur@meta.data["mito"]=mn[names(seur@active.ident)]
}

mito=grep("^Rp[s,l]",rownames(seur@assays$RNA@counts))
if(length(mito)>0)
{
top=Matrix::colSums(seur@assays$RNA@counts[mito,])
mn=top/bot
seur@meta.data["ribo"]=mn[names(seur@active.ident)]
}


}

print("load!")


print("Return!")

return(seur)

}



