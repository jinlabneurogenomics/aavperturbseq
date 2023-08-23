# code from Sean Simmons and Kira Perzel Mendell


library(tidyr)
library(dplyr)
library(clusterProfiler)
library(fgsea)

##
##Performs GSEA with given genes in universe univ
##Requires:
##-array of scores, with genes as names or dataframe df, with one columns for genes (specified by genes), one columns for Scores (specificed by Score
##-
RunEnrichment<-function(scores=c(),df=NULL,genes="Genes",Score="ES",Terms=NULL,CleanTerms=NULL,maxSize=500)
{
if(is.null(Terms))
{
print("Load Mouse GO terms")

}


if(length(scores)<5)
{
print("Get Scores from data frame")
scores=as.numeric(df[,Score])
names(scores)=df[,genes]
print(head(scores))
}

print("Clean up go terms")
tab=Terms #For later

if(is.null(CleanTerms))
{
Terms2<-lapply(unique(Terms[,"pathway"]),function(x){return(intersect(names(scores),Terms[Terms[,"pathway"]==x,"Gene"]))})

names(Terms2)=unique(Terms[,"pathway"])

Terms=Terms2
}
else
{
Terms=lapply(CleanTerms,function(x){return(intersect(names(scores),x))})
}


print("Run Enrichment!")
out<-data.frame(fgsea(Terms,scores,minSize=15,maxSize=maxSize,nproc=1))
#out<-data.frame(fgsea(Terms,scores,minSize=15,maxSize=maxSize))
#out<-data.frame(fgsea(Terms,scores,nperm=10000,minSize=15,maxSize=500))
if(is.null(tab))
{
out=out[order(out[,"pval"]),]
return(out)
}
print("Clean up output!")
tab=tab[,intersect(c("pathway","description","ont"),colnames(tab))]
tab=tab[!duplicated(tab[,"pathway"]),]
out<-left_join(out,tab)
out=out[order(out[,"pval"]),]
return(out)

}


getKegg<-function(organism="mmu",orgDB=NULL)
{

if(is.null(orgDB))
{
orgDB="org.Mm.eg.db"

if(organism=="hsa")
{
orgDB="org.Hs.eg.db"
}

if(organism=="rno")
{
orgDB="org.Rn.eg.db"
}
}
print("Download Kegg")
tab<-download_KEGG(organism)

tab<-left_join(tab[[1]],tab[[2]],by="from")

print("Get Nice Names")
colnames(tab)=c("pathway","Gene_ENT","description")

tab=data.frame(tab)

res<-bitr(tab[,"Gene_ENT"],fromType="ENTREZID",toType="SYMBOL",OrgDb=orgDB)


colnames(res)=c("Gene_ENT","Gene")

res=data.frame(res)
tab=left_join(tab,res)

tab["ont"]="KEGG"

return(tab)

}


RunEnrichment_KEGG<-function(scores = c(), df = NULL, genes = "Genes", Score = "ES",organism="mmu",orgDB=NULL)
{
print("Download Kegg")
tab<-getKegg(organism,orgDB)
print("Enrich")

print(head(tab))
out<-RunEnrichment(scores=scores,df=df,genes=genes,Score=Score,Terms=tab)

return(out)
}


getGO<-function(organism,orgDB)
{
    if(is.null(orgDB))
    {
    library(org.Mm.eg.db)
    orgDB=org.Mm.eg.db

    if(organism=="hsa")
    {
    library(org.Hs.eg.db)
    orgDB=org.Hs.eg.db
    }

    if(organism=="rno")
    {
    library(org.Rn.eg.db)
    orgDB=org.Rn.eg.db
    }
    }

    goterms <- AnnotationDbi::Ontology(GO.db::GOTERM)
    go2gene=AnnotationDbi::mapIds(orgDB, keys=names(goterms), column="SYMBOL",keytype="GOALL", multiVals='list')
    go2gene=go2gene[-which(sapply(go2gene, function(x){is.na(x)[1]}))]
    go2gene=lapply(go2gene,function(x){unique(x)})
    return(go2gene)

}

RunEnrichment_GO_new<-function(scores = c(), df = NULL, genes = "Genes", Score = "ES",organism="mmu",orgDB=NULL, maxSize=500, tab=NULL)
{
if (is.null(tab)) {
print("Download GO")
tab<-getGO(organism,orgDB)
}

print("Enrich")
##lapplyprint(head(tab))
out<-RunEnrichment(scores=scores,df=df,genes=genes,Score=Score,CleanTerms=tab,maxSize=maxSize)
print("Get ONTs")
onts=go2ont(out[,"pathway"])
rownames(onts)=onts[,1]
out["ont"]=onts[out[,"pathway"],2]
print("Get terms")
terms=go2term(out[,"pathway"])
rownames(terms)=terms[,1]
out["desc"]=terms[out[,"pathway"],2]
return(out)
}



RunEnrichment_GO<-function(scores = c(), df = NULL, genes = "Genes", Score = "ES",organism="mmu")
{
tab=NULL
Terms=NULL
if(organism=="mmu")
{
print("Load table")
load("/stanley/levin_dr/ssimmons/SingleCell2/GOEnrich/GO.output.Robj")
print("Load Terms")
load("/stanley/levin_dr/ssimmons/SingleCell2/GOEnrich/GO.Clean.Robj")
}
if(is.null(tab)){return(NULL)}

out<-RunEnrichment(scores=scores,df=df,genes=genes,Score=Score,Terms=tab,CleanTerms=Terms)

return(out)

}


RunEnrichment_All<-function(scores = c(), df = NULL, genes = "Genes", Score = "ES",organism="mmu")
{
out1<-RunEnrichment_GO(scores = c(), df = NULL, genes = "Genes", Score = "ES",organism="mmu")
out2<-RunEnrichment_KEGG(scores = c(), df = NULL, genes = "Genes", Score = "ES",organism="mmu")

out=rbind(out1,out2)

return(out)

}



