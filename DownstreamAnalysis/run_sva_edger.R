# sva + edgeR script


RunSVAedgeR <- function(seur, use.var.genes=F, min.reads=10, min.norm.exp=.1, name='test') {
require(Matrix.utils)
require(edgeR)
require(stringr)
require(sva)
require(EnhancedVolcano)
dat = seur@assays$RNA@data
genes = c()
if (use.var.genes) {
seur <- FindVariableFeatures(seur)
genes = VariableFeatures(seur)
} else {
counts = Matrix::rowSums(seur@assays$RNA@counts, na.rm = TRUE)
genes1=names(counts)[which(counts>min.reads)]
rmeans = Matrix::rowMeans(dat, na.rm = TRUE)
genes2=names(rmeans)[which(rmeans>min.norm.exp)]
genes = intersect(genes1, genes2)
}
# pseudobulk
dat.bulk = as.matrix(t(aggregate.Matrix(t(dat[genes,]), groupings=seur$sample, fun='sum')))
meta.dat = data.frame(orig.ident=substr(colnames(dat.bulk), 1, 3), 
		      assignment=gsub("Ch[0-9]_", "", colnames(dat.bulk)), 
		      sample=colnames(dat.bulk))
meta.dat$assignment = factor(meta.dat$assignment, levels=c('GFP', g))
# LRT models
mod = model.matrix(~assignment+orig.ident, data=meta.dat)  # full model
mod0 = model.matrix(~orig.ident, data=meta.dat)  # null model
# SVA
#n.sv = num.sv(as.matrix(dat.bulk), mod)  # number of SVs
n.sv=1
modSv = mod
#if (n.sv > 0) {
#if (n.sv > 1) {print("Only include one surrogate variable");n.sv=1}
svobj = sva(dat.bulk, mod, mod0, n.sv=n.sv)
svs = svobj$sv
colnames(svs) = paste0("SV", 1:ncol(svs))
modSv = cbind(mod, svs)  # add SVs to model design
#}
# edgeR
dge = DGEList(counts=dat.bulk, samples=meta.dat$sample, group=meta.dat$assignment)
dge = estimateDisp(dge, modSv)
fit = glmFit(dge, modSv)
res = glmLRT(fit, coef=2)  # coef=2
print(summary(decideTests(res)))
pval = res$table$PValue
# p-value histogram
pdf(paste0("hist.", name, '.pdf'))
print(hist(pval, main=name, ylab='nGene'))
dev.off()
# adjust p-value
padj = p.adjust(pval, method='BH')
res$table$padj = padj
# volcano plot
pdf(paste0("volcano.", name, ".pdf"))
print(EnhancedVolcano(res$table,
    lab = rownames(res$table),
    x = 'logFC',
    y = 'padj',
    pCutoff=0.05,
    FCcutoff=0.1,
    subtitle=name))  #bquote(italic(.(name)))
dev.off()
# save result table
write.table(res$table, paste0("edgeR_LRT_with_sva.", name, ".tsv"), sep='\t', quote=F, row.names=T, col.names=NA)
pval.genes = rownames(res$table)[which(pval<0.05)]
padj.genes = rownames(res$table)[which(padj<0.05)]
return(list(pval.genes, padj.genes))
}


