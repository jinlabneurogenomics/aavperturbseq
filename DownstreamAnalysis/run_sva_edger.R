# sva + edgeR script


RunSVAedgeR <- function(dat.bulk, meta.dat, form="~assignment+orig.ident", form2="~orig.ident", name='test', ref='GFP') {
require(edgeR)
require(stringr)
require(sva)
require(EnhancedVolcano)
# LRT models
mod = model.matrix(as.formula(form), data=meta.dat)  # full model
mod0 = model.matrix(as.formula(form2), data=meta.dat)  # null model
# SVA
#n.sv = num.sv(as.matrix(dat.bulk), mod)  # number of SVs
n.sv=1
modSv = mod
#if (n.sv > 0) {
#if (n.sv > 1) {print("Only include one surrogate variable");n.sv=1}
tryCatch(
	{
	svobj = sva(dat.bulk, mod, mod0, n.sv=n.sv)
	svs = svobj$sv
	colnames(svs) = paste0("SV", 1:ncol(svs))
	modSv = cbind(mod, svs)  # add SVs to model design
	},
	error=function(cond) {
		message("sva failed revert to using the provided form")
	}
)
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


