# conda activate hdWGCNA

print("Set up!")
library(qs)
library(Seurat)
library(dplyr)
library(ggplot2)
source("~/kwanho/src/seurat_tools.R")


# load data
seur <- qread("seur.for.downstream.Jul28.qs")
Idents(seur) = 'CellType'
seur <- subset(seur, idents='Inhib_Meis2', invert=T)

########################################
# Proportional change
########################################
setwd("proportion_change")
print("Propeller!")
source("propeller.R")
library(glmnet)
library(stringr)
library(reshape2)

seur$pert = seur$assignment
seur$pert = gsub("_[123]", "", seur$pert)
seur$pert = factor(seur$pert, levels=c('GFP','NonTarget1','NonTarget2','SafeTarget','Foxg1','Nr2f1','Tbr1','Tcf4'))
levels(seur$pert) = c('Control','Control','Control','Control','Foxg1','Nr2f1','Tbr1','Tcf4')

guide_order = c('GFP','NonTarget1','NonTarget2','SafeTarget',
	'Foxg1_1','Foxg1_2','Foxg1_3',
	'Nr2f1_1','Nr2f1_2','Nr2f1_3',
	'Tbr1_1','Tbr1_2','Tbr1_3',
	'Tcf4_1','Tcf4_2','Tcf4_3')
celltype_order = c('CR',
	'Excit_Upper','Excit_L2 IT ENTl',
	'Excit_L5IT','Excit_L5 PT CTX','Excit_L5NP_CTX',
	'Excit_L6IT','Excit_L6CT_CTX','Excit_Car3','Excit_L6b CTX','Excit_L6b/CT ENT',
	'Inhib_Sst','Inhib_Lhx6+Sst-','Inhib_Id2')

seur$assignment = factor(seur$assignment, levels=guide_order)
seur$CellType = factor(seur$CellType, levels=celltype_order)
seur$sample = paste0(seur$orig.ident, '_', seur$assignment)
x = expand.grid(levels(seur$orig.ident), levels(seur$assignment))
samp.order = paste(x[,1], x[,2], sep='_')
seur$sample = factor(seur$sample, levels=samp.order)

# Propeller t-test with covariate
# remove cell types with <200 cells
Idents(seur) <- 'CellType'
rm.cts = names(table(seur$CellType))[which(table(seur$CellType)<200)]
seur2 = subset(seur, idents=rm.cts, invert=T)
seur2@meta.data = droplevels(seur2@meta.data)
Idents(seur2) <- 'assignment'
dir.create("DemuxEM")
setwd('DemuxEM')
ref = 'NonTarget2'  # one of control guides with the most guide RNA counts
# each guide vs ref
dir.create("1v1")
for (g in setdiff(guide_order, ref)) {
print(g)
sseur <- subset(seur2, idents=c(g, ref))
sseur@meta.data = droplevels(sseur@meta.data)
props <- getTransformedProps(sseur$CellType, sseur$sample)
batch = str_split(colnames(props$Proportions), '_', simplify=T)[,1]
treat = gsub("Ch[1-5]_", "", colnames(props$Proportions))
meta = data.frame(batch, treat)
meta$treat = relevel(as.factor(meta$treat), ref=ref)
meta$batch = as.factor(meta$batch)
res = RunPropeller(props, meta, form="~0+treat+batch", out.prefix=paste0("1v1/propeller_res_",g,"_v_", ref))
#res = RunPropeller(props, meta, form="~0+treat+batch+ngRNA", out.prefix=paste0("1v1/propeller_res_",g,"_v_", ref))
sig_ct = rownames(res)[which(res$FDR<0.05)]
df = melt(props$Proportions)
df$guide = gsub("Ch[1-5]_", "", df$sample)
tab=df %>% group_by(clusters, guide) %>% summarise(sd=sd(value), mean=mean(value), n=n()) %>% mutate( se=sd/sqrt(n)) %>% mutate( ic=se * qt((1-0.05)/2 + .5, n-1)) 
sig.index = which(levels(tab$clusters) %in% sig_ct)
pdf(paste0("1v1/bar_plot_with_error_bar_", g, "_v_", ref, ".pdf"), height=4, width=8)
print(MyBar(tab, sig.index) + ggtitle(g))
dev.off()
}

# collect res
fl = list.files('1v1', full.names=T)
fl = fl[grepl('tsv', fl)]
res.list = lapply(fl, read.table, sep='\t', header=T, row.names=1)
names(res.list) = gsub(paste0("propeller_res_|_v_",ref,".tsv"), "", basename(fl))
for (i in 1:length(res.list)) {
nam = names(res.list)[i]
tab = res.list[[i]]
tab = tab %>% tibble::rownames_to_column('CellType')
tab[['guide']] = nam
colnames(tab)[grep("\\.\\.\\.", colnames(tab))] = "limma_coef"
nam.col = setdiff(colnames(tab), c(paste0('PropMean.treat', ref), 'limma_coef','PropRatio','Tstatistic','P.Value','FDR','CellType','guide'))
idx.nam.col = which(colnames(tab) == nam.col)
colnames(tab)[idx.nam.col] = "PropMean.treatGuide"
col.order = c('CellType','guide',paste0('PropMean.treat', ref),'PropMean.treatGuide','PropRatio','limma_coef','Tstatistic','P.Value','FDR')
tab = tab[,col.order]
res.list[[i]] = tab
}
res.all = do.call(rbind, res.list)
rownames(res.all) = paste0(res.all$guide, '_', res.all$CellType)
res.all = res.all %>% arrange(FDR)
write.table(res.all, paste0("propeller_res_combine_1v1_to_", ref, ".tsv"), sep='\t', quote=F, col.names=NA, row.names=T)


# summary plot
library(ComplexHeatmap)
rdbu = c('#67001F', '#B2182B', '#D6604D', '#F4A582', '#FDDBC7', '#FFFFFF','#D1E5F0', '#92C5DE', '#4393C3', '#2166AC', '#053061')
cols = rev(colorRampPalette(rdbu)(100))
#res.all$guide = factor(res.all$guide, levels=rev(setdiff(guide_order, ref)))
res.all$guide = factor(res.all$guide, levels=setdiff(guide_order, ref))
res.all$CellType = factor(res.all$CellType, levels=celltype_order)
res.all$PropRatio[res.all$PropRatio==0] = 0.00001
res.all = res.all %>% mutate(log2Ratio = log(PropRatio, base=2))
res.all$log2Ratio[which(res.all$log2Ratio>2)] = 2
res.all$log2Ratio[which(res.all$log2Ratio<(-2))] = -2
#val = dcast(res.all, guide~CellType, value.var='limma_coef')
val = dcast(res.all, guide~CellType)
val = val %>% tibble::column_to_rownames("guide")
target = str_split(rownames(val),'_',simplify=T)[,1]
target[c(1,2,3)] = 'Control'
guideAnno = data.frame(Target=target)
guideAnno$Target = factor(guideAnno$Target, levels=c("Control","Foxg1","Nr2f1","Tbr1","Tcf4"))
guideAnnoCols = list(Target=brewer.pal(5, 'Set2'))
names(guideAnnoCols$Target) = levels(guideAnno$Target)
rowHA = HeatmapAnnotation(df=guideAnno, which='row', col=guideAnnoCols)
nc_tab = table(seur2$CellType, seur2$assignment)
nc_tab = nc_tab[colnames(val), rownames(val)]
seur_ref = subset(seur2, idents=ref)
ref_ncell = table(seur_ref$CellType, seur_ref$assignment)[colnames(val), ref]
nc_tab = cbind(ref_ncell, nc_tab)
colHA = columnAnnotation(nCells = anno_barplot(nc_tab, beside=T, which='column', attach=T, bar_width=.9, height=unit(1,'in'), 
			 gp=gpar(fill=c(ref_ncell='gray', rep(guideAnnoCols$Target,each=3)))))
tab_props = getTransformedProps(seur2$CellType, seur2$assignment)$Proportions
cell_props = cbind(ref=tab_props[,ref], tab_props[,setdiff(colnames(tab_props), ref)])
guide_colors = c(ref='gray', rep(guideAnnoCols$Target,each=3))
colHA2 = columnAnnotation(cell_props = anno_barplot(cell_props, beside=T, which='column', attach=T, bar_width=.9, height=unit(1,'in'), gp=gpar(fill=guide_colors)))
mat_sig = dcast(res.all, guide~CellType, value.var="FDR") %>% tibble::column_to_rownames("guide")
mat_sig = mat_sig<0.05
mat_sig[mat_sig==T] = '*'
mat_sig[mat_sig==F] = ''
p=Heatmap(val, cluster_rows=F, cluster_columns=F,
	name="Effect size",
        row_title="Guide", row_names_side='left',
	column_title="Cell type",
        left_annotation=rowHA,
	top_annotation=colHA2,
        show_column_names=T,
	column_title_side='bottom',
        row_split=guideAnno$Target, gap=unit(2,'mm'),
	cell_fun = function(j, i, x, y, w, h, col) {
		if (mat_sig[i,j] =='*') {
			grid.circle(x, y, r=min(unit.c(w, h)*0.5), gp=gpar(fill='black', col='white'))
		}
	})
pdf('heatmap_prop_change_summary_DemuxEM.pdf', height=10, width=12)
draw(p, padding=unit(c(1, 1, 2, 2), "cm"))
dev.off()




########################################
# DEG 
########################################
# SVA + edgeR
setwd("../DEG")
print("DEG!")
library(Matrix.utils)
source("run_sva_edger.R")

seur$sample = paste0(seur$orig.ident, '_', seur$assignment)
ref = 'NonTarget2'
min.reads=10
min.norm.exp=.1

guides = setdiff(levels(seur$assignment), ref)
cts = levels(seur$CellType)

pval.genes = c()
padj.genes = c()

Idents(seur) = 'assignment'
for (g in guides) {
print(g)
sseur <- subset(seur, idents=c(g, ref))
Idents(sseur) = 'CellType'
for (ct in cts) {
print(ct)
ssseur <- subset(sseur, idents=ct)
if (ncol(ssseur) < 50) {print('Not enough cells!');next}
dat = ssseur@assays$RNA@data
# filtering genes to test for DE
counts = Matrix::rowSums(ssseur@assays$RNA@counts, na.rm = TRUE)
genes1=names(counts)[which(counts>min.reads)]
rmeans = Matrix::rowMeans(dat, na.rm = TRUE)
genes2=names(rmeans)[which(rmeans>min.norm.exp)]
genes = intersect(genes1, genes2)
# pseudobulk
dat.bulk = as.matrix(t(aggregate.Matrix(t(dat[genes,]), groupings=ssseur$sample, fun='sum')))
# prep metadata
meta.dat = data.frame(orig.ident=substr(colnames(dat.bulk), 1, 3),
                      assignment=gsub("Ch[0-9]_", "", colnames(dat.bulk)),
                      sample=colnames(dat.bulk))
meta.dat$assignment = factor(meta.dat$assignment, levels=c(ref, g))
ssseur@meta.data = droplevels(ssseur@meta.data)
ngRNA = ssseur@meta.data[, c('nCount_Crispr', 'sample')] %>% group_by(sample) %>% summarise(mean=mean(nCount_Crispr)) %>% mutate(log_mean=log(mean, base=10)) %>% pull(log_mean)
gl = RunSVAedgeR(dat.bulk, meta.dat, form="~assignment+orig.ident+ngRNA", form2="~orig.ident+ngRNA", ref=ref, name=paste0(g, '.', gsub(" |/", "_", ct)))
pval.genes = unique(c(pval.genes, gl[[1]]))
padj.genes = unique(c(padj.genes, gl[[2]]))
}
}

saveRDS(pval.genes, "all_pval_sig_genes.rds")
saveRDS(padj.genes, "all_padj_sig_genes.rds")

# combined plot
library(EnhancedVolcano)
library(patchwork)
fs = list.files(pattern="edgeR")
deg.res = list()
for (f in fs) {
cur.g = str_split(f, '\\.', simplify=T)[,2]
cur.ct = str_split(f, '\\.', simplify=T)[,3]
cur.tab = read.table(f, sep='\t', header=T, row.names=1)
cur.tab$CellType = cur.ct
cur.tab$guide = cur.g
deg.res[[paste0(cur.g, '_', cur.ct)]] = cur.tab
}
deg.all = do.call(rbind, deg.res)

plist = list()
for (g in setdiff(levels(seur$assignment), ref)) {
for (ct in levels(seur$CellType)[gsub(' |/', '_', levels(seur$CellType)) %in% names(table(deg.all$CellType))]) {
name = paste0(g, '_', ct)
cur.deg = deg.all %>% filter(CellType==ct, guide==g)
if (nrow(cur.deg)==0) {plist[[name]] = NULL; next}
p=EnhancedVolcano(cur.deg,
    lab = rownames(cur.deg),
    x = 'logFC',
    y = 'padj',
    pCutoff=0.05,
    FCcutoff=0.1,
    subtitle=name)  #bquote(italic(.(name)))
plist[[name]] = p
}}



# GSEA (since not many genes come up as significant, decided to run GSEA instead of GO)
source("~/kwanho/src/Enrich_FGSEA_new.R")
library(org.Mm.eg.db)

fs = list.files(pattern=".tsv")
orgdb = getGO('mmu',org.Mm.eg.db)
for (ct in cts) {
#for (ct in cts[9:length(cts)]) {
print(ct)
sfs = grep(gsub(" |/", "_", ct), fs, value=T)
comb.gsea.file = paste0("GSEA_res_GO_combined_", gsub(" |/", "_", ct), ".rds")
if (!file.exists(comb.gsea.file)) {
gsea.list = list()
if (length(sfs)<1) next
for (f in sfs) {
	guide_assign = str_split(f, "\\.", simplify=T)[,2]
	outnam = paste0("GSEA_res_GO.", guide_assign, ".", gsub(" |/", "_", ct), ".rds")
	if (!file.exists(outnam)) {
		df = read.table(f, sep='\t', header=T, row.names=1) %>% tibble::rownames_to_column("gene")
		res = RunEnrichment_GO_new(df=df, genes="gene", Score="logFC", organism="mmu", tab=orgdb)
		res$guide = guide_assign
		saveRDS(res, outnam)
	} else {
		res <- readRDS(outnam)
	}
	gsea.list[[guide_assign]] = res
}
tab = do.call(rbind, gsea.list)
saveRDS(tab, paste0("GSEA_res_GO_combined_", gsub(" |/", "_", ct), ".rds"))
} else {
tab = readRDS(comb.gsea.file)
}
desc.to.plot = tab %>% filter(padj < 0.05) %>% group_by(guide) %>% top_n(n=50, wt=abs(NES)) %>% pull(desc) %>% unique()
tab.plot = tab %>% filter(desc %in% desc.to.plot)
tab.plot$guide = factor(tab.plot$guide, levels=guides)
if (nrow(tab.plot)>0) {
pdf(paste0("dotplot_GSEA_res_GO_", gsub(" |/", "_", ct), ".pdf"), height=length(table(tab.plot$desc)), width=9)
print(ggplot(tab.plot, aes(x=guide, y=desc, color=NES, size=-log10(padj))) +
        geom_point() +
        geom_point(data=tab.plot[tab.plot$padj<0.05,], pch=21, fill=NA, color='black', stroke=2) +
        scale_colour_gradient2(low='blue2', mid='white', high='red2') +
        theme_classic() +
	ggtitle(ct) +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              axis.text=element_text(size=12),
	      plot.title=element_text(size=20)))
dev.off()
}
}


