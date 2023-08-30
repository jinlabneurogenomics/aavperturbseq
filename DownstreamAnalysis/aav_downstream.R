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
dir.create("min200_cov_gRNA_batch")
setwd('min200_cov_gRNA_batch')
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
#ngRNA = sseur@meta.data[, c('nCount_Crispr', 'sample')] %>% mutate(log_nGuides=log(nCount_Crispr, base=10)) %>% mutate(scale_nGuides=log_nGuides-mean(log_nGuides)) %>% group_by(sample) %>% summarise(mean=mean(scale_nGuides)) %>% pull(mean)
ngRNA = sseur@meta.data[, c('nCount_Crispr', 'sample')] %>% mutate(log_nGuides=log(nCount_Crispr, base=10)) %>% group_by(sample) %>% summarise(mean=mean(log_nGuides)) %>% pull(mean)
meta = data.frame(batch, treat, ngRNA)
meta$treat = relevel(as.factor(meta$treat), ref=ref)
meta$batch = as.factor(meta$batch)
#res = RunPropeller(props, meta, form="~0+treat+batch", out.prefix=paste0("1v1/propeller_res_",g,"_v_", ref))
res = RunPropeller(props, meta, form="~0+treat+batch+ngRNA", out.prefix=paste0("1v1/propeller_res_",g,"_v_", ref))
sig_ct = rownames(res)[which(res$adj.P.Val<0.05)]
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


# Plot
#PropsPlot(seur, "CellType", "assignment", name.prefix='plot_cell_type_props_per_guide')
#PropsPlot(seur, "CellType", "pert", name.prefix='plot_cell_type_props_per_gene')

# summary plot
library(ComplexHeatmap)
rdbu = c('#67001F', '#B2182B', '#D6604D', '#F4A582', '#FDDBC7', '#FFFFFF','#D1E5F0', '#92C5DE', '#4393C3', '#2166AC', '#053061')
cols = rev(colorRampPalette(rdbu)(100))
#res.all$guide = factor(res.all$guide, levels=rev(setdiff(guide_order, ref)))
res.all$guide = factor(res.all$guide, levels=setdiff(guide_order, ref))
res.all$CellType = factor(res.all$CellType, levels=celltype_order)
#res.all$PropRatio[res.all$PropRatio==0] = 0.00001
#res.all = res.all %>% mutate(log2Ratio = log(PropRatio, base=2))
#res.all$log2Ratio[which(res.all$log2Ratio>2)] = 2
#res.all$log2Ratio[which(res.all$log2Ratio<(-2))] = -2
val = dcast(res.all, guide~CellType, value.var='limma_coef')
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
ngRNA = seur2@meta.data[, c('nCount_Crispr', 'assignment', 'CellType')] %>% mutate(log_nGuides=log(nCount_Crispr, base=10)) %>% group_by(assignment,CellType) %>% summarise(mean=mean(log_nGuides))
guide_counts = dcast(ngRNA, CellType~assignment) %>% tibble::column_to_rownames("CellType")
guide_counts = cbind(ref=guide_counts[,ref],guide_counts[,setdiff(colnames(guide_counts), ref)])
guide_colors = c(ref='gray', rep(guideAnnoCols$Target,each=3))
colHA2 = columnAnnotation(cell_props = anno_barplot(cell_props, beside=T, which='column', attach=T, bar_width=.9, height=unit(1,'in'), gp=gpar(fill=guide_colors)),
			  guide_mean_log_counts = anno_barplot(guide_counts, beside=T, which='column', attach=T, bar_width=.9, height=unit(1,'in'),
                          gp=gpar(fill=guide_colors)))
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
		#grid.text(mat_sig[i,j], x, y)
		if (mat_sig[i,j] =='*') {
			grid.circle(x, y, r=min(unit.c(w, h)*0.5), gp=gpar(fill='black', col='white'))
			#grid.circle(x, y, r=min(unit.c(w, h)*0.5), gp=gpar(fill='black'))
		}
	})
#pdf('heatmap_prop_change_summary_with_nCells.pdf', height=10, width=10)
pdf('heatmap_prop_change_summary_with_props_limma_coef.pdf', height=10, width=12)
draw(p, padding=unit(c(1, 1, 2, 3), "cm"))
dev.off()

# Split up by lineage (Exn v IN) and do prop change testing


########################################
# Most affected cell populations
########################################
setwd("../../affected_cells")
print("HiDDEN!")
library(RColorBrewer)
source("mod.hidden.mult.R")

# UMAP cell types
pdf("umap_cell_type.pdf", height=9, width=9)
MyDimPlot(seur, alpha=1, group.by='CellType', pt.size=1, label=T, repel=T, label.size=5, label.box=T) + NoLegend()
dev.off()

# UMAP pert
pdf("umap_pert.pdf", height=10, width=9)
MyDimPlot(seur, alpha=1, group.by='pert', pt.size=1, legend.text.size=16)
dev.off()

# UMAP guide
pdf("umap_guide.pdf", height=10, width=9)
MyDimPlot(seur, alpha=1, group.by='assignment', pt.size=1, legend.text.size=16)
dev.off()

# Hidden method (multinomial version)
# guide level
ret = HiddenGLMNET(seur, 'assignment')
saveRDS(ret, "hidden.result.guide.rds")

pred = ret$Prediction
colnames(pred) = paste0("hidden.pred.", colnames(pred))
seur <- AddMetaData(seur, as.data.frame(pred))

lim = max(abs(max(pred)), abs(min(pred))) * c(-1,1)
feas=c(paste0('hidden.pred.',levels(seur$assignment))[5:16], paste0('hidden.pred.',levels(seur$assignment))[1:4])

plot_feature2(seur, feas,
        title="Most affected cells per condition",
        nc=3,
        limits=lim,
        size=4,
        dev='pdf',
        alpha=0.7,
        filename="hidden.pred.guide.color_v2.pdf")

plot_feature2(seur, feas,
        title="Most affected cells per condition",
        nc=3,
        #cols=colorRampPalette(colors=c("#004b88", "gray90", "red3"))(50),
        center_cols=T,
        size=4,
        dev='pdf',
        alpha=0.7,
        filename="hidden.pred.guide.color_v1.pdf")


# Correlation between each pair of HiDDEN scores
library(Hmisc)
library(corrplot)

res1 <- rcorr(pred, type='pearson')
res2 <- rcorr(pred, type='spearman')
#rdbu = c('#67001F', '#B2182B', '#D6604D', '#F4A582', '#FDDBC7', '#FFFFFF','#D1E5F0', '#92C5DE', '#4393C3', '#2166AC', '#053061')
#cols=rev(colorRampPalette(rdbu)(100))

pdf("heatmap.corr.pearson.pdf")
corrplot(res1$r, type = "upper", order = "original", tl.col = "black", tl.srt = 45,
        title='Pearson correlation of HiDDEN scores', mar=c(0,0,1,0))
dev.off()

pdf("heatmap.corr.spearman.pdf")
corrplot(res2$r, type = "upper", order = "original", tl.col = "black", tl.srt = 45,
        title='Spearman correlation of HiDDEN scores', mar=c(0,0,1,0))
dev.off()



# Milo (NN method)



########################################
# DEG 
########################################
# SVA + edgeR
setwd("../DEG")
print("DEG!")
library(Matrix.utils)
source("run_sva_edger.R")

dir.create('res_082823')
setwd("res_082823")

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
#rownames(cur.deg) = str_split(rownames(cur.deg), '\\.', simplify=T)[,2]
p=EnhancedVolcano(cur.deg,
    lab = rownames(cur.deg),
    x = 'logFC',
    y = 'padj',
    pCutoff=0.05,
    FCcutoff=0.1,
    subtitle=name)  #bquote(italic(.(name)))
plist[[name]] = p
}}
#pdf(paste0("volcano_combined.", name, ".pdf"), height=30, width=30)
#wrap_plots(plist, ncol=nlevels(seur$assignment)-1, byrow=T)
#dev.off()



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


########################################
# Gene module
########################################
# stm
setwd('../../gene_module/stm')
library(clusterProfiler)
library(gtools)
source("RunSTM.R")

#form = ~assignment + pert  # + orig.ident
form = ~1
#dir.create('form_guide+target')
#setwd('form_guide+target')
dir.create('form_1')
setwd('form_1')

seur$CellType_Gen = factor(seur$CellType_Gen, levels=c('CR','Excitatory','Inhibitory'))
levels(seur$CellType_Gen)[1] = 'Excitatory'

for (ct in levels(seur$CellType_Gen)) {
k = ifelse(ct=='Excitatory', 20, 10)

ct_str = gsub(" |/", "_", ct)
if (!dir.exists(ct_str)) dir.create(ct_str)  # create & change directory
setwd(ct_str)

Idents(seur) <- 'CellType_Gen'  # runSTM assums it
out.nam = paste0("stm_res_", ct_str, ".rds")
if (!file.exists(out.nam)) {
res = runSTM(seur, ct, form=form, numTopics=k, getMeta=T)
saveRDS(res, out.nam)
} else {
res = readRDS(out.nam)
}
meta = res[[2]]
res = res[[1]]

sseur <- subset(seur, idents=ct)
sc = res$theta
rownames(sc) = colnames(sseur)
colnames(sc) = paste0('STM_', 1:k)
sseur <- AddMetaData(sseur, metadata=as.data.frame(sc))

# stats
library(lme4)
library(lmerTest)
dat = as.data.frame(sc)
dat$guide = sseur$assignment[rownames(dat)]
dat$channel = sseur$orig.ident[rownames(dat)]
dat$ngRNA = sseur$nCount_Crispr
# clean table for stats
tab = reshape(dat, varying=1:k, direction='long', sep='_')
colnames(tab)[which(colnames(tab)=='time')] = "module"
colnames(tab)[which(colnames(tab)=='STM')] = "value"
tab$module = factor(paste0("STM_", tab$module), levels=colnames(sc))
#tab = melt(dat, variable.name='module')
saveRDS(tab, paste0("table_stm_theta_", ct, ".rds"))
ref='NonTarget2'
#df.padj = data.frame(matrix(NA,nrow=nlevels(tab$guide)-1, ncol=k),row.names=setdiff(levels(tab$guide), 'GFP'))
#ncomp = nlevels(tab$guide)-1+nlevels(tab$channel)-1
#rownam = c(setdiff(levels(tab$guide), 'GFP'), setdiff(levels(tab$channel), 'Ch1'))
ncomp = nlevels(tab$guide)-1
rownam = setdiff(levels(tab$guide), ref)
df.padj = data.frame(matrix(NA,nrow=ncomp, ncol=k),row.names=rownam)
colnames(df.padj) = colnames(sc)
for (cur.mod in levels(tab$module)) {
	print(cur.mod)
	stab = tab %>% filter(module == cur.mod) %>% mutate(log_ngRNA=log(ngRNA,base=10), scale_ngRNA=log_ngRNA-mean(log_ngRNA))
	#fit = lm(value~guide+channel, data=stab)
	fit = lmer(value~guide+scale_ngRNA+(1|channel:guide), data=stab)
	fit.res = summary(fit)
	#cur.pval = fit.res$coefficients[grep("guide|channel", rownames(fit.res$coefficients)), 4]
	cur.pval = fit.res$coefficients[grep("guide", rownames(fit.res$coefficients)), 5]
	df.padj[,cur.mod] = p.adjust(cur.pval, method='BH')
}
write.table(df.padj, paste0("stats_lmer_padj_stm_", ct, ".tsv"), sep='\t', quote=F, row.names=T, col.names=NA)
#df.padj = read.table(paste0("stats_lmer_padj_stm_", ct, ".tsv"), sep='\t', header=T, row.names=1)

# dotplot stats stm
tab.plot = tab %>% group_by(guide, module) %>% summarise(n=n(),mean_theta=mean(value),std=sd(value))
plot.tab = NULL
for (stm.mod in levels(tab.plot$module)) {
	cur.tab = tab.plot %>% filter(module==stm.mod)
	gfp.idx = which(cur.tab$guide==ref)
	gfp.val = cur.tab$mean_theta[gfp.idx]
	cur.tab$effect_size = (cur.tab$mean_theta - gfp.val) / cur.tab$std
	if (is.null(plot.tab)) {
		plot.tab = cur.tab
	} else {
		plot.tab = rbind(plot.tab, cur.tab)
	}
}
tab.plot = plot.tab
tab.plot$padj = NA
for (i in 1:nrow(tab.plot)) {
	if (tab.plot$guide[i] == ref) {
		tab.plot$padj[i] = 1
		next
	}
	tab.plot$padj[i] = df.padj[as.character(tab.plot$guide[i]), as.character(tab.plot$module)[i]]
}
tab.plot$padj[tab.plot$padj<2e-16] = 2e-16  # cap min p-value
tab.plot$guide = factor(tab.plot$guide, levels=rev(levels(tab$guide)))
wid = ifelse(ct=='Excitatory', 7,5)
rdbu = c('#67001F', '#B2182B', '#D6604D', '#F4A582', '#FDDBC7', '#FFFFFF','#D1E5F0', '#92C5DE', '#4393C3', '#2166AC', '#053061')
cols = rev(colorRampPalette(rdbu)(100))
tab.plot = tab.plot %>% filter(guide != ref)
pdf(paste0("dot_plot_stm_module_FC_lmer_pvalues_", ct, ".pdf"), height=5, width=wid)
print(tab.plot %>%
	ggplot(aes(x=module, y=guide, color=effect_size, size=-log10(padj))) +
	geom_point() +
	#scale_color_gradientn(colors=cols) +
	scale_colour_gradient2(low="#004b88", mid="gray90", high="red3") +
	geom_hline(yintercept=seq(3,12,3)+0.5, linetype='dashed') +
	geom_point(data=tab.plot[tab.plot$padj<0.05,], pch=1, fill=NA, color='black', stroke=1) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text=element_text(size=12)))
dev.off()

# violin plot theta
MyViolinPlot(sseur, features=colnames(sc), group.by="assignment", nc=1, 
	filename=paste0("violin_stm_theta_", ct, ".pdf"), height.multiplier=1.25)

# feature plot theta
plot_feature2(sseur, features=colnames(sc), title='STM theta', nc=4, filename=paste0("feature_plot_theta_", ct, ".pdf"))

# plot topic props (The proportion of words in a document that belong to a topic; i.e. find most prevalent module in a cell type)
topicNames = paste0("Topic ", 1:k)
pdf(paste0("topic_props_", ct_str, ".pdf"))
par(bty="n",col="black",lwd=5)
plot.STM(res,type="summary", n=3, main=ct)
plot.STM(res,type="summary",topic.names=topicNames, custom.labels='', main=ct)
dev.off()

# cloud plot of words
pdf(paste0('cloud_', ct_str, '.pdf'), height=14, width=15)
par(mfrow=c(3, 4), cex.main=2.5)
for (i in 1:10) {
print(cloud(res, topic=i, scale=c(3.5, .15), colors=brewer.pal(8, "Dark2"), max.words=50))
print(title(main=paste0('Topic ', i)))
}
dev.off()

# topic correlation
cor.cutoff = 0.1
topic_corr = topicCorr(res, method='simple', cutoff=cor.cutoff)
# graph
pdf("graph_corr_topics.pdf")
plot(topic_corr, vertex.size=25, vertex.color='gray90', main=paste0('correlation cutoff: ', cor.cutoff))
dev.off()
# heatmap
pdf("heatmap_corr_topics.pdf")
corrplot(topic_corr$cor, type = "upper", tl.col = "black", tl.srt = 45,
        title='Correlation of STM topics', mar=c(0,0,1,0))
corrplot(topic_corr$poscor, type = "upper", tl.col = "black", tl.srt = 45,
        title=paste0('Correlation of STM topics (r>', cor.cutoff, ')'), mar=c(0,0,1,0))
dev.off()

# Plot GO terms for the topics of interest
n=100  # number of genes to include per topic
top.genes = GetTopGenes(res, n=n)
tg.frex = top.genes$frex %>% dplyr::select(-Topic)
mod = c(as.data.frame(t(tg.frex)))
saveRDS(mod, paste0("gene_list_top", n, "_", ct, ".rds"))
# select interesting topics
topic.idx = c(2,3,4,6,7,9,13,17,19)  # Excitatory
if (ct == "Inhibitory") topic.idx = c(2,5,6,7,9,10)
mod = mod[topic.idx]
# Seurat module score
sseur <- MyModuleScore(sseur, mod, filename=paste0("metadata_Seurat_module_score_of_STM_top", n, "_genes_", ct, ".rds"))
# violin plot module score
MyViolinPlot(sseur, features=names(mod), group.by="assignment", nc=1,
        filename=paste0("violin_Seurat_module_score_of_STM_top", n, "_genes_", ct, ".pdf"), height.multiplier=1.25, hline.intercept=0)
# GO terms associated with topics
mn = rowMeans(sseur@assays$RNA@data)
bg.genes = names(mn)[mn>0.1]  # background genes: remove very lowly expressed genes
tl = list()
for (i in 1:length(mod)) {
genes = mod[[i]]
mod.nam = names(mod)[i]
print(mod.nam)
go.file = paste0("GO_res_frex_", mod.nam, "_", ct, ".rds")
if (!file.exists(go.file)) {
enrich_go <- enrichGO(gene=genes, universe=bg.genes, OrgDb=org.Mm.eg.db, keyType="SYMBOL", ont="ALL", pvalueCutoff=1, readable=F)
saveRDS(enrich_go, go.file)
} else {
enrich_go = readRDS(go.file)
}
print(dim(enrich_go))
if (nrow(enrich_go)>0) {
hei = 4+ceiling(ifelse(nrow(enrich_go)>45, 45, nrow(enrich_go))/3)
print(hei)
pdf(paste0("dotplot_GO_frex_", mod.nam, "_", ct, ".pdf"), width=8, height=hei)
print(dotplot(enrich_go, showCategory=45)+theme(axis.text.y=element_text(size=16))+ggtitle(paste0(ct, " ", names(mod)[i])))
dev.off()
tab = enrich_go@result
tab$topic = mod.nam
#write.table(tab, paste0("table_GO_res_frex_topic",i,".tsv"), sep='\t', quote=F, row.names=F, col.names=T)
tl[[i]] = tab
}
}
tab = do.call(rbind, tl)
write.table(tab, paste0("GO_res_combined_", ct, ".tsv"), sep='\t', quote=F, row.names=F, col.names=T)
desc.to.plot = tab %>% filter(p.adjust < 0.05) %>% group_by(topic) %>% 
	top_n(n=10, wt=Count) %>% pull(Description) %>% unique()
#topic_ord = paste0('Topic', unique(unlist(apply(topic_corr$posadj, 1, function(x) which(x>0)))))
topic_ord = mixedsort(names(table(tab$topic)))
tab.plot = tab %>% filter(Description %in% desc.to.plot) %>% mutate(xlab=factor(topic, levels=topic_ord))
outname = paste0("dotplot_GO_comparison_", ct, ".pdf")
pdf(outname, height=4+length(table(tab.plot$Description))/3, width=13)
print(ggplot(tab.plot, aes(x=xlab, y=Description, color=Count, size=-log10(p.adjust))) +
        geom_point() +
        geom_point(data=tab.plot[tab.plot$p.adjust<0.05,], pch=21, fill=NA, color='black', stroke=2) +
        scale_colour_viridis(direction=-1) +
	scale_y_discrete(labels=function(x) str_wrap(x, width=80)) +
	scale_size_continuous(range=c(1, ceiling(max(-log10(tab.plot$p.adjust)))), breaks=c(2, 4, 6, 8), name="p.adjust") +
	xlab("Topic") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              axis.text=element_text(size=20)))
dev.off()

# Collapse GO terms
source("~/kwanho/src/go_reduce.R")
colnames(tab)[1] = 'go_type'
colnames(tab)[2] = 'go_id'
go.list = split(tab, tab$topic)[topic_ord]
reduc.list = list()
for (i in 1:length(go.list)) {
df = go.list[[i]]
score.arr = -log10(df$p.adjust)
names(score.arr) = df$go_id
df.reduc = go_reduce(df, orgdb = "org.Mm.eg.db", threshold = 0.7, scores=score.arr)  # higher threshold = less terms
idx=c(); for (j in 1:nrow(df)) {if(df.reduc$go_id[j]==df.reduc$parent_id[j]) idx=c(idx,j)}
n.terms = table(df.reduc$parent_id)
sdf = df.reduc[idx,]
for (gi in sdf$go_id) {
cur.row = which(sdf$go_id == gi)
sdf$Description[cur.row] = paste0(sdf$Description[cur.row], " (", n.terms[gi], ")")
}
library(multienrichjam)
plot.enr <- enrichDF2enrichResult(enrichDF = sdf, keyColname = "go_id", geneColname = "geneID", 
				  pvalueColname = "p.adjust", descriptionColname = "Description", pvalueCutoff=1)
hei = 4+ceiling(ifelse(nrow(plot.enr)>45, 45, nrow(plot.enr))/3)
pdf(paste0("dotplot_GO_frex_reduced_", names(go.list)[i], "_", ct, ".pdf"), width=10, height=hei)
print(dotplot(plot.enr, showCategory=45, label_format=60)+
	theme(axis.text.y=element_text(size=16))+
	ggtitle(paste0(ct, " ", names(go.list)[i])))
dev.off()
reduc.list[[names(go.list)[i]]] = sdf
}
tab = do.call(rbind, reduc.list)
write.table(tab, paste0("GO_res_reduced_", ct, ".tsv"), sep='\t', quote=F, row.names=F, col.names=T)
desc.to.plot = tab %>% filter(p.adjust < 0.05) %>% group_by(topic) %>%
        top_n(n=10, wt=Count) %>% pull(parent_term) %>% unique()
tab.plot = tab %>% filter(parent_term %in% desc.to.plot) %>% mutate(xlab=factor(topic, levels=topic_ord))
outname = paste0("dotplot_GO_reduced_", ct, ".pdf")
pdf(outname, height=4+length(table(tab.plot$parent_term))/3, width=13)
print(ggplot(tab.plot, aes(x=xlab, y=parent_term, color=Count, size=-log10(p.adjust))) +
        geom_point() +
        geom_point(data=tab.plot[tab.plot$p.adjust<0.05,], pch=21, fill=NA, color='black', stroke=2) +
        scale_colour_viridis(direction=-1) +
        scale_y_discrete(labels=function(x) str_wrap(x, width=80)) +
        scale_size_continuous(range=c(1, ceiling(max(-log10(tab.plot$p.adjust)))), breaks=c(2, 4, 6, 8), name="p.adjust") +
        xlab("Topic") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
              axis.text=element_text(size=20)))
dev.off()

setwd('../')
}  # end stm per cell type

