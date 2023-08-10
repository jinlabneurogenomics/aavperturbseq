library(tidyverse)
library(xlsx)
setwd('/filepath/to/working/dir')

bcs <- read.csv('bc_list.csv', row.names = 1)

meta <- read.xlsx('BC_counts_summary.xlsx', sheetIndex=2)
meta['expt'] <- paste(meta[['experiment']],meta[['time']],meta[['condition']], sep = '_')

mat <- t(subset(read.xlsx('BC_counts_summary.xlsx', sheetIndex=1, colClasses = "numeric")[-47,], select=3:101))
colnames(mat) <- meta[['sample']]
rownames(mat) <- bcs[rownames(mat),1]
mat <- mat[rowSums(mat) > 10,]

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = meta,
                              design = ~ expt)
dds$expt <- relevel(dds$expt, ref = "lib_0_lib")

dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
res %>% as.data.frame() %>% arrange(-log2FoldChange)

res_df <- data.frame()
samples <- c('HT_24_sort','HT_48_sort','mm_24_sort','mm_48_sort')

library("EnhancedVolcano")
for(sample in samples){
  res <- results(dds, alpha = 0.05, contrast = c("expt",sample,"lib_0_lib")) %>% 
    as.data.frame() %>% rownames_to_column('BC')
  res['sample'] <- sample
  res_df <- rbind(res_df,res)
  EnhancedVolcano(res,
                  lab = res$BC,
                  x = 'log2FoldChange',
                  y = 'padj',
                  xlim = c(-5,5),
                  pCutoff = 0.05,
                  title = sample)
  ggsave(paste0('figures/volcanoes/',sample,'_volcano.pdf'), width = 5, height = 5)
}
write.csv(res_df, file = 'sampleDE.csv')

res_df2 <- data.frame()
contrasts <- list(c('HT_24_sort','HT_48_sort'),c('mm_24_sort','mm_48_sort'),c('mm_24_sort','HT_24_sort'),c('mm_48_sort','HT_48_sort'))
for (contrast in contrasts){
  print(contrast)
  res <- results(dds, alpha = 0.05, contrast = c("expt",contrast[2],contrast[1])) %>% 
    as.data.frame() %>% rownames_to_column('BC')
  res['sample'] <- paste(contrast[2],contrast[1],sep='/')
  res_df2 <- rbind(res_df2,res)
  EnhancedVolcano(res,
                  lab = res$BC,
                  x = 'log2FoldChange',
                  y = 'padj',
                  xlim = c(-5,5),
                  pCutoff = 0.05,
                  title = sample)
  ggsave(paste0('figures/volcanoes/',contrast[2],'_vs_',contrast[1],'_volcano.pdf'), width = 5, height = 5)
}
write.csv(res_df2, file = 'timeConditionDE.csv')
save(dds, res_df, res_df2, file = 'sampleDE.Robj')

load('sampleDE.Robj')
library("pheatmap")
# select <- order(rowMeans(counts(dds,normalized=TRUE)),
#                 decreasing=TRUE)
df <- as.data.frame(colData(dds)[,c("time","condition","experiment")])

vsd <- varianceStabilizingTransformation(dds)
pheatmap(assay(vsd)[,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, filename = 'figures/heatmap.pdf', width = 10, height = 15)

pheatmap(assay(vsd)[,sort(colData(vsd)[colData(vsd)$time %in% c(0,24) & colData(vsd)$experiment %in% c('lib','HT') & colData(vsd)$condition != 'unsort','sample'])], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, filename = 'figures/heatmaps/heatmap_HT24h.pdf', width = 6, height = 12)
pheatmap(assay(vsd)[,sort(colData(vsd)[colData(vsd)$time %in% c(0,48) & colData(vsd)$experiment %in% c('lib','HT') & colData(vsd)$condition != 'unsort','sample'])], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, filename = 'figures/heatmaps/heatmap_HT48h.pdf', width = 6, height = 12)
pheatmap(assay(vsd)[,sort(colData(vsd)[colData(vsd)$time %in% c(0,24) & colData(vsd)$experiment %in% c('lib','mm') & colData(vsd)$condition != 'unsort','sample'])], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, filename = 'figures/heatmaps/heatmap_mm24h.pdf', width = 6, height = 12)
pheatmap(assay(vsd)[,sort(colData(vsd)[colData(vsd)$time %in% c(0,48) & colData(vsd)$experiment %in% c('lib','mm') & colData(vsd)$condition != 'unsort','sample'])], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, filename = 'figures/heatmaps/heatmap_mm48h.pdf', width = 6, height = 12)
pheatmap(assay(vsd)[,sort(colData(vsd)[colData(vsd)$time %in% c(0,24) & colData(vsd)$condition != 'unsort','sample'])], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, filename = 'figures/heatmaps/heatmap_24h.pdf', width = 8, height = 12)
pheatmap(assay(vsd)[,sort(colData(vsd)[colData(vsd)$time %in% c(0,48) & colData(vsd)$condition != 'unsort','sample'])], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df, filename = 'figures/heatmaps/heatmap_48h.pdf', width = 8, height = 12)

colData(vsd)$by_time <- paste(colData(vsd)[['time']],colData(vsd)[['experiment']],colData(vsd)[['condition']], sep = '_')
colnames(vsd) <- colData(vsd)$by_time
pheatmap(assay(vsd)[,order(colData(vsd)$time, colData(vsd)$sample)], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=as.data.frame(colData(vsd)[,c("time","condition","experiment")]), 
         filename = 'figures/heatmaps/heatmap_bytime.pdf', width = 10, height = 15)
sub_ht <- vsd[,colData(vsd)$experiment %in% c('lib','HT') & colData(vsd)$condition != 'unsort']
pheatmap(assay(sub_ht)[,order(colData(sub_ht)$time, colData(sub_ht)$sample)], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=as.data.frame(colData(sub_ht)[,c("time","condition","experiment")]), 
         filename = 'figures/heatmaps/heatmap_HTbytime.pdf', width = 10, height = 15)
sub_mm <- vsd[,colData(vsd)$experiment %in% c('lib','mm') & colData(vsd)$condition != 'unsort']
pheatmap(assay(sub_mm)[,order(colData(sub_mm)$time, colData(sub_mm)$sample)], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=as.data.frame(colData(sub_mm)[,c("time","condition","experiment")]), 
         filename = 'figures/heatmaps/heatmap_mmbytime.pdf', width = 10, height = 15)

vsd <- vsd[,colData(vsd)$condition != 'unsort']

library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$sample
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, 
         filename = 'figures/sampleDistHeatmap.pdf', width = 5, height = 5)

pcaData <- plotPCA(vsd, intgroup=c("time", "experiment", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=as.factor(condition), shape=paste(experiment,time,sep = '_'))) +
  geom_point(size=2) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_classic()
ggsave('figures/pca.pdf', width = 8, height = 8)
