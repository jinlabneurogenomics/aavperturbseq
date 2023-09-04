library(tidyverse)
library(magrittr)
library(ChIPseeker)

library(TxDb.Mmusculus.UCSC.mm39.refGene)
txdb <- TxDb.Mmusculus.UCSC.mm39.refGene

setwd('/filepath/to/working/dir')
reads <- readxl::read_xlsx('junction_data.xlsx', sheet = 1) 

bed <- cbind(reads$chr1_name, reads$pos1, reads$pos1)
write.table(bed, file = 'junctions.bed', sep = '\t', col.names = FALSE, row.names = FALSE)

peak <- readPeakFile('junctions.bed')
peakAnno <- annotatePeak('junctions.bed', tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")

pdf('annoPie.pdf')
plotAnnoPie(peakAnno)
dev.off()

pdf('annoVenn.pdf')
vennpie(peakAnno)
dev.off()

anno_df <- as.data.frame(peakAnno)
colnames(anno_df) <- c('chr1_name','start', 'pos1', colnames(anno_df)[4:length(colnames(anno_df))])
df <- merge(reads,anno_df, by = c('chr1_name','pos1'), all = TRUE)
write.csv(df, 'chipseeker.csv')

anno <- df$annotation
anno[str_detect(anno, 'Intron')] <- 'Intron'
table(anno) %>% as.data.frame() %>% write.csv('junction_counts.csv')
