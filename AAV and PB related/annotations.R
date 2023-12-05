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
peak <- peak[!duplicated(peak)]
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

## Plot background distribution (
background <- readPeakFile('insertSequences.bed') # BED file of all TTAATTAA seqences in genome
backgroundAnno <- annotatePeak(background, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Mm.eg.db")

pdf('backgroundPie.pdf')
plotAnnoPie(backgroundAnno)
dev.off()

pdf('backgroundVenn.pdf')
vennpie(backgroudAnno)
dev.off()

anno_df <- as.data.frame(backgroundAnno)
anno <- anno_df$annotation
anno[str_detect(anno, 'Intron')] <- 'Intron'
anno[str_detect(anno, 'Exon')] <- 'Exon'
table(anno) %>% as.data.frame() %>% write.csv('background_counts.csv')

## Comparison between insertions and background
df <- read.csv('background_counts.csv')
df2 <- read.csv('insert_counts.csv')
df <- full_join(df, df2, by = 'anno', suffix = c("_background", "_insert"))
df %>% write.csv('joined_counts.csv')

back_sum <- sum(df$Freq_background)
insert_sum <- sum(df$Freq_insert, na.rm = TRUE)
df %<>% mutate(Percent_background = Freq_background / back_sum, 
               Percent_insert = Freq_insert / insert_sum) %>%
  dplyr::select(!starts_with("X_")) %>%
  pivot_longer(!anno)

df %>% 
  filter(name %in% c('Percent_background', 'Percent_insert'), anno != 'Downstream (<=300bp)') %>%
  ggplot() +
  geom_bar(mapping = aes(x = anno, y = value, fill = name), stat = 'identity', position = 'dodge') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('fraction of insertions')
ggsave('insert_vs_background.pdf')

df2 %>% ggplot() +
  geom_bar(mapping = aes(x = anno, y = Freq), stat = 'identity') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave('insert.pdf')
