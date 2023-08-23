library(Seurat)
library(speckle)
library(limma)

######
# DONT USE THIS
######
# issues with duplicateCorrelation
# assumptions that are probably best to avoid (the use of duplicateCorrelation is probably the main reason that pseudocell has issues--it doesn't fit a standard mixed model), namely that each variable (in this case each cell type) has the same amount of batch effect (in some sense), which is not likely to be a reasonable assumption
# Also, it does not correctly estimate p-values (it does not take into account the uncertainty in the estimation of the random effects when calculating p-values, so leads to smaller p-values than it would get normally).
# run propeller with random effect
RunPropellerRand <- function(props, batch, treat, out.prefix="propeller_result") {
# design matrix specification: (1) group information is always first and (2) there is NO intercept
mm.treat = model.matrix(~treat)
dupcor = duplicateCorrelation(props$TransformedProps, design=mm.treat, block=batch)
fit1 = lmFit(props$TransformedProps, design=mm.treat, block=batch, correlation=dupcor$consensus)
fit1 <- eBayes(fit1)
print(summary(decideTests(fit1)))
res=topTable(fit1, number=Inf)
print(res)
write.table(res, paste0(out.prefix, ".tsv"), sep='\t', quote=F, row.names=T, col.names=NA)
return(res)
}



# Add batch as fixed effect
RunPropeller <- function(props, batch, treat, print.res=T, out.prefix="propeller_result") {
meta = data.frame(batch, treat)
meta$treat = relevel(as.factor(meta$treat), ref='GFP')
meta$batch = as.factor(meta$batch)
design <- model.matrix(~ 0 + treat + batch)
#print(design)
contrast.cmd = paste0("makeContrasts(", paste0('treat',levels(meta$treat)[2], '-treat',levels(meta$treat)[1]), ", levels=design)")
my.contrast = eval(parse(text=contrast.cmd))
#print(my.contrast)
# robust = TRUE
# robust empirical Bayes shrinkage of the variances is performed which mitigates the effects of outlying observations
# trend = FALSE
# as we don't expect a mean-variance trend after performing our variance stabilising transformation. 
# There may also be an error when `trend` is set to TRUE because there are often not enough data points to estimate the trend.
res = propeller.ttest(prop.list=props, design=design, contrast=my.contrast, robust=TRUE, trend=FALSE, sort=TRUE)
if (print.res) print(res)
write.table(res, paste0(out.prefix, ".tsv"), sep='\t', quote=F, row.names=T, col.names=NA)
return(res)
}



library(dplyr)
library(ggplot2)
# sample column
# treat column
# cell_type column
MyBar <- function(df, sig.index) {
xface = rep('plain',nlevels(df$clusters))
xface[sig.index] = 'bold'
xsize = rep(8, nlevels(df$clusters))
xsize[sig.index] = 12
p=df %>%
  ggplot(aes(x=clusters, y=mean, fill=guide))+
  geom_bar(stat = "identity", position = position_dodge())+
  #geom_text(aes(label = paste(round(mean,2),"%")), position = position_dodge(0.9), vjust =-2, size = 3)+
  geom_errorbar(aes(ymin=mean-ic, ymax=mean+ic), width=0.3, position=position_dodge(0.9))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(face=xface, size=xsize))
  #scale_y_continuous(limits = c(-2,60))
return(p)
}

