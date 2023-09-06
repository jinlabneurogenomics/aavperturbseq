library(Seurat)
library(speckle)
library(limma)

my.propeller.ttest <- function(prop.list=prop.list, design=design,
                            contrasts=contrasts, robust=robust, trend=trend,
                            sort=sort)
{
    prop.trans <- prop.list$TransformedProps
    prop <- prop.list$Proportions

    fit <- lmFit(prop.trans, design)
    fit.cont <- contrasts.fit(fit, contrasts=contrasts)
    fit.cont <- eBayes(fit.cont, robust=robust, trend=trend)

    # Get mean cell type proportions and relative risk for output
    # If no confounding variable included in design matrix
    if(length(contrasts)==2){
        fit.prop <- lmFit(prop, design)
        z <- apply(fit.prop$coefficients, 1, function(x) x^contrasts)
        RR <- apply(z, 2, prod)
    }
    # If confounding variables included in design matrix exclude them
    else{
        new.des <- design[,contrasts!=0]
        fit.prop <- lmFit(prop,new.des)
        new.cont <- contrasts[contrasts!=0]
        z <- apply(fit.prop$coefficients, 1, function(x) x^new.cont)
        RR <- apply(z, 2, prod)
    }

    fdr <- p.adjust(fit.cont$p.value, method="BH")

    out <- data.frame(PropMean=fit.prop$coefficients, PropRatio=RR, limmaCoef=fit.cont$coefficients,
                    Tstatistic=fit.cont$t[,1], P.Value=fit.cont$p.value[,1],
                    FDR=fdr)
    if(sort){
        o <- order(out$P.Value)
        out[o,]
    }
    else out
}




# Add batch as fixed effect
RunPropeller <- function(props, meta, print.res=T, out.prefix="propeller_result", form="~ 0 + treat + batch") {
design <- model.matrix(as.formula(form))
#print(design)
contrast.cmd = paste0("makeContrasts(", paste0('treat',levels(meta$treat)[2], '-treat',levels(meta$treat)[1]), ", levels=design)")
print(contrast.cmd)
my.contrast = eval(parse(text=contrast.cmd))
#print(my.contrast)
# robust = TRUE
# robust empirical Bayes shrinkage of the variances is performed which mitigates the effects of outlying observations
# trend = FALSE
# as we don't expect a mean-variance trend after performing our variance stabilising transformation. 
# There may also be an error when `trend` is set to TRUE because there are often not enough data points to estimate the trend.
res = my.propeller.ttest(prop.list=props, design=design, contrasts=my.contrast, robust=TRUE, trend=FALSE, sort=TRUE)
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

