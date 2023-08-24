# Massively parallel in vivo Perturb-seq
Zheng et al 2023

## Perturb-seq clustering and QC

See PerturbSeqAnalysis and the associated documentation for the code used in upstream analysis of the Perturb-seq data (starting with CellRanger output through clustering and cell type identification).

## Perturb-seq downstream analysis
**aav_downstream.R**: Main driver script that calls each of the below  
**propeller.R**: Method to detect cell type proportion changes  
**mod.hidden.mult.R**: Multinomial version of [HiDDEN](https://github.com/tudaga/LabelCorrection/tree/main) to identify degrees of perturbation effect in single cells  
**run_sva_edger.R**: Identify DEGs using `edgeR` with 1 surrogate variable from `sva`  
**Enrich_FGSEA_new.R**: Identify enriched GO terms  
**RunSTM**: Group genes into modules by structural topic modeling (`stm`)  


## Other related analysis

### AAV serotype primary screen: AAV barcode bulk analysis
fastq_barcodemapping.py: Count number of barcodes for each AAV variant  
deseq2.R: Identify significantly enriched AAV variants


### AAV serotype secondary screen: AAV barcode-single cell analysis
XIN


### HypPB insertion site genome-wide analysis
annotations.R: annotate insertion sites
