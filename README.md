# Massively parallel in vivo Perturb-seq
Zheng et al 2023

## Perturb-seq related analysis

FILL ME IN HERE



### AAV serotype secondary screen: AAV barcode-single cell analysis
XIN

## Non-single cell RNAseq analysis

### AAV serotype primary screen: AAV barcode bulk analysis
fastq_barcodemapping.py: Count number of barcodes for each AAV variant  
deseq2.R: Identify significantly enriched AAV variants

### HypPB insertion site genome-wide analysis
annotations.R: annotate insertion sites


## Downstream analysis
**aav_downstream.R**: Main driver script that calls each of the below
**propeller.R**: Method to detect cell type proportion changes
**mod.hidden.mult.R**: Multinomial version of [HiDDEN](https://github.com/tudaga/LabelCorrection/tree/main) to identify degrees of perturbation effect in single cells
**run_sva_edger.R**: Identify DEGs using `edgeR` with 1 surrogate variable from `sva`
**Enrich_FGSEA_new.R**: Identify enriched GO terms
**RunSTM**: Group genes into modules by structural topic modeling (`stm`)
