# Massively parallel in vivo Perturb-seq
Zheng et al 2023

Leveraging AAV's versatile tropism and labeling capacity, we expanded the scale of in vivo CRISPR screen with single-cell transcriptomic phenotyping across embryonic to adult brains and peripheral nervous systems. Through extensive tests of 86 AAV serotypes, combined with transposon systems, we substantially amplified labeling and accelerated in vivo gene delivery from weeks to days. We performed a proof-of-principle in utero screen and identified pleiotropic effects of Foxg1, featured by its tight regulation of distinct networks essential for cell fate specification of Layer 6 corticothalamic neurons. Notably, our platform can label >6% of cerebral cells, surpassing the current state-of-the-art efficacy at <0.1% (by lentivirus), to achieve analysis of over 30,000 cells in one experiment, thus enabling massively parallel in vivo Perturb-seq. Compatible with various phenotypic measurements (single-cell or spatial multi-omics), it presents a flexible approach to interrogate gene function across cell types in vivo, translating gene variants to their causal functions.

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

### HypPB insertion site genome-wide analysis
annotations.R: annotate insertion sites
