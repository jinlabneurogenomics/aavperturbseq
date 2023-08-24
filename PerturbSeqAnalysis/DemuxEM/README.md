## DemuxEM on Perturb-seq guide RNA

DemuxEM is a tool built to assign cells/nuclei to samples using hashing barcodes. The same underlying math, however, should be applicable to assigning cells/nuclei to perturbations using guide barcodes in a Perturb-seq or Crop-seq experiment. The code in this directory is meant to enable this.

# How to prepare computational environment

To run this pipeline, need python installed, with the packages pandas and pegasusio installed. In addition, need demuxEM to be installed (the demuxEM binary should either be placed in the bin directory in this repo or needs a pointer to the demuxEM command). Also need to have nextflow downloaded. Finally, need to have the context of this directory downloaded to somewhere on your local machine (call the directory codedir)

# Running the code

To run the code run:

```
nextflow $codedir/RunDemux.nf [args]
```

The arguments you can use are:
`--inDir:` The name of a CellRanges outs directory, assumes include Crispr counts in analysis.

`--demuxCode:` A pointer to the demuxEM binary. By default will look in the bin subdirectory of this directory.

`--outdir:` The name of the out directory, by default is a directory DemuxEM in the CellRanger outs directory.

# Output

The output will be a file named Results.csv in the specified output directory. This is a table with the assignment for each cell.

