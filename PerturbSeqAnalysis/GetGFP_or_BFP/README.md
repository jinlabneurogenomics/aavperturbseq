## Counting genes in 10X

In many cases when running 10X experiments we want to count the expression of transgenes and the like that are not part of the standard reference (such as GFP, BFP, etc). It is possible to add such sequences to the 10X reference. However, an alternative is to instead extract these sequences from the premapped 10X bam file. There are downsides to this approach (it is possible for such reads to be mis-aligned to a region of the reference genome, for example), but there are advantages as well (easier and faster than rerunning the entire 10X pipeline if you don't think to include the sequence ahead of time, etc). As such, we create a simple nextflow pipeline to extract this information.

# Steps in pipeline

The pipeline is fairly simple. The pipeline takes in a CellRanger bam file and a fasta file with the sequences we are looking for. First, it extracts unmapped reads from the CellRanger bam file. It changes this into paired fastq file, extracts the UMI and CBC with umitools. The resulting fastq file is mapped to the reference of interest with minimap2. Finally, umitools is used to count the number of UMI in each cell mapping to each contig in the reference fasta. 

# How to run

Need to have umitools, samtools, minimap2, and cellranger installed and on your path (will add script to do this with conda or similiar). Will also need to download nextflow (see the nextflow website). Then it is as simple as running:

```
nextflow $pathtocode/GetCounts.GFP.UMITools.nf [args]
```

where pathtocode is the path to the directory containing GetCounts.GFP.UMITools.nf. The arguments one can give are:

`--fa:` The fasta with the sequences you are interested in looking for. One entry per sequence. Required.

`--bam:` The CellRanger (or similiar) bam from the single cell experiment. Required.

`--outdir:` The name of the output directory. Set to output by default.

