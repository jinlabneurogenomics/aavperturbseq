## Extract Indels

This code (CodeToUse.sh) extracts a region of interest from the 10X bam file (usually the region targetted by a guide), splits the bam up by perturbation, then indexes the results. This requires sinto and samtools to be installed and on the PATH. Takes 4 arguments (must be given in order):

`outdir:` The name of the output directory (will create if doesn't exist).

`reg:` The region to target (similiar syntax to samtools).

`dir:` The directory containing the 10X channels of interest in it. At the moment the code assumes 5 subdirectories (Ch1, ..., Ch5) that each contain the CellRanger output for a different channel, but can easily be modified.

`cellstopert:` A tab seperated file with two columns, the first being cell barcodes, the second being the guide each cell barcode is assigned to.

Can then run `source CodeToUse.sh $outdir $reg $dir $cellstopert` to get one bam file per perturbation in the outdir.

## Count Indels

GetInsertDelete.py is used to count the number of indels per cell in a particular region in a particular bam file. Designed to be run after the CodeToUse.sh step, but can be run on any bam file. For a given bam file and region can run:

```
import GetInsertDelete as GID
bamfile="/path/to/bam"
chrom="chromosome_of_region"
start=start_position_region
end=end_position_region
dat=GID.GetInsertDel(bamfile,chrom,start,end)
```

Can also run on all bams in a directory with GetInsertDel_directory (pass dir, a string pointing to a directory, instead of a bamfile, the dir string can contain regular expressions allowing multiple directories at once), and on multiple regions with GetInsertDel_regions (pass it the directory and a list of regions, where each region is a list with 4 entries: chromosomse of region, start of region, end position of region, and name of region (the name is up to you)). More details will be added. 





