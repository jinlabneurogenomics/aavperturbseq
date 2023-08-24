## Extract Indels

This code (CodeToUse.sh) extracts a region of interest from the 10X bam file (usually the region targetted by a guide), splits the bam up by perturbation, then indexes the results. This requires sinto and samtools to be installed and on the PATH. Takes 4 arguments (must be given in order):

`outdir:` The name of the output directory (will create if doesn't exist).

`reg:` The region to target (similiar syntax to samtools).

`dir:` The directory containing the 10X channels of interest in it. At the moment the code assumes 5 subdirectories (Ch1, ..., Ch5) that each contain the CellRanger output for a different channel, but can easily be modified.

`cellstopert:` A tab seperated file with two columns, the first being cell barcodes, the second being the guide each cell barcode is assigned to.

Can then run `source CodeToUse.sh $outdir $reg $dir $cellstopert` to get one bam file per perturbation in the outdir.

##Count Indels

GetInsertDelete.py is used to count the number of indels in a particular region in a particular bam file. Designed to be run after the CodeToUse.sh step, but can be run on any bam file. Will add details.

