outdir=$1 #outdir--use the name of the gene
reg=$2 #region to extract--in other words the coordinate of the gene I care about
dir=$3 #Directory with CellRanger output in it
cellstopert=$4
mkdir $outdir
echo Extract region
for num in {1..5}; do echo $num;bam=${dir}/possorted_genome_bam.bam;samtools view -Sb $bam $reg > $outdir/reads.${num}.bam;done
echo Combine multiple 10X channels
samtools cat $outdir/r*bam > $outdir/combined.bam
echo Sort
samtools sort $outdir/combined.bam > $outdir/combined.sort.bam
echo Index
samtools index $outdir/combined.sort.bam
cd $outdir
echo Split by perturbation
sinto filterbarcodes -b combined.sort.bam -c $cellstopert  ##splits cells from each perturbation into a different bam, labels.txt has the information about which cell goes with which perturbation
echo Index split bams
for f in *bam; do samtools index $f;done
