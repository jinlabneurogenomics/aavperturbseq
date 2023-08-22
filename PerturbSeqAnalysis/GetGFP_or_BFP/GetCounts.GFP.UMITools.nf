params.fa="/stanley/levin_dr_storage/ssimmons/XinAAVPert/Bams/RunGFPBFP/seq.fa" //fasta with seqs to use
params.bam //input bam from CellRanger
//params.t2g //input t2g.txt
params.outdir="output" // The final output directory

workflow{
    unmap=GetUnmapped(params.bam)
    umi=ProcessUMI(unmap)
    mapped=MapToFasta(umi,params.fa)
    dedup=GeneAnnUMI(mapped)
    quant=Quantify(dedup)
}

process GetUnmapped //Extracts unmapped reads from CellRanger bam and makes fastq
{
    input:
    path "input.bam"

    output:
    path "output"

    '''
    echo Get unmapped
    samtools view -f 4 -b input.bam > unmap.bam
    echo Make fastq
    cellranger bamtofastq unmap.bam output
    cat output/*/*R1_001.fastq.gz > output/umap.R1.fq.gz
    cat output/*/*R2_001.fastq.gz > output/umap.R2.fq.gz
    rm unmap.bam
    rm output/*/*gz
    '''
}

process ProcessUMI //Extracts UMI and CBC from read 1 and adds to read 2
{
    input:
    path "output"

    output:
    path "umi.fq.gz"

    '''
    umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --stdin output/umap.R1.fq.gz --read2-in output/umap.R2.fq.gz --read2-stdout --stdout umi.fq.gz
    '''

}

process MapToFasta //Maps reads to a fasta using minimap2
{

    input:
    path "umi.fq.gz"
    path "seq.fa"

    output:
    path "mapped.bam"

    '''
    minimap2 -ax sr seq.fa umi.fq.gz | samtools view -Sb > mapped.bam
    '''
}

process GeneAnnUMI //Annotates bam file with Gene (in this case chromosome) and sorts it
{
    input: 
    path "mapped.bam"

    output:
    path "ann.sort.bam"

    '''
    samtools view -H mapped.bam > header.txt
    samtools view -F 4 mapped.bam | awk '{print \$0"\tXT:Z:"\$3}' | cat header.txt - | samtools view -Sb > ann.bam
    samtools sort ann.bam > ann.sort.bam
    rm ann.bam
    rm header.txt
    '''
}

process Quantify //Counts the number of UMIs per cell for each gene
{
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path "ann.sort.bam"

    output:
    path "counts.tsv.gz"

    '''
    samtools index ann.sort.bam
    umi_tools count --per-gene --gene-tag=XT --per-cell -I ann.sort.bam -S counts.tsv.gz
    '''
}
