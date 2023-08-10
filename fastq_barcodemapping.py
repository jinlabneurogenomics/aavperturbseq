import sys
import os
import pandas as pd
import Levenshtein

def process(lines=None):
### Read fastq in and store as dictionary
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

def hamming(string1, string2):
    return sum(c1 != c2 for c1, c2 in zip(string1,string2))

def correct_primer(reference, string, reference2 = 'ACGG'):    
### reference: correct primer sequence
### string: read sequence
### reference2: start of sequence after first BC
    substring = string[:len(reference)]
    if hamming(reference, substring) == 1:
        if string[64:].startswith(reference2): #make sure not deletion
            return reference
    return None

def error_correct(barcode,poss_bar,bc_name,maxDist=2):
### Select closest barcode and barcode name to input barcode 
    dist=[Levenshtein.distance(barcode,bar_test) for bar_test in poss_bar]
    if min(dist)>maxDist:
        return '', ''
    correct_ind=dist.index(min(dist))
    correct_bar=poss_bar[correct_ind]
    correct_bcname=bc_name[correct_ind]
    return correct_bar, correct_bcname


barcode_file = 'AAV Capsid library -Final98XU_complete.xlsx'
barcodes = pd.read_excel(barcode_file).drop(28)
barcodes.columns = ['BCs','Name','BC1','spacer','BC2']

fastq_dir = 'fastqs/'
files = os.listdir(fastq_dir)

reads = []
discard = []
primer_seq1 = 'GACGAGTCGGATCTCCCTTTGGGCCGCCTCCCCGCATCGATACCGTCAATTG'
            # 'GACGAGTCGGATCTCCCTTTGGGCCGCCTCCCCGCATCGATACCGTCGACCT' 4978 reads most common error seq
# align = 'gacgagtcggatctccctttgggccgcctccccgcatcgataccgtCAATTGNNNNNNNNNNNNACGGAAATACGATGTCGGGANNNNNNNNNNNNGAGCtcgctgatcagcctcgactgtgccttctagttgccagccatctgttgtttgcccctcccccgtgccttccttgaccctggaaggtgccactcccactgtcctttcctaataaaatgaggaaattgcatcgc'.upper()
# primer_seq2 = 'ACGGAAATACGATGTCGGGA'
for fn in files:
    with open(fastq_dir + fn, 'r') as fh:
        lines = []
        for line in fh:
            lines.append(line.rstrip())
            if len(lines) == 4: #each fastq read is 4 lines
                record = process(lines)
                if record['sequence'].startswith(primer_seq1):
                    reads.append(record)
                else:
                    correction = correct_primer(primer_seq1, record['sequence'])
                    if correction is not None:
                        record['sequence'] = correction + record['sequence'][len(correction):]
                    record['hamming'] = hamming(primer_seq1, record['sequence'][:len(primer_seq1)])
                    if record['hamming'] == 0:
                        reads.append(record)
                    else:
                        discard.append(record)
                lines = []

counts = {bc: 0 for bc in barcodes['BCs']}
for seq in reads:
    seq['BC1_seq'] = seq['sequence'][52:64]
    seq['BC2_seq'] = seq['sequence'][84:96]
    seq['BC1_seq_correct'], seq['BC1_name'] = error_correct(seq['BC1_seq'], list(barcodes['BC1']), list(barcodes['BCs']))
    seq['BC2_seq_correct'], seq['BC2_name'] = error_correct(seq['BC2_seq'], list(barcodes['BC2']), list(barcodes['BCs']))
    # print(seq['BC1_name'],seq['BC2_name'])
    if seq['BC1_name'] == seq['BC2_name']:
        seq['BC_name'] = seq['BC1_name']
        if seq['BC1_name'] != '':
            counts[seq['BC_name']] += 1

counts_df = pd.DataFrame({'BC_name': counts.keys(), 'Count': counts.values()})
counts_df.to_csv('BC_counts.csv')