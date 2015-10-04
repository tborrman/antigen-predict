#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import re

def main():
	parser = argparse.ArgumentParser(description= 'Deconvolute 3QIU seq data into rounds of selection from specific library')
	parser.add_argument('-f', help='fastq file (must be only single reads or single mates from paired end reads)', dest='fastq', type=str, required=True)
	parser.add_argument('-b', help=' 6 nucleotide dna barcode for selection round and library', type=str, dest='barcode', required=True)

	args=parser.parse_args()
	fastq=args.fastq
	barcode=args.barcode

	sel = 0 
	tot = 0

	sel_seq = barcode + 'CTGTTATTGCTAGCGTTTTAGCA'
	
	OUT = open('3QIU_process_step1.fastq', 'w+')
	for seq_record in SeqIO.parse(fastq, 'fastq', IUPAC.ambiguous_dna):
		tot += 1
		if tot % 100000 == 0:
			print 'On read:' + str(tot)
		if re.search(sel_seq, str(seq_record.seq)):
			sel += 1
			OUT.write(seq_record.format('fastq'))
	print 'Total reads processed: ' + str(tot)
	print 'Total reads selected: ' + str(sel)
	print 'Percentage: ' + str((float(sel)/float(tot)) * 100) + '%'
	OUT.close()


if __name__ == '__main__':
	main()
	
