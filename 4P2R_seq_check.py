#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import numpy as np

def main():
	parser=argparse.ArgumentParser(description='Build substitution matrix for 4P2R selections output by 4P2R_process_step1.py to check against Birnbaum results')
	parser.add_argument('-f', help='fastq file from step1 output', dest='fastq', required=True, type=str)

	args=parser.parse_args()
	fastq=args.fastq

	read_counter = 0
	aa_dict = {"A":0,"C":1,"D":2,"E":3,"F":4,"G":5,"H":6,"I":7,"K":8,"L":9,"M":10,"N":11,"P":12,"Q":13,"R":14,"S":15,"T":16,"V":17,"W":18,"Y":19}
	sub_matrix = np.zeros((20,13), dtype=float)

	for seq_record in SeqIO.parse(fastq, 'fastq', IUPAC.ambiguous_dna):
		read_counter += 1
		if read_counter % 10000 == 0:
			print 'Parsing read: ' + str(read_counter)

		# Counted peptide position in codon seq by hand (37-75) 0-based inclusive
		codon_pep = seq_record.seq[37:76]
		pep_seq = str(codon_pep.translate())

		# Remove peptides with unknown amino acids or stop codons
		if 'X' in pep_seq or '*' in pep_seq:
			read_counter -= 1
			continue
		else:	
			for res_pos in range(len(pep_seq)):
				aa_pos = aa_dict[pep_seq[res_pos]]
				sub_matrix[aa_pos, res_pos] += 1

	percent_matrix = sub_matrix/read_counter

	print 'Total reads processed: ' + str(read_counter)
	
	np.savetxt('4P2R_sub_matrix_seq_check.txt', percent_matrix, delimiter='\t', fmt='%1.4f')


if __name__ == '__main__':
	main()