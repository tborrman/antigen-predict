#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Alphabet import IUPAC

def main():
	parser=argparse.ArgumentParser(description='Count redundancy from fastq file output by 3QIB_process_step1.py')
	parser.add_argument('-f', help= 'fastq file from step1 output', dest='fastq',type=str, required=True)

	args=parser.parse_args()
	fastq=args.fastq

	# Dictionary of seqs and counts
	out_dict = {}

	# Counter
	read_counter = 0

	for seq_record in SeqIO.parse(fastq, 'fastq', IUPAC.ambiguous_dna):
		read_counter += 1
		if read_counter % 10000 == 0:
			print 'Parsing seq: ' + str(read_counter)
		# Counted peptide position in codon seq by hand (37-75) 0-based inclusive
		codon_pep = seq_record.seq[37:76]
		pep_seq = str(codon_pep.translate())

		if 'X' in pep_seq or '*' in pep_seq:
			read_counter -= 1
			continue
		else:
			# Check if peptide has been found already
			if pep_seq in out_dict:
				out_dict[pep_seq] +=1
			# If not found add peptide as key
			else:
				out_dict[pep_seq] = 1

	# Print results to output file
	OUT = open('3QIB_redundant_count.txt' , 'w+')
	for key in out_dict.keys():
		OUT.write(key +  '\t' + str(out_dict[key]) + '\n')
	OUT.close()

	print 'Number of peptide processed: ' + str(read_counter) + '\n'
	





if __name__ == '__main__':
	main()
