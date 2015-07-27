#!/usr/bin/env python

import argparse
import numpy as np
import math

def main():
	parser=argparse.ArgumentParser(description='build substitution matrix based on rosetta3.5 scores')
	parser.add_argument('-p', help='number of parts output was split into', type=int, dest='parts', required=True)
	parser.add_argument('-l', help='wt peptide sequence length', type=int, dest='length', required=True)
	parser.add_argument('-pMHC', help='calculate dG for pMHC complex', dest='pMHC', action='store_true')

	args=parser.parse_args()
	parts=args.parts
	length=args.length
	pMHC=args.pMHC

	if pMHC:
		TCR_part_file = 'TCR_score.sc'
	else:
		TCR_part_file = 'TCR_MHC_score.sc'

	aa_dict = {"A":0,"C":1,"D":2,"E":3,"F":4,"G":5,"H":6,"I":7,"K":8,"L":9,"M":10,"N":11,"P":12,"Q":13,"R":14,"S":15,"T":16,"V":17,"W":18,"Y":19}

	# Calculate list of dG_bind scores
	dG_bind = []
	TCR_PART = open(TCR_part_file, 'r')
	TCR_PART_score = float(TCR_PART.readline().split()[1])
	
	for part_num in range(1, parts + 1):
		PEP= open('pep_score_' + str(part_num) + '.sc', 'r')
		PEP_list = PEP.readlines()
		COM= open('com_score_' + str(part_num) + '.sc', 'r')
		COM_list = COM.readlines()
		for i in range(len(PEP_list)):
			PEP_score = float(PEP_list[i].split()[1])
			COM_score = float(COM_list[i].split()[1])
			dG_score = COM_score - (TCR_PART_score + PEP_score)
			dG_bind.append(dG_score)
		PEP.close()
		COM.close()
	min_score = min(dG_bind)
	max_score = max(dG_bind)

	dG_cutoffs = np.arange(math.ceil(max_score), math.ceil(min_score), -0.5)

	# Build substitution matrices
	for dG_cutoff in dG_cutoffs:
		bind_counter = 0
		pep_counter = 0
		# Initialize substitution matrix
		sub_matrix = np.zeros((20, length), dtype=float)
		# Open peptide sequence file
		for part_num in range(1, parts +1):
			seq_fh = open('redundant_out_' + str(part_num) + '.txt', 'r')
			seq_lines = seq_fh.readlines()
			for line in seq_lines:
				seq_line = line.split()
				seq = seq_line[0]
				score = dG_bind[pep_counter]
				pep_counter += 1
				if score <= dG_cutoff:
					bind_counter +=1
					# Update substitution matrix
					for res_pos in range(0, len(seq)):
						aa_pos = aa_dict[seq[res_pos]]
						sub_matrix[aa_pos,res_pos] += 1
		# Get percentages
		percent_matrix = sub_matrix/bind_counter

		# Output to file
		np.savetxt('sub_matrix_ros3.5_' + str(dG_cutoff) + '.txt', percent_matrix, delimiter='\t', fmt= '%1.4f')
	print 'bind_counter\n'
	print bind_counter
	print 'pep_counter\n'
	print pep_counter

if __name__ == '__main__':
	main()