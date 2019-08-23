#!/usr/bin/env python

import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Create heatmap matrices for nr data')
parser.add_argument('-f', help='peptide score file (ex. scores_below_-33.0)', type=str, dest='f', required=True) 
args = parser.parse_args()

# Constants
aa_dict = {"A":0,"C":1,"D":2,"E":3,"F":4,"G":5,"H":6,"I":7,"K":8,"L":9,"M":10,"N":11,"P":12,"Q":13,"R":14,"S":15,"T":16,"V":17,"W":18,"Y":19}

def heatmap(pep_file):
	'''
	Create heatmap of amino acid frequencies
	'''
	IN = open(pep_file, 'r')
	line=IN.readline()
	length = len(line.split()[0])
	IN.close()

	# Initialize substitution matrix
 	sub_matrix = np.zeros((20, length), dtype=float)
 	IN = open(pep_file, 'r')
 	total_peptides = 0
 	for line in IN:
 		peptide = line.split()[0]
		# Update substitution matrix
		for res_pos in range(0, length):
			aa_pos = aa_dict[peptide[res_pos]]
			sub_matrix[aa_pos,res_pos] += 1
		total_peptides += 1
		if total_peptides % 1000000 == 0:
			print 'On row: ' + str(total_peptides)
	IN.close()

 	# Calculate percentage
 	percent_matrix = sub_matrix/total_peptides
 	return percent_matrix		

def main():
	nr_heatmap = heatmap(args.f)
	np.savetxt(args.f + '_score_peptides_submatrix.txt', nr_heatmap, delimiter='\t', fmt= '%1.4f')

if __name__ == '__main__':
	main()