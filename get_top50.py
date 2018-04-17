#!/usr/bin/env python

import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Output top50 most abundant and top50 best scoring peptides for input TCR (should be run in [PDB]/ folder)')
parser.add_argument('-s', help='scoring function (ex. standard, bind, bind_pep_alone, atlas)', type=str, dest='s', required=True)
args = parser.parse_args()

# Constants
aa_dict = {"A":0,"C":1,"D":2,"E":3,"F":4,"G":5,"H":6,"I":7,"K":8,"L":9,"M":10,"N":11,"P":12,"Q":13,"R":14,"S":15,"T":16,"V":17,"W":18,"Y":19}


def heatmap(peptides):
	'''
	Create heatmap of amino acid frequencies
	'''
	length = len(peptides[0])
	# Initialize substitution matrix
	sub_matrix = np.zeros((20, length), dtype=float)
	for peptide in peptides:
		# Update substitution matrix
		for res_pos in range(0, length):
			aa_pos = aa_dict[peptide[res_pos]]
			sub_matrix[aa_pos,res_pos] += 1
	# Calculate percentage
	percent_matrix = sub_matrix/len(peptides)
	return percent_matrix


def main():

	# Store preselection and round 4 reads in dictionary
	PRE = open('preselection/all/score_table.txt', 'r')
	RND4 = open('round_4/all/score_table.txt', 'r')

	scores = {}
	header = PRE.readline().split()
	for line in PRE:
		scores[line.split()[0]] = [float(line.split()[1]), float(line.split()[header.index(args.s)])]
	RND4.readline()
	for line in RND4:
		scores[line.split()[0]] = [float(line.split()[1]), float(line.split()[header.index(args.s)])]
	PRE.close()
	RND4.close()

	# Store only round 4 reads in dictionary
	RND4 = open('round_4/all/score_table.txt', 'r')
	rnd4_scores = {}
	RND4.readline()
	for line in RND4:
		rnd4_scores[line.split()[0]] = [float(line.split()[1]), float(line.split()[header.index(args.s)])]
	RND4.close()

	# Sort scoring
	sorted_scoring = sorted(scores.keys(), key = lambda k: scores[k][1])
	# Sort abundant
	sorted_abundant = sorted(rnd4_scores.keys(), key = lambda k: rnd4_scores[k][0], reverse=True)
	
	# Make top scoring heatmaps

	bound = 50
	scoring_heatmap = heatmap(sorted_scoring[:bound])
	abundant_heatmap = heatmap(sorted_abundant[:bound])
	# Save heatmaps and peptides
	S_OUT = open('weblogo/' + args.s + '_score_' + str(bound) +'_peptides.fa', 'w')
	A_OUT = open('weblogo/abundant_' + str(bound) +'_peptides.fa', 'w')
	S_OUT.write('\n'.join(sorted_scoring[:bound]))
	A_OUT.write('\n'.join(sorted_abundant[:bound]))

	np.savetxt('weblogo/' + args.s + '_score_' + str(bound) + '_peptides_submatrix.txt',
		scoring_heatmap, delimiter='\t', fmt= '%1.4f')
	np.savetxt('weblogo/abundant_' + str(bound) + '_peptides_submatrix.txt', 
		abundant_heatmap, delimiter='\t', fmt= '%1.4f')

if __name__ == '__main__':
	main()
