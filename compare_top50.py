#!/usr/bin/env python

import argparse
from sets import Set
import numpy as np


parser = argparse.ArgumentParser(description='Compare top 50 experimental with top 50 predicted peptides (should be run in [PDB]/ folder)')
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
		scores[line.split()[0]] = [float(line.split()[1]), float(line.split()[header.index('bind_pep_alone')])]
	RND4.readline()
	for line in RND4:
		scores[line.split()[0]] = [float(line.split()[1]), float(line.split()[header.index('bind_pep_alone')])]
	PRE.close()
	RND4.close()

	# Store only round 4 reads in dictionary
	RND4 = open('round_4/all/score_table.txt', 'r')
	rnd4_scores = {}
	RND4.readline()
	for line in RND4:
		rnd4_scores[line.split()[0]] = [float(line.split()[1]), float(line.split()[header.index('bind_pep_alone')])]
	RND4.close()

	# Sort scoring
	sorted_scoring = sorted(scores.keys(), key = lambda k: scores[k][1])
	# Sort abundant
	sorted_abundant = sorted(rnd4_scores.keys(), key = lambda k: rnd4_scores[k][0], reverse=True)

	top50_scoring = Set(sorted_scoring[:50])
	top50_abundant = Set(sorted_abundant[:50])

	intersection = list(top50_scoring & top50_abundant)
	#intersection = list(top50_scoring & Set(sorted_abundant))
	#print len(intersection) 

	# top50 scoring peptides that are NOT in top50 abundant
	top50_abundant_comp = list(top50_scoring - top50_abundant)
	# top50 abundant peptides that are NOT in top50 scoring
	top50_scoring_comp = list(top50_abundant - top50_scoring)
	# Generate heatmaps
	intersection_h = heatmap(intersection)
	top50_abundant_comp_h = heatmap(top50_abundant_comp)
	top50_scoring_comp_h = heatmap(top50_scoring_comp)

	# Save heatmaps
	np.savetxt('compare_top50/intersection_50_bind_pep_alone_50_reads_3_peptides_submatrix.txt',
			intersection_h, delimiter='\t', fmt= '%1.4f')
	np.savetxt('compare_top50/scoring_only_50_bind_pep_alone_50_reads_47_peptides_submatrix.txt',
			top50_abundant_comp_h, delimiter='\t', fmt= '%1.4f')
	np.savetxt('compare_top50/reads_only_50_bind_pep_alone_50_reads_47_peptides_submatrix.txt',
			top50_scoring_comp_h, delimiter='\t', fmt= '%1.4f')

if __name__ == '__main__':
	main()