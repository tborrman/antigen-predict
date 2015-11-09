#!/usr/bin/env python
import numpy as np
import argparse
parser = argparse.ArgumentParser(description='Identify which substitutions have frequency 2-fold higher than expected')
parser.add_argument('-pdb', help='PDB id (ex. 3QIB)', type=str, required=True, dest='pdb')
parser.add_argument('-wt', help='wildtype peptide seq (ex. ADLIAYLKQATKG)', type=str, required=True, dest='wt')
args = parser.parse_args()

# Amino acid frequencies based on codon degeneracy using NNK
freq = np.array([2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 2.0, 1.0, 3.0, 3.0, 2.0, 2.0, 1.0, 1.0]) 
aa_list = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
NNK_probs = freq/np.sum(freq)

print '************************************************************************'
print 'Substiutions with frequency 2-fold higher than expected in both heatmaps'
print '************************************************************************'
# Open heatmap files
abundant_heatmap = np.loadtxt(args.pdb + '/heatmaps/abundant_50_peptides_submatrix.txt')
score_heatmap = np.loadtxt(args.pdb + '/heatmaps/bind_pep_alone_score_50_peptides_submatrix.txt')
for aa_pos in range(len(aa_list)):
	for pep_pos in range(len(args.wt)):
		# Skip if we're on wildype residue
		if args.wt[pep_pos] == aa_list[aa_pos]:
			continue	
		else:
			# Check if greater than 2-fold change in both heatmaps
			if (abundant_heatmap[aa_pos, pep_pos] > NNK_probs[aa_pos]*2) and (score_heatmap[aa_pos, pep_pos] > NNK_probs[aa_pos]*2):
				print 'Amino acid: ' + aa_list[aa_pos] +'\tPeptide position: ' + str(pep_pos + 1)
				

