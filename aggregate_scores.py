#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Aggregate scoring results output by score.py with read counts')
parser.add_argument('-p', help='number of parts output was split into', dest='p', type=int, required=True)
args=parser.parse_args()


def main():
	# Output file
	OUT = open('score_table.txt', 'w')
	header = 'peptide\treads\tstandard\tbind\tbind_pep_alone\tatlas\n'
	OUT.write(header)
	# Invariable scores
	STD_TCR = open('standard_score/TCR_score.sc', 'r')
	std_TCR_score = float(STD_TCR.readlines()[1].split()[1])
	STD_TCR.close()
	STD_TCR_MHC = open('standard_score/TCR_MHC_score.sc', 'r')
	std_TCR_MHC_score = float(STD_TCR_MHC.readlines()[1].split()[1])
	STD_TCR_MHC.close()
	ATLAS_TCR = open('atlas_score/TCR_score.sc', 'r')
	atlas_TCR_score = float(ATLAS_TCR.readlines()[1].split()[1])
	ATLAS_TCR.close()
	
	for part in range(1, args.p + 1):
		# Open read counts (redundant out files)
		READ_CNT=open('redundant_out_' + str(part) + '.txt', 'r')
		read_counts = READ_CNT.readlines()
		# Open score files
		# STANDARD scores
		STD_COM = open('standard_score/com_score_' + str(part) + '.sc', 'r')
		std_com_scores = STD_COM.readlines()
		STD_PEP = open('standard_score/peptide_score_' + str(part) + '.sc', 'r')
		std_pep_scores = STD_PEP.readlines()
		STD_PMHC = open('standard_score/pMHC_score_' + str(part) + '.sc', 'r')
		std_pMHC_scores = STD_PMHC.readlines()
		# ATLAS scores
		ATLAS_COM = open('atlas_score/com_score_' + str(part) + '.sc', 'r')
		atlas_com_scores = ATLAS_COM.readlines()
		ATLAS_PMHC = open('atlas_score/pMHC_score_' + str(part) + '.sc', 'r')
		atlas_pMHC_scores = ATLAS_PMHC.readlines()

		# Length checks
		list_len = len(read_counts)
		if not all(len(x) == list_len for x in [std_com_scores, std_pep_scores, std_pMHC_scores, atlas_com_scores, atlas_pMHC_scores]):
			print 'ERROR: different length of score files'
			quit()

		for i in range(list_len):
			line = read_counts[i].split()
			# STANDARD
			line.append(float(std_com_scores[i].split()[1]))
			# BIND
			COM = float(std_com_scores[i].split()[1])
			TCR = std_TCR_score
			pMHC = float(std_pMHC_scores[i].split()[1])
			bind = COM - (TCR + pMHC)
			line.append(bind)
			# BIND peptide alone
			peptide = float(std_pep_scores[i].split()[1])
			TCR_MHC = std_TCR_MHC_score
			bind_pep_alone = COM - (TCR_MHC + peptide)
			line.append(bind_pep_alone)
			# ATLAS
			COM = float(atlas_com_scores[i].split()[1])
			TCR = atlas_TCR_score
			pMHC = float(atlas_pMHC_scores[i].split()[1])
			atlas = COM - (TCR + pMHC)
			line.append(atlas)
			# Write to file
			OUT.write('\t'.join(map(str,line)) + '\n')


		READ_CNT.close()
		STD_COM.close()
		STD_PEP.close()
		STD_PMHC.close()
		ATLAS_COM.close()
		ATLAS_PMHC.close()

	OUT.close()
	
if __name__ == '__main__':
	main()