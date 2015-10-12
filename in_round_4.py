#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser(description='Check if top 50 scoring are in round 4 (should be run in [PDB]/ folder)')
parser.add_argument('-s', help='scoring function (ex. standard, bind, bind_pep_alone, atlas)', type=str, dest='s', required=True)
args = parser.parse_args()

def main():

	# Store preselection and round 4 reads in dictionary
	PRE = open('preselection/all/score_table.txt', 'r')
	RND4 = open('round_4/all/score_table.txt', 'r')

	scores = {}
	rnd4_scores = []
	header = PRE.readline().split()
	for line in PRE:
		scores[line.split()[0]] = float(line.split()[header.index(args.s)])
	RND4.readline()
	for line in RND4:
		scores[line.split()[0]] = float(line.split()[header.index(args.s)])
		rnd4_scores.append(line.split()[0])
	PRE.close()
	RND4.close()

	# Sort scoring
	sorted_scoring = sorted(scores.keys(), key = lambda k: scores[k])
	top50 = sorted_scoring[:50]
	# Check if top 50 scoring are in round 4
	counter = 0
	for peptide in top50:
		if peptide in rnd4_scores:
			counter+=1
	print str(counter) + ' of top 50 scoring peptides are in round 4 (' + str((counter/50.0) * 100) + '%).'

if __name__ == '__main__':
	main()