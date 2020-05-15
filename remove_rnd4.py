#!/usr/bin/env python
import argparse

parser=argparse.ArgumentParser(description='remove round 4 peptides from preselection library')
parser.add_argument('-p', help='preselection score table (score_table_hamming_ENPVVHFFKNIVTP_preselection.txt)', 
	type=str, required=True )
parser.add_argument('-r', help='round 4 score table (score_table_hamming_ENPVVHFFKNIVTP_round4.txt)',
	type=str, required=True)
args=parser.parse_args()

def in_round4(pep, f):
	'''
	Search for pep in round 4 library
	Args:
		pep: string query peptide sequence
		f: round 4 score table file
	Return:
		found: boolean True if found
	'''
	found = False
	IN = open(f, 'r')
	IN.readline()
	for line in IN:
		rpep = line.split()[0]
		if pep == rpep:
			found = True
			IN.close()
			return found
	IN.close()
	return found

def main():

	PRE = open(args.p, 'r')
	# Remove header
	PRE.readline()
	# Parse preselection library
	for i, line in enumerate(PRE):
		# if i % 1000 == 0:
		# 	print "on row " + str(i)
		pre_peptide = line.split()[0]
		# Search in round 4 library
		if not in_round4(pre_peptide, args.r):
			print line.strip()
	PRE.close()
	
if __name__ == '__main__':
	main()