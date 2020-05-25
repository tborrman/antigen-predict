#!/usr/bin/env python
import argparse
from Bio.SubsMat import MatrixInfo 
import sys


parser = argparse.ArgumentParser(description='Record BLOSUM62 similarity score between input WT peptide and each peptide in score_table')
parser.add_argument('-wt', help='wildtype peptide', type=str, dest='wt', required=True)
parser.add_argument('-f', help='score table (ex. unique_pre_rnd4_score_table.txt)', type=str, dest='f', required=True)
args = parser.parse_args()


def check_length(str1, str2):
	if len(str1) != len(str2):
		print 'ERROR: strings of unequal length'
		sys.exit()


def blosum_score(str1, str2):
	'''
	Args: 
		str1: string wildtype peptide
		str2: string query peptide
	Return:
		bs: float BLOSUM62 similarity score 
			between str1, str2
	'''

	blosum_dict = MatrixInfo.blosum62
	bs = 0.0

	for i in range(len(str1)):
		aa_pair = (str1[i], str2[i])
		rev_aa_pair = aa_pair[::-1] 
		if aa_pair in blosum_dict.keys():
			bs += blosum_dict[aa_pair]
		elif rev_aa_pair in blosum_dict.keys():  
			bs += blosum_dict[rev_aa_pair]
		else: 
			print 'ERROR: unknown amino acid tuple'
			sys.exit()
	return bs

def main():

	IN = open(args.f, 'r')
	header = IN.readline()
	print header.rstrip() + '\t' + 'blosum62'
	wt = args.wt
	for line in IN:
		query = line.split()[0]
		check_length(wt, query)
		bs = blosum_score(wt, query)
		print line.rstrip() + '\t' + str(bs)
	IN.close()



if __name__ == '__main__':
 	main() 
 	