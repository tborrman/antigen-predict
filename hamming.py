#!/usr/bin/env python
import argparse
parser = argparse.ArgumentParser(description='Record Hamming distance between input WT peptide and each peptide in score_table')
parser.add_argument('-wt', help='wildtype peptide', type=str, dest='wt', required=True)
parser.add_argument('-f', help='score table (ex. unique_pre_rnd4_score_table.txt)', type=str, dest='f', required=True)
args = parser.parse_args()


def main():
	IN = open(args.f, 'r')
	header = IN.readline()
	print header.rstrip() + '\t' + 'hamming'
	wt = args.wt
	for line in IN:
		query = line.split()[0]
		check_length(wt, query)
		hd = hamming(wt, query)
		print line.rstrip() + '\t' + hd
	IN.close()

def check_length(str1, str2):
	if len(str1) != len(str2):
		print 'ERROR: strings of unequal length'
		quit()


def hamming(str1, str2):
	hd = 0
	for i in range(len(str1)):
		if str1[i] != str2[i]:
			hd +=1
	return str(hd)


if __name__ == '__main__':
 	main() 
 	