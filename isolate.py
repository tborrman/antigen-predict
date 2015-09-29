#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(description='Isolate chains from PDB by chain name')
parser.add_argument('-p', help='PDB file', type=str, required=True, dest='pdb')
parser.add_argument('-c', help='chain/chains to isolate separated by commas (ex: D,E)', type=str, required=True, dest='chains')
args = parser.parse_args()

def main():
	chains = args.chains.split(',')
	IN = open(args.pdb, 'r')
	
	for line in IN:
		row = line.split()
		if row[0] =='ATOM' and row[4] in chains:
			print line.rstrip('\n')
	IN.close()
	
if __name__ == '__main__':
	main()