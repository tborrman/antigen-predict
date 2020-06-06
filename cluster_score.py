#!/usr/bin/env python

import argparse
import subprocess

def main():
	parser = argparse.ArgumentParser(description= 'Controller for submitting score.py to cluster')
	parser.add_argument('-p', help= 'number of processes to run', dest = 'p', type =int, required=True)
	parser.add_argument('-pdb', help='pdb ID (3QIB)', dest='pdb', type=str, required=True)
	parser.add_argument('-start', help='peptide residue start', dest='start', type=int, required=True)
	parser.add_argument('-w', help='weights file for ATLAS trained scoring (ATLAS_weights.wts)', dest='w')

	args=parser.parse_args()
	proc=args.p
	pdb_ID=args.pdb
	start=args.start
	w=args.w

	# Submit rosetta3.5 jobs to cluster
	
	for part in range(1, proc + 1):
		process = subprocess.Popen(['bsub', '-q', 'long', '-W', '36:00', '-n', '1', '-e', 'job_' + str(part) + '_err', '-J', 'rosetta_' + str(part),'-R', 'rusage[mem=5000]', 'score.py', 
			'-in', 'redundant_out', '-part', str(part), '-pdb', pdb_ID, '-start', str(start), '-w', w])  

if __name__ == '__main__':
		main()