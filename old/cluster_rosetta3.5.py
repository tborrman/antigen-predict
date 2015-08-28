#!/usr/bin/env python

import argparse
import subprocess

def main():
	parser = argparse.ArgumentParser(description= 'Controller for submitting run_rosetta3.5.py to cluster')
	parser.add_argument('-p', help= 'number of processes to run', dest = 'p', type =int, required=True)
	parser.add_argument('-pdb', help='pdb ID (3QIB)', dest='pdb', type=str, required=True)
	parser.add_argument('-start', help='peptide residue start', dest='start', type=int, required=True)
	parser.add_argument('-pMHC', help='calculate dG for pMHC complex', dest='pMHC', action='store_true')
	parser.add_argument('-w', help='weights file for scoring function (standard.wts)', dest='w')

	args=parser.parse_args()
	proc=args.p
	pdb_ID=args.pdb
	start=args.start
	pMHC=args.pMHC
	w=args.w

	if pMHC:
		# Submit rosetta3.5 jobs to cluster
		if w:
			for part in range(1, proc + 1):
				process = subprocess.Popen(['bsub', '-q', 'short', '-n', '2', '-W', '240', '-e', 'job_' + str(part) + '_err', '-m', 'blades', '-J', 'rosetta_' + str(part),'-R', 'rusage[mem=5000]', 'run_rosetta3.5.py', 
					'-in', 'redundant_out', '-part', str(part), '-pdb', pdb_ID, '-start', str(start), '-pMHC', '-w', w])  
		else:
			for part in range(1, proc + 1):
				process = subprocess.Popen(['bsub', '-q', 'short', '-n', '2', '-W', '240', '-e', 'job_' + str(part) + '_err', '-m', 'blades', '-J', 'rosetta_' + str(part),'-R', 'rusage[mem=5000]', 'run_rosetta3.5.py', 
					'-in', 'redundant_out', '-part', str(part), '-pdb', pdb_ID, '-start', str(start), '-pMHC'])  
	else:
		
		# Submit rosetta3.5 jobs to cluster
		if w:
			for part in range(1, proc + 1):
				process = subprocess.Popen(['bsub', '-q', 'short', '-n', '2', '-W', '240', '-e', 'job_' + str(part) + '_err', '-m', 'blades', '-J', 'rosetta_' + str(part),'-R', 'rusage[mem=5000]', 'run_rosetta3.5.py', 
					'-in', 'redundant_out', '-part', str(part), '-pdb', pdb_ID, '-start', str(start), '-w', w])
		else:
			for part in range(1, proc + 1):
				process = subprocess.Popen(['bsub', '-q', 'short', '-n', '2', '-W', '240', '-e', 'job_' + str(part) + '_err', '-m', 'blades', '-J', 'rosetta_' + str(part),'-R', 'rusage[mem=5000]', 'run_rosetta3.5.py', 
					'-in', 'redundant_out', '-part', str(part), '-pdb', pdb_ID, '-start', str(start)])

if __name__ == '__main__':
		main()