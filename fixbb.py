#!/usr/bin/env python
import argparse
import subprocess
import os


parser = argparse.ArgumentParser(description='get chi-1 angles for all peptide residues' +
	' after modeling all peptides in TCR-pMHC')
parser.add_argument('-in', help='input score table (3QIB_top50_abundant_peptides_round4.txt)', 
	type=str, required=True, dest="in_file")
parser.add_argument('-start', help='starting residue position for peptide', dest='start', type=int, required=True)
parser.add_argument('-pdb', help='pdb ID (3QIB)', dest='pdb', type=str, required=True)

args = parser.parse_args()
in_file=args.in_file
start=args.start
pdb=args.pdb

# Path to Rosetta3.5 
ros_path = '/home/tb37w/rosetta3.5/rosetta-3.5/'


def fixbb(pdb, resfile, label):
	'''
	Design mutations using Rosetta's fixed backbone application
	'''
	fixbb_cmd = [ros_path + 'rosetta_source/bin/fixbb.linuxgccrelease', '-database',
	ros_path + 'rosetta_database/', '-s', pdb , '-resfile', resfile, '-suffix', '_' + label, 
	'-extrachi_cutoff', '1', '-ex1', '-ex2', '-ex3']
	process = subprocess.Popen(fixbb_cmd)
	process.wait()



def main():

	# Parse peptides
	FH = open(in_file, 'r')
	# Skip header
	FH.readline()
	line_index = 0
	for line in FH:
		line_index += 1
		peptide=line.split()[0]
		# Create resfile
		RESFILE=open('resfile_' + str(line_index), 'w')
		RESFILE.write('NATRO\n')
		RESFILE.write('start\n')
		for i in range(len(peptide)):
			RESFILE.write(str(start + i)  + ' C PIKAA ' + peptide[i] + '\n')
		RESFILE.close()
		
		###############################################################
		# DESIGN
		###############################################################
		
		# Make label for structure
		label = str(line_index)
		resfile = 'resfile_' + label
		# Design complex
		fixbb(pdb + '.trunc.pdb', resfile, label)
		os.remove('score_' + label + '.sc')
		

if __name__ == '__main__':
	main()