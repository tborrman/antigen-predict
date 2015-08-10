#!/usr/bin/env python

import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='Perform scoring methods via Rosetta 3.5 on Birnbaum selections')
parser.add_argument('-in', help= 'input_files prefix (redundant_out)', dest = 'in_file', type=str, required=True)
parser.add_argument('-part', help= 'part number for parallel processing', dest = 'p', type =int, required=True)
parser.add_argument('-pdb', help='pdb ID (3QIB)', dest='pdb', type=str, required=True)
parser.add_argument('-start', help='starting residue position for peptide', dest='start', type=int, required=True)
parser.add_argument('-w', help='weights file for ATLAS trained scoring', dest='w', required=True)

args=parser.parse_args()
in_file=args.in_file
part=args.p
pdb=args.pdb
start=args.start
w=args.w

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

def score(pdb, label, out_dir):
	'''
	Score PDB complex using Rosetta3 standard scoring function
	'''
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	score_cmd = [ros_path + 'rosetta_source/bin/score.linuxgccrelease', '-database',
	ros_path + 'rosetta_database/', '-s', pdb , '-out:file:scorefile', 
	out_dir + '/' + label + '.sc', '-extrachi_cutoff', '12', '-ex1', '-ex2', '-ex3']
	process = subprocess.Popen(score_cmd)
	process.wait()

def weights_score(pdb, label, weights, out_dir):
	'''
	Score PDB complex using specified weights file
	'''
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	score_cmd = [ros_path + 'rosetta_source/bin/score.linuxgccrelease', '-database',
	ros_path + 'rosetta_database/', '-s', pdb , '-out:file:scorefile', 
	out_dir + '/' + label + '.sc', '-extrachi_cutoff', '12', '-ex1', '-ex2', '-ex3',
	'-score:weights', weights]
	process = subprocess.Popen(score_cmd)
	process.wait()
	
def isolate(pdb, label):
	'''
	Isolate peptide and pMHC
	'''
	
	IN = open(pdb, 'r')
	OUT = open('pMHC_' + label + '.pdb', 'w')
	for line in IN:
		row = line.split()
		if row[0] == 'ATOM':
			if row[4] == 'A' or row[4] == 'B' or row[4] == 'C':
				OUT.write(line)
	IN.close()
	OUT.close()

	IN = open(pdb, 'r')
	OUT = open('peptide_' + label + '.pdb', 'w')
	for line in IN:
		row = line.split()
		if row[0] == 'ATOM':
			if row[4] == 'C':
				OUT.write(line)
	IN.close()
	OUT.close()

def main():

	# Open Birnbaum selections
	line_index = 0
	FH = open(in_file + '_' + str(part) +'.txt', 'r')
	for line in FH:
		line_index += 1
		peptide=line.split()[0]
		# Create resfile
		RESFILE=open('resfile_' + str(part) + '_' + str(line_index), 'w')
		RESFILE.write('NATRO\n')
		RESFILE.write('start\n')
		for i in range(len(peptide)):
			RESFILE.write(str(start + i)  + ' C PIKAA ' + peptide[i] + '\n')
		RESFILE.close()
		
		###############################################################
		# DESIGN
		###############################################################
		
		# Make label for structure
		label = '_'.join([str(part), str(line_index)])
		resfile = 'resfile_' + label
		# Design complex
		fixbb(pdb + '.trunc.pdb', resfile, label)
		os.remove('score_' + label + '.sc')
		os.remove('resfile_' + label)
		
		##############################################################
		# SCORING
		##############################################################

		# STANDARD (standard Rosetta scoring function)
		out_dir = 'standard_score'
		# Score the designed complex
		designed_complex = pdb + '.trunc_' + label + '_0001.pdb'
		score(designed_complex, 'com_' + label, out_dir)
		
		# Score pMHC and peptide individually
		# Isolate peptide and pMHC
		isolate(designed_complex, label) 
		# Score pMHC
		score('pMHC_' + label + '.pdb', 'pMHC_' + label, out_dir)
		# Score peptide
		score('peptide_' + label + '.pdb', 'peptide_' + label, out_dir)

		# ATLAS (optimized ATLAS scoring function)
		out_dir = 'atlas_score'
		weights_score(designed_complex, 'com_' + label, w, out_dir)
		weights_score('pMHC_' + label + '.pdb', 'pMHC_' + label, w, out_dir)

		# Remove temp files
		os.remove('pMHC_' + label + '.pdb')
		os.remove('peptide_' + label + '.pdb')
		os.remove(designed_complex)

	#########################################################
	# AGGREGATE RESULTS
	#########################################################
	# Standard directory
	os.chdir('standard_score')
	pMHC_score = open('pMHC_score_' + str(part) + '.sc', 'w')
	peptide_score = open('peptide_score_' + str(part) + '.sc', 'w')
	com_score = open('com_score_' + str(part) + '.sc', 'w')
	for i in range(1,line_index + 1):
		pMHC=open('pMHC_' + str(part) + '_' + str(i) + '.sc', 'r')
		pMHC_lines=pMHC.readlines()
		pMHC_score.write(pMHC_lines[1])
		pMHC.close()
		COM=open('com_' + str(part) + '_' + str(i) + '.sc', 'r')
		com_lines=COM.readlines()
		com_score.write(com_lines[1])
		COM.close()
		PEPTIDE=open('peptide_' + str(part) + '_' + str(i) + '.sc', 'r')
		peptide_lines=PEPTIDE.readlines()
		peptide_score.write(peptide_lines[1])
		PEPTIDE.close()
	pMHC_score.close()
	peptide_score.close()
	com_score.close()
	# Remove used temp files
	for i in range(1, line_index + 1):
		os.remove('pMHC_' + str(part) + '_' + str(i) + '.sc')
		os.remove('com_' + str(part) + '_' + str(i) + '.sc')
		os.remove('peptide_' + str(part) + '_' + str(i) + '.sc')
	# ATLAS directory
	os.chdir('../atlas_score')
	pMHC_score = open('pMHC_score_' + str(part) + '.sc', 'w')
	com_score = open('com_score_' + str(part) + '.sc', 'w')
	for i in range(1,line_index + 1):
		pMHC=open('pMHC_' + str(part) + '_' + str(i) + '.sc', 'r')
		pMHC_lines=pMHC.readlines()
		pMHC_score.write(pMHC_lines[1])
		pMHC.close()
		COM=open('com_' + str(part) + '_' + str(i) + '.sc', 'r')
		com_lines=COM.readlines()
		com_score.write(com_lines[1])
		COM.close()
	pMHC_score.close()
	com_score.close()
	# Remove used temp files
	for i in range(1, line_index + 1):
		os.remove('pMHC_' + str(part) + '_' + str(i) + '.sc')
		os.remove('com_' + str(part) + '_' + str(i) + '.sc')

	FH.close()

if __name__ == '__main__':
	main()
