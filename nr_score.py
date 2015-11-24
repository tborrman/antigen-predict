#!/usr/bin/env python

import argparse
import sys
import os
import subprocess

# Flush STOUT continuously
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

parser = argparse.ArgumentParser(description='Calculate scores for nr database peptides')
parser.add_argument('-pf', help='peptide file (ex. peptide_00000)', type=str, dest='pf', required=True)
args = parser.parse_args()

# Path to Rosetta3.5 
ros_path = '/home/tb37w/rosetta3.5/rosetta-3.5/'

# Starting peptide position
start = 85
# PDB
pdb = '1YMM'
# TCR-MHC score for 1YMM
tcr_mhc_score = 7517.872
# For suppressing Rosetta output
FNULL = open(os.devnull, 'w')

def fixbb(pdb, resfile, label):
	'''
	Design mutations using Rosetta's fixed backbone application
	'''
	fixbb_cmd = [ros_path + 'rosetta_source/bin/fixbb.linuxgccrelease', '-database',
	ros_path + 'rosetta_database/', '-s', pdb , '-resfile', resfile, '-suffix', '_' + label, 
	'-extrachi_cutoff', '1', '-ex1', '-ex2', '-ex3']
	process = subprocess.Popen(fixbb_cmd, stdout=FNULL)
	process.wait()

def score(pdb, label):
	'''
	Score PDB complex using Rosetta3 standard scoring function
	'''
	score_cmd = [ros_path + 'rosetta_source/bin/score.linuxgccrelease', '-database',
	ros_path + 'rosetta_database/', '-s', pdb , '-out:file:scorefile', label + '.sc']
	process = subprocess.Popen(score_cmd, stdout=FNULL)
	process.wait()

def isolate(pdb, label):
	'''
	Isolate peptide
	'''
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
	# Open score output file
	index = args.pf[-5:]
	OUT=open('scores_' + index, 'w')
	# Parse input
	IN = open(args.pf, 'r')
	line_index = 0
	for line in IN:
		line_index += 1
		if line_index%10 == 0:
			print 'Peptide #: '+str(line_index) 
		peptide = line.strip()
		# Create resfile
		RESFILE = open('resfile_' + index + '_' + str(line_index), 'w')
		RESFILE.write('NATRO\n')
		RESFILE.write('start\n')
		for i in range(len(peptide)):
			RESFILE.write(str(start + i)  + ' C PIKAA ' + peptide[i] + '\n')
		RESFILE.close()

		###############################################################
		# DESIGN
		###############################################################
		
		# Make label for structure
		label = '_'.join([str(index), str(line_index)])
		resfile = 'resfile_' + label
		# Design complex
		fixbb(pdb + '.trunc.pdb', resfile, label)
		os.remove('score_' + label + '.sc')
		os.remove('resfile_' + label)

		##############################################################
		# SCORING
		##############################################################

		# Score the designed complex
		designed_complex = pdb + '.trunc_' + label + '_0001.pdb'
		score(designed_complex, 'com_' + label)

		# Score peptide individually
		# Isolate peptide 
		isolate(designed_complex, label) 
		# Score peptide
		score('peptide_' + label + '.pdb', 'peptide_' + label)

		# Remove temp files
		os.remove('peptide_' + label + '.pdb')
		os.remove(designed_complex)

		#########################################################
		# AGGREGATE RESULTS
		#########################################################

		COM=open('com_' + label + '.sc', 'r')
		# Skip header
		COM.readline()
		com_score = float(COM.readline().split()[1])
		COM.close()
		os.remove('com_' + label + '.sc')
		
		PEP=open('peptide_' + label + '.sc', 'r')
		PEP.readline()
		pep_score = float(PEP.readline().split()[1])
		PEP.close()
		os.remove('peptide_' + label + '.sc')

		bind_score = com_score - (tcr_mhc_score + pep_score)

		OUT.write(peptide +'\t' + str(bind_score) + '\n')
	
	IN.close()
	OUT.close()

if __name__ == '__main__':
	main()