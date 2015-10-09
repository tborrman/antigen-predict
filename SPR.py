#!/usr/bin/env python
import argparse
import os
import subprocess

parser=argparse.ArgumentParser(description='Design and score peptides from Birnbaum SPR data')
parser.add_argument('-pep', help='peptide sequence file (peptides.txt)', type=str, dest='pep', required=True)
parser.add_argument('-pdb', help='pdb ID (3QIB)', type=str, dest='pdb', required=True)
parser.add_argument('-start', help='starting residue position for peptide', type=int, dest='start', required=True)
args = parser.parse_args()

# Path to Rosetta3.5 
ros_path = '/home/tb37w/rosetta3.5/rosetta-3.5/'


def main():
	
	# Open peptides file 
	PEP = open(args.pep, 'r')
	OUT = open('peptides_scores.txt', 'w')
	for line in PEP:
		peptide = line.split()[0]
		peptide_id = line.split()[1]
		# Create resfile
		RESFILE=open('resfile_' + peptide_id, 'w')
		RESFILE.write('NATRO\n')
		RESFILE.write('start\n')
		for i in range(len(peptide)):
			RESFILE.write(str(args.start + i)  + ' C PIKAA ' + peptide[i] + '\n')
		RESFILE.close()

		###############################################################
		# DESIGN
		###############################################################
		
		# Make label for structure
		label = peptide_id
		resfile = 'resfile_' + label
		# Design complex
		fixbb(args.pdb + '.trunc.pdb', resfile, label)
		os.remove('score_' + label + '.sc')
		
		##############################################################
		# SCORING
		##############################################################

		# Score the designed complex
		designed_complex = args.pdb + '.trunc_' + label + '_0001.pdb'
		score(designed_complex, 'com_' + label)

		# Score peptide individually
		# Isolate peptide
		isolate(designed_complex, label)
		# Score peptide
		score('peptide_' + label + '.pdb', 'peptide_' + label) 

		#########################################################
		# AGGREGATE RESULTS
		#########################################################
		TCR_MHC_file = open('TCR_MHC_score.sc', 'r')
		COM_file= open('com_' + label + '.sc', 'r')
		PEP_file = open('peptide_' + label + '.sc', 'r')
		TCR_MHC_file.readline()
		COM_file.readline()
		PEP_file.readline()

		TCR_MHC_score = float(TCR_MHC_file.readline().split()[1])
		COM_score = float(COM_file.readline().split()[1])
		PEP_score = float(PEP_file.readline().split()[1])

		bind_score = COM_score - (TCR_MHC_score + PEP_score)

		OUT.write(line.rstrip() + '\t' + str(bind_score) + '\n')
		TCR_MHC_file.close()
		COM_file.close()
		PEP_file.close()
	PEP.close()
	OUT.close()




def fixbb(pdb, resfile, label):
	'''
	Design mutations using Rosetta's fixed backbone application
	'''
	fixbb_cmd = [ros_path + 'rosetta_source/bin/fixbb.linuxgccrelease', '-database',
	ros_path + 'rosetta_database/', '-s', pdb , '-resfile', resfile, '-suffix', '_' + label, 
	'-extrachi_cutoff', '1', '-ex1', '-ex2', '-ex3']
	process = subprocess.Popen(fixbb_cmd)
	process.wait()

def score(pdb, label):
	'''
	Score PDB complex using Rosetta3 standard scoring function
	'''
	score_cmd = [ros_path + 'rosetta_source/bin/score.linuxgccrelease', '-database',
	ros_path + 'rosetta_database/', '-s', pdb , '-out:file:scorefile', 
	label + '.sc']
	process = subprocess.Popen(score_cmd)
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


if __name__ == '__main__':
	main()