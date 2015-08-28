#!/usr/bin/env python

import argparse
import os
import subprocess

def main():
	parser = argparse.ArgumentParser(description='Run Rosetta 3.5 on Birnbaum selections')
	parser.add_argument('-in', help= 'input_files prefix (redundant_out)', dest = 'in_file', type=str, required=True)
	parser.add_argument('-part', help= 'part number for parallel processing', dest = 'p', type =int, required=True)
	parser.add_argument('-pdb', help='pdb ID (3QIB)', dest='pdb', type=str, required=True)
	parser.add_argument('-start', help='starting residue position for peptide', dest='start', type=int, required=True)
	parser.add_argument('-pMHC', help='calculate dG for pMHC complex', dest='pMHC', action='store_true')
	parser.add_argument('-w', help='weights file for scoring function (standard.wts)', dest='w')

	args=parser.parse_args()
	in_file=args.in_file
	part=args.p
	pdb=args.pdb
	start=args.start
	pMHC=args.pMHC
	w=args.w

	if pMHC:
		pep_pdb_file = pdb + '_pMHC.trunc.pdb'
		
	else:
		pep_pdb_file = pdb + '_peptide.trunc.pdb'

	# Process array
	processes = []
	# Output files
	# Peptide scores
	PEP_SCORE = open('pep_score_' + str(part) + '.sc', 'w')
	# Complex scores
	COM_SCORE = open('com_score_' + str(part) + '.sc', 'w')

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
		if w:
			# Score peptide/pMHC
			process = subprocess.Popen(['/home/tb37w/rosetta3.5/rosetta-3.5/rosetta_source/bin/fixbb.linuxgccrelease', '-database', 
			'/home/tb37w/rosetta3.5/rosetta-3.5/rosetta_database/', '-s', pep_pdb_file, '-resfile', 'resfile_' + str(part) + '_' + str(line_index),
			'-out:file:scorefile','pep_score_' + str(part) + '_' + str(line_index) + '.sc', '-out:file:score_only', '-extrachi_cutoff', '12', 
			'-ex1', '-ex2', '-ex3', '-constant_seed', '-jran', '14', '-score:weights', w])
			processes.append(process)
			# Score complex
			process = subprocess.Popen(['/home/tb37w/rosetta3.5/rosetta-3.5/rosetta_source/bin/fixbb.linuxgccrelease', '-database', 
			'/home/tb37w/rosetta3.5/rosetta-3.5/rosetta_database/', '-s', pdb + '.trunc.pdb', '-resfile', 'resfile_' + str(part) + '_' + str(line_index),
			'-out:file:scorefile','com_score_' + str(part) + '_' + str(line_index) + '.sc', '-out:file:score_only', '-extrachi_cutoff', '12', 
			'-ex1', '-ex2', '-ex3', '-constant_seed', '-jran', '14', '-score:weights', w])
			processes.append(process)
		else:
			# Score peptide/pMHC
			process = subprocess.Popen(['/home/tb37w/rosetta3.5/rosetta-3.5/rosetta_source/bin/fixbb.linuxgccrelease', '-database', 
			'/home/tb37w/rosetta3.5/rosetta-3.5/rosetta_database/', '-s', pep_pdb_file, '-resfile', 'resfile_' + str(part) + '_' + str(line_index),
			'-out:file:scorefile','pep_score_' + str(part) + '_' + str(line_index) + '.sc', '-out:file:score_only', '-extrachi_cutoff', '12', 
			'-ex1', '-ex2', '-ex3', '-constant_seed', '-jran', '14'])
			processes.append(process)
			# Score complex
			process = subprocess.Popen(['/home/tb37w/rosetta3.5/rosetta-3.5/rosetta_source/bin/fixbb.linuxgccrelease', '-database', 
			'/home/tb37w/rosetta3.5/rosetta-3.5/rosetta_database/', '-s', pdb + '.trunc.pdb', '-resfile', 'resfile_' + str(part) + '_' + str(line_index),
			'-out:file:scorefile','com_score_' + str(part) + '_' + str(line_index) + '.sc', '-out:file:score_only', '-extrachi_cutoff', '12', 
			'-ex1', '-ex2', '-ex3', '-constant_seed', '-jran', '14'])
			processes.append(process)

		for proc in processes:
			proc.wait()
		processes = []
		
	# for proc in processes:
	# 	proc.wait()

	# Aggregate scores
	for i in range(1,line_index + 1):
		PEP=open('pep_score_' + str(part) + '_' + str(i) + '.sc', 'r')
		pep_lines=PEP.readlines()
		PEP_SCORE.write(pep_lines[2])
		PEP.close()
		COM=open('com_score_' + str(part) + '_' + str(i) + '.sc', 'r')
		com_lines=COM.readlines()
		COM_SCORE.write(com_lines[2])
		COM.close()
	PEP_SCORE.close()
	COM_SCORE.close()
	FH.close()		

	# Remove used temp files
	for i in range(1, line_index + 1):
		os.remove('resfile_' + str(part) + '_' + str(i))
		os.remove('pep_score_' + str(part) + '_' + str(i) + '.sc')
		os.remove('com_score_' + str(part) + '_' + str(i) + '.sc')




if __name__ == '__main__':
	main()
