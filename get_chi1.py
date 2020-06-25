#!/usr/bin/env python
import argparse
from pyrosetta import *
from rosetta import *
from pyrosetta import toolbox
import os
init()

parser=argparse.ArgumentParser(description='get chi-1 angles for residues of peptide in modeled TCR-pMHC complex')
parser.add_argument('-in', help='input score table (3QIB_top50_abundant_peptides_round4.txt)',
        type=str, required=True, dest="in_file")
parser.add_argument('-start', help='starting residue position for peptide', dest='start', type=int, required=True)
parser.add_argument('-pdb', help='pdb ID (3QIB)', dest='pdb', type=str, required=True)

args = parser.parse_args()
in_file=args.in_file
start=args.start
pdb=args.pdb

def get_peptide_length(file):
	'''
	Get peptide length
	'''
	with open(file, 'r') as f:
		f.readline()
		l = len(f.readline().split()[0]) 
	return l

def calculate_chi(my_pose, res_num):
	'''
	Calculate chi angle and rotamer states for residue
	'''
	new_vec = rosetta.utility.vector1_unsigned_long()
	testy = rosetta.core.pack.dunbrack.rotamer_from_chi(my_pose.residue(res_num), new_vec)
	print my_pose.residue(res_num).name()
	chi_angles =  my_pose.residue(res_num).chi()
	print chi_angles
	rotamer_states = map(int, list(new_vec))
	print map(int, list(new_vec))
	return chi_angles, rotamer_states

def get_pose(pdb):
	'''
	Return Pyrosetta pose
	'''

	toolbox.cleanATOM(pdb)
	pose = pose_from_pdb(pdb[:-4] + '.clean.pdb')
	return pose


def main():

	pep_length = get_peptide_length(in_file)
	# Parse peptides
	FH = open(in_file, 'r')
	OUTchi = open(in_file[:-4] + '_chi1.txt', 'w')
	OUTrotamer = open(in_file[:-4] + '_rotamer1.txt', 'w')
	# New header
	header = FH.readline()
	OUTchi.write(header.strip() + '\t' + 
		'\t'.join(map(str, range(start, start + pep_length))) + '\n')
	OUTrotamer.write(header.strip() + '\t' + 
		'\t'.join(map(str, range(start, start + pep_length))) + '\n')
	line_index = 0
	for line in FH:
		line_index += 1
		OUTchi.write('\t'.join(line.split()[:6]) + '\t')
		OUTrotamer.write('\t'.join(line.split()[:6]) +'\t')
		p = get_pose(pdb + '.trunc_' + str(line_index) + '_0001.pdb')
		for i in range(pep_length):
			res_pos = start + i
			absolute_res_num = p.pdb_info().pdb2pose('C', res_pos)
			chis, rotamers = calculate_chi(p, absolute_res_num)
			if len(chis) == 0:
				chi1 = 'NA'	
				rotamer1 = 'NA'
			else:
				chi1 = '{:.2f}'.format(float(list(chis)[0]))
				rotamer1 = str(rotamers[0])
			if i < (pep_length-1):
				OUTchi.write(chi1 + '\t')
				OUTrotamer.write(rotamer1 + '\t')
			else:
				OUTchi.write(chi1 + '\n')
				OUTrotamer.write(rotamer1 + '\n')
	FH.close()
	OUTchi.close()	
	OUTrotamer.close()
			
			



if __name__ == '__main__':
	main()


