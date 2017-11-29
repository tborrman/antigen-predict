#!/usr/bin/env python
import re

peptides = ['LVMPWVVYL', 'NLQGSPVYV',
	'AVQGSTAYL', 'LGYGFVNYI',
	'LLQGWVMYV', 'TMGGYCGYL',
	'DLKGFLSYL', 'SLHGYKKYL',
	'CLFGYDAYV', 'LGYGFVNYV',
	'LLFGFASLV', 'VVFGFLNLV']

A6_lysis = [8.7, 51.9,
	6.8, 58.5,
	58.6, 74.3,
	46.4, 28.3,
	10.7, 57.2,
	17.9, 6.9]

B7_lysis = [3.6, 4.0,
	3.7, 20.6,
	4.2, 6.2,
	4.6, 3.8,
	3.5, 18.4,
	0.7, 5.8]

A6_file = open('/home/tb37w/project/Research/TCR/uniprot/A6/netMHC/all/dG_bind_netMHC_table.txt', 'r')
B7_file = open('/home/tb37w/project/Research/TCR/uniprot/B7/netMHC/all/dG_bind_netMHC_table.txt', 'r')

A6_str = A6_file.read()
B7_str = B7_file.read()
OUT = open('hausmann_predict_table.txt', 'w')
OUT.write('peptide\tA6_Rosetta_dG_bind\tA6_lysis\tB7_Rosetta_dG_bind\tB7_lysis\tnetMHC_affinity\n')

for i, p in enumerate(peptides):
	A6_line = re.search(p + '.*', A6_str).group()
	A6_spltline = A6_line.split()
	A6_Rosetta_dG_bind = A6_spltline[1]
	netMHC_affinity = A6_spltline[2]
	B7_line = re.search(p + '.*', B7_str).group()
	B7_spltline = B7_line.split()
	B7_Rosetta_dG_bind = B7_spltline[1]
	OUT.write('\t'.join([p, A6_Rosetta_dG_bind, str(A6_lysis[i]), 
		B7_Rosetta_dG_bind, str(B7_lysis[i]), netMHC_affinity]) + '\n')
A6_file.close()
B7_file.close()
OUT.close()




	
	
