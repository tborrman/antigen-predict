#!/usr/bin/env python
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Calculate bind_pep_alone score from score files')
parser.add_argument('-c', help='score file for complex (ex. com_score.sc)', type=str, required=True, dest='com')
parser.add_argument('-t', help='score file for TCR/MHC (ex. TCR_MHC_score.sc)', type=str, required=True, dest='tcr_mhc')
parser.add_argument('-p', help='score file for peptide (ex. pep_score.sc)', type=str, required=True, dest='pep')
args = parser.parse_args()

def main():
	COM = open(args.com, 'r')
	TCR_MHC = open(args.tcr_mhc, 'r')
	PEP = open(args.pep, 'r')

	header = '\t'.join(COM.readline().split()[1:-1])
	print header
	TCR_MHC.readline()
	PEP.readline()

	com_scores = np.array(map(float, COM.readline().split()[1:-1]))
	tcr_mhc_scores = np.array(map(float, TCR_MHC.readline().split()[1:-1]))
	pep_scores = np.array(map(float, PEP.readline().split()[1:-1]))

	bind_pep_alone = com_scores - (tcr_mhc_scores + pep_scores)
	bind_out = map('{:.3f}'.format, bind_pep_alone)

	print '\t'.join(bind_out)

	COM.close()
	TCR_MHC.close()
	PEP.close()




if __name__ == '__main__':
	main()