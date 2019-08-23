#!/usr/bin/env python 

import argparse
import pandas as pd

parser=argparse.ArgumentParser(description='select peptides in top percent of dG_bind and netMHC results')
parser.add_argument('-i', help= 'dG_bind_netMHC_table (ex. dG_bind_netMHC_table_netMHC_sort.txt)', type=str, required=True)
parser.add_argument('-p', help='top percent (ex. .05)', type=float, required=True)
args=parser.parse_args()



def main():

	X = pd.read_csv(args.i, sep='\t')
	total_pep =  X.shape[0]
	top_pep = int(total_pep*args.p)
	dG_sort = X.sort_values('Rosetta_dG_bind')
	dG_thresh = dG_sort.iloc[top_pep-1, 1]
	netMHC_sort = X.sort_values('netMHC_affinity')
	netMHC_thresh = netMHC_sort.iloc[top_pep-1, 2]

	Y = X[(X['Rosetta_dG_bind'] < dG_thresh) & (X['netMHC_affinity'] < netMHC_thresh)]
	Y.to_csv('best_peptides_top_' + str(args.p), sep='\t', index=False)

if __name__ == '__main__':
	main()