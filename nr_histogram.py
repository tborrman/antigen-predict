#!/usr/bin/env python

import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Manually create histogram heights to plot in R from peptide scores")
parser.add_argument('-p', help='peptide score file', type=str, dest="p", required=True)
args=parser.parse_args()



def main():
	histogram = dict.fromkeys(['{:.2f}'.format(x) for x in np.arange(-4000, 4000, 0.01)], 0)
	
	IN = open(args.p, 'r')
	for i, line in enumerate(IN):
		 histogram['{:.2f}'.format(float(line.split()[1]))] += 1
		 if i%1000 == 0:
		 	print 'On row: ' + str(i)

		 	
					
	IN.close()
	OUT = open('peptide_histogram.txt', 'w')
	sorted_keys = ['{:.2f}'.format(x) for x in sorted(map(float,histogram.keys()))]
	for key in sorted_keys:
		OUT.write(key + '\t' + str(histogram[key]) + '\n')
	OUT.close()







if __name__ == '__main__':
	main()