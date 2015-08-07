#!/usr/bin/env python

import argparse

def main():
	parser=argparse.ArgumentParser(description='Split redundant read count file into parts')
	parser.add_argument('-in', help='redundant peptide count file', type=str, dest= 'infile', required=True)
	parser.add_argument('-n', help='number of split output files to create (NOTE: remainder files created if not divisible)', type=int, dest='n', required=True)

	args=parser.parse_args()

	# Count number of lines  
	IN = open(args.infile, 'r')
	pep_list = IN.readlines()  
	tot_pep = len(pep_list)
	IN.close()

	if tot_pep < args.n:
		print 'ERROR: number of files to create is greater than number of peptides'
		quit()

	part_pep = tot_pep/args.n

	# Write files to directory
	start = 0
	stop = part_pep
	for out_num in range(1, args.n + 1):
		OUT = open('redundant_out_' + str(out_num) + '.txt', 'w+')
		for part in range(start, stop):
			OUT.write(pep_list[part])
		start = start + part_pep
		stop = stop + part_pep
		OUT.close()

	if tot_pep%args.n != 0:
		remainder = tot_pep%args.n
		rem_peps = pep_list[-remainder:]
		counter = 0
		for x in range(0, len(rem_peps), 30):
			counter += 1
			chunk = rem_peps[x:x+30]
			OUT = open('redundant_out_' + str(args.n + counter) + '.txt', 'w+')
			for line in chunk:
				OUT.write(line)
			OUT.close()



if __name__ == '__main__':
	main()
