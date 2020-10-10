#!/usr/bin/env python
import argparse
import proteindatabank


def main():
	parser = argparse.ArgumentParser(description='Renumber ATOM field from 1..N with N total ATOMS')
	parser.add_argument('-i', help='input PDB file', type=str, required=True)
	args = parser.parse_args()

	myPDB = proteindatabank.PDB(args.i)
	myPDB.renumber_ATOM()


if __name__ == '__main__':
	main()