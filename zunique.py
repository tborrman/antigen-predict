#!/usr/bin/env python
import argparse
from Bio import SeqIO
from Bio.Alphabet import IUPAC
import os
import subprocess

parser=argparse.ArgumentParser(description="Zhiping's indexing solution for counting unique peptides")
parser.parse_args()


aa_letters = IUPAC.ExtendedIUPACProtein().letters
for aa1 in aa_letters:
	for aa2 in aa_letters:
		# Make indexed directories with first and second amino acids
		if not os.path.exists(aa1+aa2):
			os.makedirs(aa1+aa2)
		os.chdir(aa1+aa2)
		# Write script for cluster
		SRC = open(aa1+aa2+"_hash_count.py", 'w')
		SRC.write(
			"from Bio import SeqIO\n"+
			"import sys\n"+
			"import os\n"+
			"# Flush STOUT continuously\n"+
			"sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)\n"+
			"OUT = open('"+aa1+aa2+"_peptides.txt', 'w')\n"+
			"counter = 0\n"+
			"hash_count = {}\n"
			"# Loop through records\n"+
			"for seq_record in SeqIO.parse('../../nr.fasta', 'fasta'):\n"+
			"	if counter%1000000 == 0:\n"+
			"		print 'On seq ' + str(counter) + '\\n'\n"+
			"	for i in range(len(seq_record.seq) - 14 + 1):\n"+
			"		if str(seq_record.seq[i:i+2]) == '" + aa1 + aa2 +"':\n"+
            "			hash_count[str(seq_record.seq[i:i+14])] = None\n"+
        	"	counter +=1\n"+
			"print 'Number of sequences checked = ' + str(counter) + '\\n'\n"+	
			"# Write unique peptides to file\n"+
			"for key in hash_count.keys():\n"+
			"	OUT.write(key + '\\n')\n"+
			"OUT.close()\n"
			)
		SRC.close()
		# Run script on cluster
		process = subprocess.Popen(['bsub', '-q', 'long', '-W', '144:00', '-R', 'rusage[mem=100000]', '-o', aa1+aa2+'_peptides.out',
			'-e',  aa1+aa2+'_peptides.err', '-J', aa1+aa2, 'python', aa1+aa2+'_hash_count.py'])
		process.wait()	 
		os.chdir("../")





		







