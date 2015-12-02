#!/bin/bash

#BSUB -q long
#BSUB -W 336:00
#BSUB -R rusage[mem=5000]
#BSUB -J "myarray[1-15716]%1000"
#BSUB -oo out.%I
#BSUB -eo err.%I

MY_FILE=$(printf "peptide_%05d" $LSB_JOBINDEX)
./nr_score.py -pf $MY_FILE

