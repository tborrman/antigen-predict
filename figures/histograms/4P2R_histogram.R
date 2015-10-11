#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description="Make histograms for Rosetta 3.5 scoring results")
parser$add_argument("-in", help="score table output by aggregate_scores.py (score_table.txt)", type="character", dest="infile", required=TRUE)
parser$add_argument("-r", help="selection round (ex. 3)", type="integer", dest="round", required=TRUE)

args <- parser$parse_args()
wkdir <- getwd()

score_table <-read.table(paste(wkdir, "/", args$infile, sep=""), header=TRUE, sep="\t")
header <- colnames(score_table)

png(paste(wkdir, "/", header[3], "_histogram.png", sep=""))
par(mar=c(5,5,4,2) + 0.1)
hist(score_table$standard, breaks=seq(-1000, 7000, by= 10), col='forestgreen', lty="blank", freq = FALSE, ylab = 'Density', xlab= "Rosetta 3.5 standard score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(-100,3000), ylim=c(0,0.022))
dev.off()

png(paste(wkdir, "/", header[4], "_histogram.png", sep=""))
par(mar=c(5,5,4,2) + 0.1)
hist(score_table$bind, breaks=seq(-500, 2500, by= 3), col='forestgreen', lty="blank",freq = FALSE, ylab = 'Density', xlab= "Rosetta 3.5 bind score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(-25,700), ylim=c(0, 0.18))
dev.off()

png(paste(wkdir, "/", header[5], "_histogram.png", sep=""))
par(mar=c(5,5,4,2) + 0.1)
hist(score_table$bind_pep_alone, breaks=seq(-200, 5000, by= 5), col='forestgreen', lty="blank", freq = FALSE, ylab = 'Density', xlab= "Rosetta 3.5 bind peptide alone score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(-100, 2500), ylim=c(0,0.07))
dev.off()

png(paste(wkdir, "/", header[6], "_histogram.png", sep=""))
par(mar=c(5,5,4,2) + 0.1)
hist(score_table$atlas, breaks=seq(-200, 300, by= 0.5), col='forestgreen', lty="blank", freq = FALSE, ylab = 'Density', xlab= "Rosetta 3.5 ATLAS score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(-50,75), ylim=c(0,0.4))
dev.off()

