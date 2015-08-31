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
hist(score_table$standard, breaks=seq(0, 8000, by= 10), col='lightpink', lty="blank", freq = FALSE, ylab = 'Density', xlab= "Rosetta 3.5 standard score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(0,4000), ylim=c(0,0.022))
dev.off()

png(paste(wkdir, "/", header[4], "_histogram.png", sep=""))
par(mar=c(5,5,4,2) + 0.1)
hist(score_table$bind, breaks=seq(-200, 2000, by= 5), col='lightpink', lty="blank",freq = FALSE, ylab = 'Density', xlab= "Rosetta 3.5 bind score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(-200,2000), ylim=c(0, 0.09))
dev.off()

png(paste(wkdir, "/", header[5], "_histogram.png", sep=""))
par(mar=c(5,5,4,2) + 0.1)
hist(score_table$bind_pep_alone, breaks=seq(-200, 5000, by= 5), col='lightpink', lty="blank", freq = FALSE, ylab = 'Density', xlab= "Rosetta 3.5 bind peptide alone score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(-200,3500), ylim=c(0,0.04))
dev.off()

png(paste(wkdir, "/", header[6], "_histogram.png", sep=""))
par(mar=c(5,5,4,2) + 0.1)
hist(score_table$atlas, breaks=seq(-100, 200, by= 1), col='lightpink', lty="blank", freq = FALSE, ylab = 'Density', xlab= "Rosetta 3.5 ATLAS score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(-90,200), ylim=c(0,0.2))
dev.off()

