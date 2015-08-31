#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description="Make histograms for Rosetta 3.5 scoring results")
parser$add_argument("-in", help="score table output by aggregate_scores.py (score_table.txt)", type="character", dest="infile", required=TRUE)
parser$add_argument("-r", help="selection round (ex. 3)", type="integer", dest="round", required=TRUE)

args <- parser$parse_args()
wkdir <- getwd()

# Top 5 most abundant peptides for 1YMM
top5_peptides <- c("AHRKIHFFKAQIGR", "WAQNLHFFASMPCL", "NVPVLHFYRGLIVS", "MEAVIHFYRGLQCR", "RGISVVHFWRSELC" )

score_table <-read.table(paste(wkdir, "/", args$infile, sep=""), header=TRUE, sep="\t")
header <- colnames(score_table)

pep1 <- score_table[which(score_table$peptide == top5_peptides[1]),]
pep2 <- score_table[which(score_table$peptide == top5_peptides[2]),]
pep3 <- score_table[which(score_table$peptide == top5_peptides[3]),]
pep4 <- score_table[which(score_table$peptide == top5_peptides[4]),]
pep5 <- score_table[which(score_table$peptide == top5_peptides[5]),]

colors <- c("red", "orange", "green", "blue", "purple")

png(paste(wkdir, "/", header[3], "_histogram_rank.png", sep=""))
par(mar=c(5,5,4,2) + 0.1)

hist(score_table$standard, breaks=seq(7200, 11000, by= 10), col='aquamarine', lty="blank", freq = FALSE, ylab = 'Density', xlab= "Rosetta 3.5 standard score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(7200,10000), ylim=c(0,0.015))

abline(v=pep1$standard, col=colors[1], lwd=1.5)
abline(v=pep2$standard, col=colors[2], lwd=1.5)
abline(v=pep3$standard, col=colors[3], lwd=1.5)
abline(v=pep4$standard, col=colors[4], lwd=1.5)
abline(v=pep5$standard, col=colors[5], lwd=1.5)

legend(x="topright", top5_peptides, col=colors, lty= rep(1,5) , lwd= rep(1.5,5), bty='n')

dev.off()

png(paste(wkdir, "/", header[4], "_histogram_rank.png", sep=""))
par(mar=c(5,5,4,2) + 0.1)

hist(score_table$bind, breaks=seq(-10, 100, by= 0.5), col='aquamarine', lty="blank", freq = FALSE, ylab = 'Density', xlab= "Rosetta 3.5 bind score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(-10,90), ylim=c(0, 0.19))

abline(v=pep1$bind, col=colors[1], lwd=1.5)
abline(v=pep2$bind, col=colors[2], lwd=1.5)
abline(v=pep3$bind, col=colors[3], lwd=1.5)
abline(v=pep4$bind, col=colors[4], lwd=1.5)
abline(v=pep5$bind, col=colors[5], lwd=1.5)

legend(x="topright", top5_peptides, col=colors, lty= rep(1,5) , lwd= rep(1.5,5), bty='n')

dev.off()

png(paste(wkdir, "/", header[5], "_histogram_rank.png", sep=""))
par(mar=c(5,5,4,2) + 0.1)

hist(score_table$bind_pep_alone, breaks=seq(-500, 4000, by= 5), col='aquamarine', lty="blank", freq = FALSE, ylab = 'Density', xlab= "Rosetta 3.5 bind peptide alone score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(-100,2200), ylim=c(0,0.019))

abline(v=pep1$bind_pep_alone, col=colors[1], lwd=1.5)
abline(v=pep2$bind_pep_alone, col=colors[2], lwd=1.5)
abline(v=pep3$bind_pep_alone, col=colors[3], lwd=1.5)
abline(v=pep4$bind_pep_alone, col=colors[4], lwd=1.5)
abline(v=pep5$bind_pep_alone, col=colors[5], lwd=1.5)

legend(x="topright", top5_peptides, col=colors, lty= rep(1,5) , lwd= rep(1.5,5), bty='n')

dev.off()

png(paste(wkdir, "/", header[6], "_histogram_rank.png", sep=""))
par(mar=c(5,5,4,2) + 0.1)

hist(score_table$atlas, breaks=seq(-70, 0, by= 0.2), col='aquamarine', lty="blank", freq = FALSE, ylab = 'Density', xlab= "Rosetta 3.5 ATLAS score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(-55, -10), ylim=c(0,0.15))

abline(v=pep1$atlas, col=colors[1], lwd=1.5)
abline(v=pep2$atlas, col=colors[2], lwd=1.5)
abline(v=pep3$atlas, col=colors[3], lwd=1.5)
abline(v=pep4$atlas, col=colors[4], lwd=1.5)
abline(v=pep5$atlas, col=colors[5], lwd=1.5)

legend(x="topright", top5_peptides, col=colors, lty= rep(1,5) , lwd= rep(1.5,5), bty='n')

dev.off()

