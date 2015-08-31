#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description="Make histograms for Rosetta 3.5 scoring results")
parser$add_argument("-in", help="score table output by aggregate_scores.py (score_table.txt)", type="character", dest="infile", required=TRUE)
parser$add_argument("-r", help="selection round (ex. 3)", type="integer", dest="round", required=TRUE)

args <- parser$parse_args()
wkdir <- getwd()

# Top 5 most abundant peptides for 3QIB
top5_peptides <- c("ADLVAFFKEASKR", "ATHVAFLKAATKK", "AAQVAFLKAATKA", "ATHVAFLKAATKA", "AAQVAFLKAATKK")

score_table <-read.table(paste(wkdir, "/", args$infile, sep=""), header=TRUE, sep="\t")
header <- colnames(score_table)

pep1 <- score_table[which(score_table$peptide == top5_peptides[1]),]
pep2 <- score_table[which(score_table$peptide == top5_peptides[2]),]
pep3 <- score_table[which(score_table$peptide == top5_peptides[3]),]
pep4 <- score_table[which(score_table$peptide == top5_peptides[4]),]
pep5 <- score_table[which(score_table$peptide == top5_peptides[5]),]


colors <- c("red", "orange", "green", "blue", "purple")

png(paste(wkdir, "/", header[3], "_histogram_rank_zoom.png", sep=""))
par(mar=c(5,5,4,2) + 0.1)
hist(score_table$standard, breaks=2000, col='lightpink',  lty="blank", ylab = 'Frequency', xlab= "Rosetta 3.5 standard score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(120,200), ylim=c(0,400))
abline(v=pep1$standard, col=colors[1], lwd=1.5)
abline(v=pep2$standard, col=colors[2], lwd=1.5)
abline(v=pep3$standard, col=colors[3], lwd=1.5)
abline(v=pep4$standard, col=colors[4], lwd=1.5)
abline(v=pep5$standard, col=colors[5], lwd=1.5)
# Get score ranking for legend
ranks = c()
standard_order <- score_table[order(score_table$standard),]
for (peptide in top5_peptides) {
  ranks <- c(ranks,which(standard_order$peptide == peptide))
}
par(family="mono")
legend(x="topright", paste(top5_peptides, ranks, sep="  "), col=colors, lty= rep(1,5) , lwd= rep(1.5,5), bty='n')
dev.off()

png(paste(wkdir, "/", header[4], "_histogram_rank_zoom.png", sep=""))
par(mar=c(5,5,4,2) + 0.1)
hist(score_table$bind, breaks=3500, col='lightpink',  lty="blank", ylab = 'Frequency', xlab= "Rosetta 3.5 bind score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(-25, -10), ylim=c(0,400))
abline(v=pep1$bind, col=colors[1], lwd=1.5)
abline(v=pep2$bind, col=colors[2], lwd=1.5)
abline(v=pep3$bind, col=colors[3], lwd=1.5)
abline(v=pep4$bind, col=colors[4], lwd=1.5)
abline(v=pep5$bind, col=colors[5], lwd=1.5)
# Get score ranking for legend
ranks = c()
bind_order <- score_table[order(score_table$bind),]
for (peptide in top5_peptides) {
  ranks <- c(ranks,which(bind_order$peptide == peptide))
}
par(family="mono")
legend(x="topright", paste(top5_peptides, ranks, sep="  "), col=colors, lty= rep(1,5) , lwd= rep(1.5,5), bty='n')
dev.off()

png(paste(wkdir, "/", header[5], "_histogram_rank_zoom.png", sep=""))
par(mar=c(5,5,4,2) + 0.1)
hist(score_table$bind_pep_alone, breaks=2000, col='lightpink',  lty="blank", ylab = 'Frequency', xlab= "Rosetta 3.5 bind peptide alone score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(-40,20), ylim=c(0,400))
abline(v=pep1$bind_pep_alone, col=colors[1], lwd=1.5)
abline(v=pep2$bind_pep_alone, col=colors[2], lwd=1.5)
abline(v=pep3$bind_pep_alone, col=colors[3], lwd=1.5)
abline(v=pep4$bind_pep_alone, col=colors[4], lwd=1.5)
abline(v=pep5$bind_pep_alone, col=colors[5], lwd=1.5)
# Get score ranking for legend
ranks = c()
bind_pep_alone_order <- score_table[order(score_table$bind_pep_alone),]
for (peptide in top5_peptides) {
  ranks <- c(ranks,which(bind_pep_alone_order$peptide == peptide))
}
par(family="mono")
legend(x="topright", paste(top5_peptides, ranks, sep="  "), col=colors, lty= rep(1,5) , lwd= rep(1.5,5), bty='n')
dev.off()

png(paste(wkdir, "/", header[6], "_histogram_rank_zoom.png", sep=""))
par(mar=c(5,5,4,2) + 0.1)
hist(score_table$atlas, breaks=300, col='lightpink',  lty="blank", ylab = 'Frequency', xlab= "Rosetta 3.5 ATLAS score",cex.main=1.8, cex.lab=1.5, 
     main = paste("Round ", args$r, sep=""), xlim=c(-65,-35), ylim=c(0,400))
abline(v=pep1$atlas, col=colors[1], lwd=1.5)
abline(v=pep2$atlas, col=colors[2], lwd=1.5)
abline(v=pep3$atlas, col=colors[3], lwd=1.5)
abline(v=pep4$atlas, col=colors[4], lwd=1.5)
abline(v=pep5$atlas, col=colors[5], lwd=1.5)
# Get score ranking for legend
ranks = c()
atlas_order <- score_table[order(score_table$atlas),]
for (peptide in top5_peptides) {
  ranks <- c(ranks,which(atlas_order$peptide == peptide))
}
par(family="mono")
legend(x="topright", paste(top5_peptides, ranks, sep="  "), col=colors, lty= rep(1,5) , lwd= rep(1.5,5), bty='n')
dev.off()
