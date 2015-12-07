#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description="Plot scores vs SPR energies")
#parser$add_argument("-pdb", help= "pdb ID (3QIB)", type="character", dest="pdb", required=TRUE)

args <- parser$parse_args()
wkdir <- getwd()

# Open files
dG_2b4_df <- read.table(paste(wkdir, "/3QIB/3QIB_dGs.txt", sep=""), header=TRUE, sep="\t")
dG_226_df <- read.table(paste(wkdir, "/3QIU/3QIU_dGs.txt", sep=""), header=TRUE, sep="\t")
dG_5cc7_df <- read.table(paste(wkdir, "/4P2R/4P2R_dGs.txt", sep=""), header=TRUE, sep="\t")
scores_2b4_df <- read.table(paste(wkdir, "/3QIB/peptides_scores.txt", sep=""), header=TRUE, sep="\t")
scores_226_df <- read.table(paste(wkdir, "/3QIU/peptides_scores.txt", sep=""), header=TRUE, sep="\t")
scores_5cc7_df <- read.table(paste(wkdir, "/4P2R/peptides_scores.txt", sep=""), header=TRUE, sep="\t")
# Collecct dGs
dGs_2b4 <- dG_2b4_df$dG
dGs_226 <- dG_226_df$dG
dGs_5cc7 <- dG_5cc7_df$dG
# Collect scores
scores_2b4 <- scores_2b4_df$score
scores_226 <- scores_226_df$score
scores_5cc7 <- scores_5cc7_df$score
# Aggregate
dGs <- c(dGs_2b4, dGs_226, dGs_5cc7)
scores <- c(scores_2b4, scores_226, scores_5cc7)
#peptides <- scores_df$id

png(paste(wkdir,"/SPR_plot.png", sep=""), width=2000, height=2000, res = 300)
	par(mar=c(5, 5, 5, 2) + 0.1, mgp=c(3,1,0),lwd=2)
	plot(dGs_2b4, scores_2b4, xlab= expression(paste("Experimentally determined ",Delta, "G (kcal/mol)", sep="")),
	cex.lab=1.5, cex.main=2.5, cex=1.5, ylab= expression(paste(Delta, "G", ""["BIND"], sep= "")), 
	pch=21, col= "black", bg="red", xlim=c(-9.5,-4.5), ylim=c(-60,100),axes=FALSE)
		points(dGs_226, scores_226, pch=22, col= "black", bg="blue",cex=1.5)
		points(dGs_5cc7, scores_5cc7, pch=24, col= "black", bg="forestgreen",cex=1.5)
	axis(1, lwd=2, cex.axis=1.3) 
	axis(2, lwd=2, cex.axis=1.3)
	box(bty="l")

	#text(dGs, scores, peptides, cex=0.7, pos=3)
	#text(-8.5, 75, paste("r = ", round(cor(dGs, scores), digits=2), sep=""), cex=1.7, pos=4)
	legend(-9.5,100, c("2B4", "226", "5cc7"), col=rep("black",3), pt.bg=c("red", "blue", "forestgreen"), pch=c(21,22,24), cex=1.5)
dev.off()

# Print correlation
print(paste("Correlation: ", cor(dGs, scores), sep=""))
