#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description="Plot scores vs SPR energies")
parser$add_argument("-pdb", help= "pdb ID (3QIB)", type="character", dest="pdb", required=TRUE)

args <- parser$parse_args()
wkdir <- getwd()

# Open files
dG_df <- read.table(paste(wkdir, "/", args$pdb, "_dGs.txt", sep=""), header=TRUE, sep="\t")
scores_df <- read.table(paste(wkdir, "/peptides_scores.txt", sep=""), header=TRUE, sep="\t")

scores <- scores_df$score
dGs <- dG_df$dG
peptides <- scores_df$id

pdf(paste(wkdir,"/", args$pdb, "_SPR_plot.pdf", sep=""))
	plot(dGs, scores, main= args$pdb, xlab= expression(paste("Experimentally determined ",Delta, "G (kcal/mol)", sep="")),
	ylab= "Rosetta score", pch=16, col= "blue", xlim=c(-9,-4), ylim=c(-60,100))
	text(dGs, scores, peptides, cex=0.7, pos=3)
dev.off()

# Print correlation
print(paste("Correlation: ", cor(dGs, scores), sep=""))
