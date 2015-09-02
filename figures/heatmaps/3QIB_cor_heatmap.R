#!/usr/bin/env Rscript
library(argparse)
library(gplots)

parser <- ArgumentParser(description="Script to compute cross-reactivity heatmaps and correlation plots on 3QIB pdb")
parser$add_argument("-s", help="scoring function (ex. standard, bind, bind_pep_alone, atlas)", type="character", dest="score", required=TRUE)
args <- parser$parse_args()
corrs <- c()
wkdir <- getwd()

# Compute correlation across positions excluding the anchor positions
my_corr <- function(abundant_df, score_df) {
  # Remove anchor columns
  clean_abundant_df <- abundant_df[-c(1,2,4,12,13)]
  clean_score_df <- score_df[-c(1,2,4,12,13)]
  return(cor(as.vector(as.matrix(clean_abundant_df)),as.vector(as.matrix(clean_score_df))))
}

# Mark WT residues
markWT <- function(){
  WT_pos <- 1:13 
  WT_aa_pos <-c(20,18,11,13,20,1,11,12,7,20,4,12,15)
  # For now let's use .47 as the distance to fit rectangles
  fitter = .47
  rect(xleft=WT_pos - fitter,xright=WT_pos + fitter, ybottom=WT_aa_pos - fitter,
       ytop=WT_aa_pos + fitter,border="black", lwd=2)
}

# Mark Anchor residues
anchor_residues <- c("gray50", "gray50", "white", "gray50", "white","white","white","white","white","white","white","gray50","gray50")

peptide_numbers <- 1:1000

for (i in peptide_numbers){
  # Get top most abundant peptides matrix
  abundant_file <- paste("abundant_", i, "_peptides_submatrix", sep="")
  abundant_df <-read.table(paste(wkdir,"/", abundant_file, ".txt", sep = ""), header=FALSE, sep="\t")
  # Get top scoring peptides matrix
  score_file <- paste(args$score, "_score_", i, "_peptides_submatrix", sep="")
  score_df <-read.table(paste(wkdir,"/", score_file, ".txt", sep = ""), header=FALSE, sep="\t")
  corrs <- c(corrs, my_corr(abundant_df, score_df))
  
  
  abundant_matrix <- data.matrix(abundant_df)
  score_matrix <- data.matrix(score_df)
  
  my_palette <- colorRampPalette(c("white", "red"))(n = 100);
  # Set row names and col names for axis reading
  names_row <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  rownames(abundant_matrix) <- names_row
  rownames(score_matrix) <- names_row
  names_col <- c(-3,-2,-1,1,2,3,4,5,6,7,8,9,10)
  colnames(abundant_matrix) <- names_col
  colnames(score_matrix) <- names_col
  
  
  pdf(paste(wkdir,"/", abundant_file, "_heatmap.pdf", sep=""))
  heatmap.2(abundant_matrix, Rowv = FALSE, Colv =FALSE, dendrogram = "none", 
            trace = "none",main= paste('2B4 TCR','\n', 'Top ', i, ' most abundant peptides', sep=""), 
            xlab = "Peptide Position", ylab = "Amino Acid", col = my_palette, cexRow = 1, cexCol = 1,
            add.expr = {markWT()}, ColSideColors=anchor_residues)
  dev.off();
  
  pdf(paste(wkdir,"/", score_file, "_heatmap.pdf", sep=""))
  heatmap.2(score_matrix, Rowv = FALSE, Colv =FALSE, dendrogram = "none", 
            trace = "none",main= paste('2B4 TCR','\n', 'Top ', i ,' scoring peptides', sep=""), 
            xlab = "Peptide Position", ylab = "Amino Acid", col = my_palette, cexRow = 1, cexCol = 1,
            add.expr = {markWT()}, ColSideColors=anchor_residues)
  dev.off();
  
}
pdf (paste(wkdir, "/", args$score, "_corr.pdf", sep=""))
  par(mar=c(5,5,4,2) + 0.1)
  plot(peptide_numbers, corrs, main= '2B4', cex.main = 1.8, cex.lab= 1.5, xlab= 'Number of peptides', ylim=c(-0.2, 1.0), ylab= "Pearson's r", pch=20, type='o', col='red')
dev.off();

# RESULTS
max_index <- which(corrs == max(corrs))
print(paste("Max corelation:", max(corrs), sep=" "))
print(paste("Number of peptides at max correlation:", peptide_numbers[max_index],sep=" "))
