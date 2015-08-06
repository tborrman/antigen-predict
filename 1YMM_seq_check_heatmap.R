#!/usr/bin/env Rscript
library(argparse)
library(gplots)

parser <- ArgumentParser(description="Make heatmap to check selection round sequence results for 1YMM against Birnbaum results")
parser$add_argument("-r", help="selection round", type="character", dest="round", required=TRUE)

args <- parser$parse_args()
wkdir <- getwd()
file <- "1YMM_sub_matrix_seq_check"

# Mark WT residues
markWT <- function(){
  WT_pos <- 1:14 
  WT_aa_pos <-c(17,9,8,3,3,14,16,16,12,9,13,3,4,8)
  # For now let's use .47 as the distance to fit rectangles
  fitter = .47
  rect(xleft=WT_pos - fitter,xright=WT_pos + fitter, ybottom=WT_aa_pos - fitter,
       ytop=WT_aa_pos + fitter,border="black", lwd=2)
}

# Mark Anchor residues
anchor_residues <- c("white", "white", "white", "white", "gray50","white","white",
                     "gray50","white","white","white","white","white", "white")

my_palette <- colorRampPalette(c("white", "mediumorchid4"))(n = 100)

dataframe <-read.table(paste(wkdir, "/", file, ".txt", sep = ""), header=FALSE, sep="\t")
matrix <- data.matrix(dataframe)
# Set row names and col names for axis reading
names_row <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
rownames(matrix) <- names_row
names_col <- c(-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10)
colnames(matrix) <- names_col

pdf(paste(wkdir, "/", file, "_heatmap.pdf", sep=""))
heatmap.2(matrix, Rowv = FALSE, Colv =FALSE, dendrogram = "none", 
          trace = "none", main=paste("Ob.1A12 TCR: ", args$round, " library", sep=""), xlab = "Peptide Position", ylab = "Amino Acid", 
          col = my_palette, cexRow = 1, cexCol = 1, add.expr = {markWT()}, ColSideColors=anchor_residues)
dev.off()
