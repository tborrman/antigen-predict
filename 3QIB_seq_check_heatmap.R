#!/usr/bin/env Rscript
library(argparse)
library(gplots)

parser <- ArgumentParser(description="Make heatmap to check selection round sequence results for 3QIB against Birnbaum results")
parser$add_argument("-r", help="selection round", type="character", dest="round", required=TRUE)

args <- parser$parse_args()
wkdir <- getwd()
file <- "3QIB_sub_matrix_seq_check"

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
anchor_residues <- c("gray50", "gray50", "white", "gray50", "white","white","white",
                     "white","white","white","white","gray50","gray50")

my_palette <- colorRampPalette(c("white", "red"))(n = 100)

dataframe <-read.table(paste(wkdir, "/", file, ".txt", sep = ""), header=FALSE, sep="\t")
matrix <- data.matrix(dataframe)
# Set row names and col names for axis reading
names_row <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
rownames(matrix) <- names_row
names_col <- c(-3,-2,-1,1,2,3,4,5,6,7,8,9,10)
colnames(matrix) <- names_col

pdf(paste(wkdir, "/", file, "_heatmap.pdf", sep=""))
heatmap.2(matrix, Rowv = FALSE, Colv =FALSE, dendrogram = "none", 
          trace = "none", main=paste("2B4 TCR: ", args$round, " library", sep=""), xlab = "Peptide Position", ylab = "Amino Acid", 
          col = my_palette, cexRow = 1, cexCol = 1, add.expr = {markWT()}, ColSideColors=anchor_residues)
dev.off()
