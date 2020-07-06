#!/usr/bin/env Rscript
library(argparse)
library(ggplot2)
library(reshape2)
library(grid)

parser <- ArgumentParser(description="Script to compute Figure 3 heatmaps for 1YMM pdb")
args <- parser$parse_args()
corrs <- c()
wkdir <- getwd()

# Get top most abundant peptides matrix
abundant_file <- paste("abundant_50_peptides_submatrix", sep="")
abundant_df <-read.table(paste(wkdir,"/", abundant_file, ".txt", sep = ""), header=FALSE, sep="\t")

# Format for ggplot heatmap
colnames(abundant_df) <- c(-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10)
amino <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
abundant_df <- cbind(amino, abundant_df)
abundant_df.m <- melt(abundant_df)
anchors <- c(rep("black", 4), "red", rep("black",2), "red", rep("black", 6))

# WT reside subset
wt_df <- abundant_df.m[c(4,32,53,78,98,107,125,145,169,192,208,238,257,273),]
# Substitution subset
#sub_df <- abundant_df.m[c(1,181,242,24,105,246,27,228,54,214,55,175,275),] # before relax
sub_df <- abundant_df.m[c(27, 54, 181, 214, 228, 242, 246, 275),] 

#png(paste(wkdir,"/", abundant_file, "_heatmap_ggplot.png", sep=""), width=2000, height=2000, res=300)
pdf(paste(wkdir,"/", abundant_file, "_heatmap_ggplot.pdf", sep=""), width=7, height=7)
print(ggplot(abundant_df.m, aes(variable, amino, fill=value)) 
	+ geom_tile()
	+ scale_fill_gradient(low="white", high="mediumorchid4", breaks=c(0,1), lim=c(0,1)) 
	+ scale_y_discrete(limits=rev(levels(abundant_df.m$amino))) 
	+ labs(x="Peptide position", y= "Amino acid", title="Top 50 most abundant peptides", fill="")
	+ theme(axis.text.y=element_text(size=15, colour="black"), axis.title.y=element_text(size=20, vjust=.2),
		axis.text.x=element_text(size=15, colour=anchors), axis.title.x=element_text(size=20, vjust=-1),
		axis.ticks = element_blank(), panel.border = element_rect(colour="black", fill=NA, size=2),
		plot.title = element_text(size=20, face="bold", vjust=2),
		plot.margin = unit(c(1,1,1,1), "cm"), # top, right, bottom, left
		legend.text = element_text(size=15), legend.key.height=unit(2.6,"cm"))
	+ geom_tile(data=wt_df, aes(variable, amino), colour="black", size=1)
	+ geom_text(data=sub_df, aes(variable, amino, label=sprintf("%.2f", value)), size=4)
	)
dev.off()

# Get top scoring peptides matrix
score_file <- paste("bind_pep_alone_score_50_peptides_submatrix", sep="")
score_df <-read.table(paste(wkdir,"/", score_file, ".txt", sep = ""), header=FALSE, sep="\t")

# Format for ggplot heatmap
colnames(score_df) <- c(-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10)
amino <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
score_df <- cbind(amino, score_df)
score_df.m <- melt(score_df)
anchors <- c(rep("black", 4), "red", rep("black",2), "red", rep("black", 6))

# WT reside subset
wt_df <- score_df.m[c(4,32,53,78,98,107,125,145,169,192,208,238,257,273),]
# Substitution subset
# sub_df <- score_df.m[c(1,181,242,24,105,246,27,228,54,214,55,175,275),] before relax
sub_df <- score_df.m[c(27, 54, 181, 214, 228, 242, 246, 275),] 

#png(paste(wkdir,"/", score_file, "_heatmap_ggplot.png", sep=""), width=2000, height=2000, res=300)
pdf(paste(wkdir,"/", score_file, "_heatmap_ggplot.pdf", sep=""), width=7, height=7)
print(ggplot(score_df.m, aes(variable, amino, fill=value)) 
	+ geom_tile()
	+ scale_fill_gradient(low="white", high="mediumorchid4", breaks=c(0,1), lim=c(0,1)) 
	+ scale_y_discrete(limits=rev(levels(score_df.m$amino))) 
	+ labs(x="Peptide position", y= "Amino acid", title="Top 50 best scoring peptides", fill="")
	+ theme(axis.text.y=element_text(size=15, colour="black"), axis.title.y=element_text(size=20, vjust=.2),
		axis.text.x=element_text(size=15, colour=anchors), axis.title.x=element_text(size=20, vjust=-1),
		axis.ticks = element_blank(), panel.border = element_rect(colour="black", fill=NA, size=2),
		plot.title = element_text(size=20, face="bold", vjust=2),
		plot.margin = unit(c(1,1,1,1), "cm"), # top, right, bottom, left
		legend.text = element_text(size=15), legend.key.height=unit(2.6,"cm"))
	+ geom_tile(data=wt_df, aes(variable, amino), colour="black", size=1)
	+ geom_text(data=sub_df, aes(variable, amino, label=sprintf("%.2f", value)), size=4)
	)
dev.off()