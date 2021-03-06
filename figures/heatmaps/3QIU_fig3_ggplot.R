#!/usr/bin/env Rscript
library(argparse)
library(ggplot2)
library(reshape2)
library(grid)

parser <- ArgumentParser(description="Script to compute Figure 3 heatmaps for 3QIU pdb")
args <- parser$parse_args()
corrs <- c()
wkdir <- getwd()

# Get top most abundant peptides matrix
abundant_file <- paste("abundant_50_peptides_submatrix", sep="")
abundant_df <-read.table(paste(wkdir,"/", abundant_file, ".txt", sep = ""), header=FALSE, sep="\t")

# Format for ggplot heatmap
colnames(abundant_df) <- c(-3,-2,-1,1,2,3,4,5,6,7,8,9,10)
amino <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
abundant_df <- cbind(amino, abundant_df)
abundant_df.m <- melt(abundant_df)
anchors <- c("black", "red", "black", "red", rep("black", 7), "red", "red")

# WT reside subset
wt_df <- abundant_df.m[c(1,23,50,68,81,120,130,149,174,181,217,229,246),]
# Substitution subset
#sub_df <- abundant_df.m[c(161,105,86,59),] # old before relax
sub_df <- abundant_df.m[c(59,86,105,161,176),]


#png(paste(wkdir,"/", abundant_file, "_heatmap_ggplot.png", sep=""), width=2000, height=2000, res=300)
pdf(paste(wkdir,"/", abundant_file, "_heatmap_ggplot.pdf", sep=""), width=7, height=7)
print(ggplot(abundant_df.m, aes(variable, amino, fill=value)) 
	+ geom_tile()
	+ scale_fill_gradient(low="white", high="blue", breaks=c(0,1), lim=c(0,1)) 
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
colnames(score_df) <- c(-3,-2,-1,1,2,3,4,5,6,7,8,9,10)
amino <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
score_df <- cbind(amino, score_df)
score_df.m <- melt(score_df)
anchors <- c("black", "red", "black", "red", rep("black", 7), "red", "red")

# WT reside subset
wt_df <- score_df.m[c(1,23,50,68,81,120,130,149,174,181,217,229,246),]
# Substitution subset
#sub_df <- abundant_df.m[c(161,105,86,59),] # old before relax
sub_df <- score_df.m[c(59,86,105,161,176),]


#png(paste(wkdir,"/", score_file, "_heatmap_ggplot.png", sep=""), width=2000, height=2000, res=300)
pdf(paste(wkdir,"/", score_file, "_heatmap_ggplot.pdf", sep=""), width=7, height=7)
print(ggplot(score_df.m, aes(variable, amino, fill=value)) 
	+ geom_tile()
	+ scale_fill_gradient(low="white", high="blue", breaks=c(0,1), lim=c(0,1)) 
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
