#!/usr/bin/env Rscript
library(argparse)
parser <- ArgumentParser(description="Make histograms for bind_pep_alone Rosetta 3.5 scoring results")

args <- parser$parse_args()
wkdir <- getwd()

png(paste(wkdir, "/TCR_histogram_rnds.png", sep=""), width=4000, height=3000, res=300)
par(mfrow=c(4,5), oma=c(6,8,5,2) + 0.1, mar=c(1,1,1,1) + 0.5, lwd=2)
# default: par(oma=c(2,2,2,2), c(5, 4, 4, 2) + 0.1,  c(3, 1, 0))

# 3QIB
score_table_0 <- read.table(paste(wkdir, "/3QIB/preselection/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_1 <- read.table(paste(wkdir, "/3QIB/round_1/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_2 <- read.table(paste(wkdir, "/3QIB/round_2/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_3 <- read.table(paste(wkdir, "/3QIB/round_3/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_4 <- read.table(paste(wkdir, "/3QIB/round_4/all/score_table.txt", sep=""), header=TRUE, sep="\t")

hist(score_table_0$bind_pep_alone, breaks=seq(-200,5000, by=50), col='red', lty="blank", xlab="",
	 	main="",xlim=c(-200,3500), cex.lab=2, cex.main=2, axes=FALSE)
par(xpd=NA)
title("Preselection", cex.main=2.5, line=1.5, font.main=1)
par(xpd=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_1$bind_pep_alone, breaks=seq(-200,5000, by=50), col='red', lty="blank", ylab="", xlab="",
		main="", xlim=c(-200,3500), cex.lab=2, cex.main=2, axes=FALSE)
par(xpd=NA)
title("Round 1", cex.main=2.5, line=1.5, font.main=1)
par(xpd=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_2$bind_pep_alone, breaks=seq(-200,5000, by=50), col='red', lty="blank", ylab="", xlab="",
		main="", xlim=c(-200,3500), cex.lab=2, cex.main=2, axes=FALSE)
par(xpd=NA)
title("Round 2", cex.main=2.5, line=1.5, font.main=1)
par(xpd=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_3$bind_pep_alone, breaks=seq(-200,5000, by=50), col='red', lty="blank", ylab="", xlab="",
		main="", xlim=c(-200,3500), cex.lab=2, cex.main=2, axes=FALSE)
par(xpd=NA)
title("Round 3", cex.main=2.5, line=1.5, font.main=1)
par(xpd=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_4$bind_pep_alone, breaks=seq(-200,5000, by=50), col='red', lty="blank", ylab="", xlab="",
		main="", xlim=c(-200,3500), cex.lab=2, cex.main=2, axes=FALSE)
par(xpd=NA)
title("Round 4", cex.main=2.5, line=1.5, font.main=1)
par(xpd=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")


# 3QIU
score_table_0 <- read.table(paste(wkdir, "/3QIU/preselection/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_1 <- read.table(paste(wkdir, "/3QIU/round_1/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_2 <- read.table(paste(wkdir, "/3QIU/round_2/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_3 <- read.table(paste(wkdir, "/3QIU/round_3/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_4 <- read.table(paste(wkdir, "/3QIU/round_4/all/score_table.txt", sep=""), header=TRUE, sep="\t")

hist(score_table_0$bind_pep_alone, breaks=seq(-200,5000, by=50), col='blue', lty="blank", xlab="",
		main="", xlim=c(-200,3500), cex.lab=2, cex.main=2, axes=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_1$bind_pep_alone, breaks=seq(-200, 5000, by= 50), col='blue', lty="blank", xlab="",
     main="", xlim=c(-200, 3500), cex.lab=2, cex.main=2,  axes=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_2$bind_pep_alone, breaks=seq(-200,5000, by=50), col='blue', lty="blank", ylab="", xlab="",
		main="", xlim=c(-200,3500), cex.lab=2, cex.main=2, axes=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_3$bind_pep_alone, breaks=seq(-200,5000, by=50), col='blue', lty="blank", ylab="", xlab="",
		main="", xlim=c(-200,3500), cex.lab=2, cex.main=2, axes=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_4$bind_pep_alone, breaks=seq(-200,5000, by=50), col='blue', lty="blank", ylab="", xlab="",
		main="", xlim=c(-200,3500), cex.lab=2, cex.main=2, axes=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")

# 4P2R
score_table_0 <- read.table(paste(wkdir, "/4P2R/preselection/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_1 <- read.table(paste(wkdir, "/4P2R/round_1/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_2 <- read.table(paste(wkdir, "/4P2R/round_2/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_3 <- read.table(paste(wkdir, "/4P2R/round_3/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_4 <- read.table(paste(wkdir, "/4P2R/round_4/all/score_table.txt", sep=""), header=TRUE, sep="\t")

hist(score_table_0$bind_pep_alone, breaks=seq(-1000,7000, by=50), col='forestgreen', lty="blank", xlab="",
		main="", xlim=c(-200,3500), cex.lab=2, cex.main=2, axes=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_1$bind_pep_alone, breaks=seq(-1000,7000, by=50), col='forestgreen', lty="blank", ylab="", xlab="",
		main="", xlim=c(-200,3500), cex.lab=2, cex.main=2, axes=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_2$bind_pep_alone, breaks=seq(-1000,7000, by=50), col='forestgreen', lty="blank", ylab="", xlab="",
		main="", xlim=c(-200,3500), cex.lab=2, cex.main=2, axes=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_3$bind_pep_alone, breaks=seq(-1000,7000, by=50), col='forestgreen', lty="blank", ylab="", xlab="",
		main="", xlim=c(-200,3500), cex.lab=2, cex.main=2, axes=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_4$bind_pep_alone, breaks=seq(-1000,7000, by=50), col='forestgreen', lty="blank", ylab="", xlab="",
		main="", xlim=c(-200,3500), cex.lab=2, cex.main=2, axes=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")

# 1YMM
score_table_0 <- read.table(paste(wkdir, "/1YMM/preselection/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_1 <- read.table(paste(wkdir, "/1YMM/round_1/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_2 <- read.table(paste(wkdir, "/1YMM/round_2/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_3 <- read.table(paste(wkdir, "/1YMM/round_3/all/score_table.txt", sep=""), header=TRUE, sep="\t")
score_table_4 <- read.table(paste(wkdir, "/1YMM/round_4/all/score_table.txt", sep=""), header=TRUE, sep="\t")

hist(score_table_0$bind_pep_alone, breaks=seq(-500, 4000, by=30), col='mediumorchid4', lty="blank", xlab="",
		main="", xlim=c(-100,2200), cex.lab=2, cex.main=2, axes=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_1$bind_pep_alone, breaks=seq(-500, 4000, by=30), col='mediumorchid4', lty="blank", ylab="", xlab="",
		main="", xlim=c(-100,2200), cex.lab=2, cex.main=2, axes=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_2$bind_pep_alone, breaks=seq(-500, 4000, by=30), col='mediumorchid4', lty="blank", ylab="", xlab="",
		main="", xlim=c(-100,2200), cex.lab=2, cex.main=2, axes=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_3$bind_pep_alone, breaks=seq(-500, 4000, by=30), col='mediumorchid4', lty="blank", ylab="", xlab="",
		main="", xlim=c(-100,2200), cex.lab=2, cex.main=2, axes=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")
hist(score_table_4$bind_pep_alone, breaks=seq(-500, 4000, by=30), col='mediumorchid4', lty="blank", ylab="", xlab="",
		main="", xlim=c(-100,2200), cex.lab=2, cex.main=2, axes=FALSE)
axis(1, lwd=2, cex.axis=1.5) 
axis(2, lwd=2, cex.axis=1.5)
box(bty="l")


title(xlab=expression(paste(Delta, "G", ""["BIND"], sep= "")), ylab="Number of peptides", outer=TRUE, cex.lab=2.5, line=3)
dev.off()