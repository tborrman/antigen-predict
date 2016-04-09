# Make histogram of peptide scores from NCBI's nr database (not complete)
library(ggplot2)
peptides <- read.table("C:/cygwin64/home/Tyler/Research/TCR/nr/birnbaum/hist_test/peptide_histogram.txt", sep="\t", header=FALSE)
colnames(peptides) <- c("dGbind", "freq")

peptides_sub <- peptides[(peptides$dGbind >= -100) & (peptides$dGbind <= 1500),]

#peptides_sub <- peptides[(peptides$dGbind >= -50) & (peptides$dGbind <= 10 ),]

# birnbaum
birnbaum <- read.table("C:/cygwin64/home/Tyler/Research/TCR/nr/birnbaum/birnbaum_nr_scores.txt", sep="\t", header=FALSE)
colnames(birnbaum) <- c("peptide", "dGbind")

sorted_birnbaum <- birnbaum[order(birnbaum$dGbind),]
top_wt_peptide <- c("GAHCIHFFKSAVCR","ENPVVHFFKINVTP")  
top_wt_score <- c(-32.5930,-4.1420)


# barplot(peptides_sub$freq, col="mediumorchid4", space=0, border=NA, names.arg = peptides_sub$peptide)

png("C:/cygwin64/home/Tyler/Research/TCR/nr/birnbaum/peptide_historgram.png", height=3000, width=4000, res=400)
p1 <- (ggplot(peptides_sub, aes(x= dGbind, y=freq)) 
 + geom_bar(stat='identity', colour="mediumorchid4") 
 + geom_vline(xintercept = top_wt_score, colour=c("red","green"), size = c(1,1))
 + labs(x= expression(paste(Delta, "G", ""["BIND"], sep= "")), y = "Frequency", size=12 )
        
 + theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour="black", size=0.5),
          axis.text.x = element_text(colour="black", size=14),
          axis.text.y = element_text(colour="black", size=14),
          axis.title = element_text(colour="black", size=20))
)
p1
dev.off()
