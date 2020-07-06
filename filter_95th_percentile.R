#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description=paste("Filter score table for peptides with read counts",
                                           " in the 95th percentile ", sep=" "))
parser$add_argument("-in", help="score table (ex. score_table.txt)", type="character", dest="i", required=TRUE)


args <- parser$parse_args()
wkdir <- getwd()

df <- read.table(paste(wkdir, "/", args$i, sep=""), header=TRUE, sep="\t") 
# Filter positives
q <- as.numeric(quantile(df$reads, 0.95))
positives <- df[df$reads >= q,]

write.table(positives, "score_table_reads_95th_percentile.txt", quote=FALSE, sep="\t", row.names=FALSE)
