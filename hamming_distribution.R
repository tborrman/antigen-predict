#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description="Output barplot of Hamming distance counts from WT peptide using preselection and round4 library")
parser$add_argument("-p", help="preselection score table output by hamming.py", type="character", required=TRUE)
parser$add_argument("-r", help="round 4 score table output by hamming.py", type="character", required=TRUE)
args <- parser$parse_args()


p_df <- read.table(args$p, header=TRUE, sep="\t")
r_df <- read.table(args$r, header=TRUE, sep="\t")
# head(p_df)
# head(r_df)
# print(p_df[p_df$peptide %in% r_df$peptide,])
# print(r_df[r_df$peptide %in% p_df$peptide,])

get_nchar <- function(d) {
  # Return number of characters of peptide used
  # in dataframe d
  return(nchar(as.character(d$peptide[1])))
}

# pep_length <- get_nchar(p_df)


get_hamming <- function(d) {
  # Combine and return hammings from preselection and round 4 libraries
  hamming <- c()
  for (i in 1:nrow(d)) {
    h <- c(d[i, "hamming.x"], d[i, "hamming.y"])
    if (all(is.na(h))) {
      print("ERROR: no hamming")
      stop()
    }
    else if (any(is.na(h))) {
      hamming <- c(hamming, h[!is.na(h)])
    }
    else {
      if (h[1] != h[2]) {
        print("ERROR: unequal hammings")
        stop()
      }
      else {
        hamming <- c(hamming, h[1])
      }
    }
  }
  return(hamming)
}

m_df <- merge(p_df, r_df, by="peptide", all=TRUE)

hamming <- get_hamming(m_df)


# print head(hamming)
# print(length(hamming))
# print(nrow(p_df) + nrow(r_df))
# print(nrow(m_df))
# print(m_df[m_df$peptide == "RDRKSGSAGMPRP",])

counts <- table(hamming)
png("hamming_distribution_pre_rnd4_zoom.png", height=2000, width=2500, res=300)
  barplot(counts, xlab= "Hamming distance from WT peptide", ylim=c(0,100))
dev.off()
png("hamming_distribution_pre_rnd4.png", height=2000, width=2500, res=300)
barplot(counts, xlab= "Hamming distance from WT peptide")
dev.off()

