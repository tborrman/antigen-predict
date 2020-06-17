#!/usr/bin/env Rscript
library(argparse)
library(ROCR)

parser <- ArgumentParser(description=paste("Produce ROC curves and calculate AUC for input TCR",
                                           "should be run in TCR/ folder", sep=" "))
parser$add_argument("-t", help="TCR (ex. 1YMM)", type="character", dest="t", required=TRUE)
parser$add_argument("-x", help="filter positives for peptides with read counts in the top x percentile", 
                    type="double", dest="x", default=0.95)

args <- parser$parse_args()
wkdir <- getwd()

plot_ROC <- function(scores, labels, TCR, path, x) {
  # Plot ROC curves for data
  # Args:
  #   scores: dataframe with scores and separate
  #           scoring methods for each column
  #   labels: dataframe with labels 
  #           (same dimensions as scores)
  #   TCR: character string of TCR
  #   path: character string full path to TCR/ directory
  #   x: numeric cutoff for read counts to identify as a positive
  mycolors <- c("dodgerblue", "green3", "red")
  pred <- prediction(scores, labels)
  perf <- performance(pred,"tpr","fpr")
  png(paste(path, "/", TCR, "/AUC/5c1/ROC_AUC_", TCR,"_", x, "_percentile_cutoff_blosum_average_rank.png", sep=""), 
      height=1800, width=1800, res=300)
  plot(perf, col=as.list(mycolors), main=TCR, lwd=2)
  abline(a=0, b=1, col="black", lty=2)
  auc.perf <- performance(pred, measure = "auc")
  legend("bottomright",paste(colnames(scores), rep(": ", 3),
                             round(as.numeric(auc.perf@y.values),2), sep=""), 
         col=mycolors, title="AUC", lty=1, lwd=2)
  dev.off()
  return()
}

plot_PR <- function(scores, labels, TCR, path, x) {
  # Plot ROC curves for data
  # Args:
  #   scores: dataframe with scores and separate
  #           scoring methods for each column
  #   labels: dataframe with labels 
  #           (same dimensions as scores)
  #   TCR: character string of TCR
  #   path: character string full path to TCR/ directory
  #   x: numeric cutoff for read counts to identify as a positive
  mycolors <- c("dodgerblue", "green3", "red")
  pred <- prediction(scores, labels)
  perf <- performance(pred, "prec", "rec")
  png(paste(path, "/", TCR, "/AUC/5c1/PR_AUC_", TCR,"_", x, "_percentile_cutoff_blosum_average_rank.png", sep=""), 
      height=1800, width=1800, res=300)
  plot(perf, col=as.list(mycolors), main=TCR, lwd=2)
  #abline(a=0, b=1, col="black", lty=2)
  auc.perf <- performance(pred, measure = "aucpr")
  legend("topright",paste(colnames(scores), rep(": ", 3),
                          round(as.numeric(auc.perf@y.values),2), sep=""), 
         col=mycolors, title="AUC", lty=1, lwd=2)
  dev.off()
  return()
}

if (args$t == "1YMM") {
  negatives <- read.table(paste(wkdir, "/", args$t, 
                                "/AUC/score_table_hamming_blosum_ENPVVHFFKNIVTP_preselection_remove_rnd4.txt",
                                sep=""), header=TRUE, sep="\t")
  positives <- read.table(paste(wkdir, "/", args$t,
                                "/AUC/score_table_hamming_blosum_ENPVVHFFKNIVTP_round4.txt",
                                sep=""), header=TRUE, sep="\t") 
} else {
  negatives <- read.table(paste(wkdir, "/", args$t, 
                                "/AUC/5c1/score_table_hamming_blosum_ANGVAFFLTPFKA_preselection_remove_rnd4.txt",
                                sep=""), header=TRUE, sep="\t")
  positives <- read.table(paste(wkdir, "/", args$t,
                                "/AUC/5c1/score_table_hamming_blosum_ANGVAFFLTPFKA_round4.txt",
                                sep=""), header=TRUE, sep="\t")
}

# add labels
labels <- rep(FALSE, nrow(negatives))
negatives <- cbind(negatives, labels)
labels <- rep(TRUE, nrow(positives))
positives <- cbind(positives, labels)

# Filter positives
q <- as.numeric(quantile(positives$reads, args$x))
positives <- positives[positives$reads >= q,]

# Merge into one table
df <- rbind(negatives, positives)


# Flip sign of scores for ROC predictions
scores <- df[c("standard", "bind", "bind_pep_alone", "hamming")]
flip_scores <- data.frame(apply(scores, 2, function(x) x*(-1) ))
all_scores <- cbind(flip_scores, df["blosum62"])

# Rank
# Largest number has largest rank
# Here, larger rank implies higher probability to be a true binder
bind_pep_alone_rank <- rank(all_scores$bind_pep_alone)
blosum_rank <- rank(all_scores$blosum62)

rank_df <- data.frame(bind_pep_alone_rank, blosum_rank)
average_rank <- apply(rank_df, 1, mean)
rank_df <- cbind(all_scores["bind_pep_alone"], 
                 all_scores["blosum62"],average_rank)

# Make labels dataframe with same dimensions as rank_df
df_lab <- data.frame(df$labels, df$labels, df$labels)

# Plot ROC
plot_ROC(rank_df, df_lab, args$t, wkdir, args$x)
plot_PR(rank_df, df_lab, args$t, wkdir, args$x)



