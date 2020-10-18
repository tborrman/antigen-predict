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
  mycolors <- c("gray", "red")
  pred <- prediction(scores, labels)
  perf <- performance(pred,"tpr","fpr")
  png(paste(path, "/AUC/ROC_AUC_modeller", TCR,"_", x, "_percentile_cutoff.png", sep=""), 
      height=1800, width=1800, res=300)
  plot(perf, col=as.list(mycolors), main=TCR, lwd=2)
  abline(a=0, b=1, col="black", lty=2)
  auc.perf <- performance(pred, measure = "auc")
  legend("bottomright",paste(colnames(scores), rep(": ", 2),
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
  mycolors <- c("gray", "red")
  pred <- prediction(scores, labels)
  perf <- performance(pred, "prec", "rec")
  png(paste(path, "/AUC/PR_AUC_modeller", TCR,"_", x, "_percentile_cutoff.png", sep=""), 
      height=1800, width=1800, res=300)
  plot(perf, col=as.list(mycolors), main=TCR, lwd=2, xlim=c(0,1), ylim=c(0,1))
  #abline(a=0, b=1, col="black", lty=2)
  auc.perf <- performance(pred, measure = "aucpr")
  legend("topright",paste(colnames(scores), rep(": ", 2),
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
  mod_negatives <- read.table(paste(wkdir, "/AUC/modeller_score_table_preselection_remove_rnd4.txt",
                                sep=""), header=TRUE, sep="\t")
  mod_positives <- read.table(paste(wkdir, "/AUC/modeller_score_table_round4.txt",
                                sep=""), header=TRUE, sep="\t")
  crystal_negatives <- read.table(paste(wkdir, "/AUC/score_table_hamming_blosum_ADLIAYLKQATKG_preselection_remove_rnd4.txt",
                                sep=""), header=TRUE, sep="\t")
  crystal_positives <- read.table(paste(wkdir, "/AUC/score_table_hamming_blosum_ADLIAYLKQATKG_round4.txt",
                                sep=""), header=TRUE, sep="\t")
}

negatives <- cbind(mod_negatives[c(1,2,5)], crystal_negatives[5])
colnames(negatives) <- c("peptide", "reads", "model_score", "crystal_score")

positives <- cbind(mod_positives[c(1,2,5)], crystal_positives[5])
colnames(positives) <- c("peptide", "reads", "model_score", "crystal_score")

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
scores <- df[c("model_score", "crystal_score")]
all_scores <- data.frame(apply(scores, 2, function(x) x*(-1) ))

# Make labels dataframe with same dimensions as rank_df
df_lab <- data.frame(df$labels, df$labels)

# Plot ROC
plot_ROC(all_scores, df_lab, args$t, wkdir, args$x)
plot_PR(all_scores, df_lab, args$t, wkdir, args$x)



