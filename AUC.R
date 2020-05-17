#!/usr/bin/env Rscript
library(argparse)
library(ROCR)

parser <- ArgumentParser(description=paste("Produce ROC curves and calculate AUC for input TCR",
                                                      "should be run in TCR/ folder", sep=" "))
parser$add_argument("-t", help="TCR (ex. 1YMM)", type="character", dest="t", required=TRUE)

args <- parser$parse_args()
wkdir <- getwd()
print (wkdir)
stop()

plot_ROC <- function(scores, labels, TCR) {
  # Plot ROC curves for data
  # Args:
  #   scores: dataframe with scores and separate
  #           scoring methods for each column
  #   labels: dataframe with labels 
  #           (same dimensions as scores)
  #   TCR: character string of TCR
  mycolors <- c("orange", "magenta", "dodgerblue", "red")
  pred <- prediction(flip_scores, df_lab)
  perf <- performance(pred,"tpr","fpr")
  png(paste("ROC_AUC_", TCR, ".png", sep=""), height=1500, width=1500, res=300)
    plot(perf, col=as.list(mycolors), main=args$t)
    abline(a=0, b=1, col="black", lty=2)
    auc.perf <- performance(pred, measure = "auc")
    legend("bottomright",paste(colnames(flip_scores), rep(": ", 4),
                             round(as.numeric(auc.perf@y.values),2), sep=""), 
         col=mycolors, title="AUC", lty=1, lwd=2)
  
  
  
  
}

t <- c()
args <- list(t)
args$t <- "1YMM"
wkdir <-  "C:/Users/tyler/Dropbox (UMass Medical School)/Research/TCR"

if (args$t == "1YMM") {
  negatives <- read.table(paste(wkdir, "/", args$t, 
                                "/AUC/score_table_hamming_ENPVVHFFKNIVTP_preselection_remove_rnd4.txt",
                                sep=""), header=TRUE, sep="\t")
  positives <- read.table(paste(wkdir, "/", args$t,
                          "/hamming/score_table_hamming_ENPVVHFFKNIVTP_round4.txt",
                                sep=""), header=TRUE, sep="\t") 
}
# add labels
labels <- rep(FALSE, nrow(negatives))
negatives <- cbind(negatives, labels)
labels <- rep(TRUE, nrow(positives))
positives <- cbind(positives, labels)

# Filter positives
#positives <- positives[positives$reads > 50,]

# Merge into one table
df <- rbind(negatives, positives)

# Flip sign of scores for ROC predictions
scores <- df[c("standard", "bind", "bind_pep_alone", "hamming")]
flip_scores <- data.frame(apply(scores, 2, function(x) x*(-1) ))

df_pred_lab <- cbind(flip_scores, df["labels"])

df_lab <- data.frame(df$labels, df$labels, df$labels, df$labels)
mycolors <- c("orange", "magenta", "dodgerblue", "red")
pred <- prediction(flip_scores, df_lab)
perf <- performance(pred,"tpr","fpr")
plot(perf, col=as.list(mycolors), main=args$t)
abline(a=0, b=1, col="black", lty=2)
auc.perf <- performance(pred, measure = "auc")
legend("bottomright",paste(colnames(flip_scores), rep(": ", 4),
                           round(as.numeric(auc.perf@y.values),2), sep=""), 
                           col=mycolors, title="AUC", lty=1, lwd=2)
      


perf <- performance(pred, "prec", "rec")
plot(perf, col=c("forestgreen", "magenta", "dodgerblue", "red"))
auc.perf <- performance(pred, measure = "aucpr")
auc.perf@y.values[[1]]


#c("forestgreen", "magenta", "dodgerblue", "red")
pred <- prediction(df_pred_lab$hamming, df_pred_lab$labels)
perf <- performance(pred,"tpr","fpr")
plot(perf, colorize=TRUE)
abline(a=0, b=1)
auc.perf <- performance(pred, measure = "auc")
auc.perf@y.values[[1]]
perf <- performance(pred, "prec", "rec")
plot(perf)
auc.perf <- performance(pred, measure = "aucpr")
auc.perf@y.values[[1]]


## TESTING ############################################################
# my_labels <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
# my_scores <- c(0.9, 0.8, 0.7, -0.1, 0.2, 0.4) *10
# 
# df <- simple_roc(my_labels, my_scores)
# plot(df$FPR, df$TPR, type="o", col="blue")
# 
# library(ROCR)
# data(ROCR.simple)
# pred <- prediction(ROCR.simple$predictions, ROCR.simple$labels)
# perf <- performance(pred,"tpr","fpr")
# plot(perf)
# 
# pred <- prediction(my_scores, my_labels)
# perf <- performance(pred, "tpr", "fpr")
# plot(perf)


# simple_roc <- function(labels, scores){
#   labels <- labels[order(scores, decreasing=TRUE)]
#   df <- data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
#   return(df)
# }
# 
# check_df <- simple_roc(df_pred_lab$labels, df_pred_lab$hamming)
# plot(check_df$FPR, check_df$TPR, type="l")
# abline(a=0, b=1)
## TESTING ############################################################
