library("ggplot2")


###############################################################
plotTheme <- theme_bw() +
             theme(
               plot.title = element_text(hjust = 0.5, size=20),
               axis.title = element_text(size=18),
               axis.text = element_text(size=12)
             )

# Get command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) <= 1) {
    stop("Usage: Rscript plotPredictedCutsInAlnForHMM.R <conservation_scores_file> <output_file>")
    }
conservation_scores_file = args[1]
if (length(args) > 1) {
  output_file = args[2]
} else {
  output_file = paste0("alnCuts.pdf")
}
if (length(args) > 2) {
  cutOff = as.double(args[3])
} else {
  cutOff = 0.2
}
###############################################################

n <- 5

# Conservation plot and finding Habc motif
scores <- as.numeric(unlist(strsplit(readLines(conservation_scores_file), " ")))
df <- data.frame(x=1:length(scores), y=scores)
df$z <- df$y
for (n in 1:n) {
  for (i in 1:length(scores)) {
    if (i >= length(scores)){
      z <- mean(c(df$z[i], df$z[i]))
    } else {
      z <- mean(c(df$z[i], df$z[i+1]))
    }
    df$z[i] <- z
  }
}

minimums <- c()
for (i in 2:(length(scores)-1)) {
  if (df$z[i-1] > df$z[i] & df$z[i] < df$z[i+1]) {
    minimums <- append(minimums, i)
  }
}
dfm <- df[minimums, ]
dfm <- dfm[dfm$z<=cutOff, ]

minimums <- c()
for (i in 2:(length(scores)-1)) {
  if ((df$z[i-1] < cutOff & df$z[i+1] > cutOff) | (df$z[i-1] > cutOff & df$z[i+1] < cutOff)) {
    minimums <- append(minimums, i)
  }
}
df02 <- df[minimums, ]
df02 <- df02[df02$z>cutOff, ]

pdf(output_file, width=20, height=5)
#ggplot(df, aes(x=df$x+n/2, y=df$z)) +
ggplot(df, aes(x=df$x, y=df$z)) +
  ylab("Conservation") +
  xlab("Alignment Position") +
  geom_line(size=1, color="red", legend = FALSE) +
  geom_line(aes(x=df$x, y=df$y), alpha=0.2, color="blue") +
  geom_hline(yintercept = cutOff, color="pink") +
  geom_vline(xintercept = dfm$x, color="orange") +
  geom_vline(xintercept = df02$x, color="darkgreen") +
  annotate(geom="text", x=df02$x, y=0, label=df02$x, hjust=-.5, vjust=-.5, angle=90, color="darkgreen") +
  annotate(geom="text", x=dfm$x, y=0, label=dfm$x, hjust=.5, vjust=-.5, angle=90, color="orange") +
  plotTheme
dev.off()





