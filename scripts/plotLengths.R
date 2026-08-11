# Plot distribution of length for sequences in 'fasta_file/args[1]'
library("ggplot2", quietly = TRUE)
###############################################################
plotTheme <- theme_bw() +
             theme(
               plot.title = element_text(hjust = 0.5, size=20),
               axis.title = element_text(size=18),
               axis.text = element_text(size=12)
             )

# Get command line arguments
args = commandArgs(trailingOnly=TRUE)
fasta_file = args[1]
if (length(args) > 1) {
  output_file = args[2]
} else {
  output_file = paste0("sequenceLengths.pdf")
}
###############################################################
file_name <- fasta_file
lengths <- c()
con = file(file_name, "r")
for (line in readLines(con)) {
  if (startsWith(line, ">")) {
    next
  }
  lengths <- append(lengths, nchar(line))
}

# Histogram and density plots
binwidth = 15
L <- data.frame(lengths=lengths)
ggh <- ggplot(L, aes(x=lengths)) +
            geom_histogram(binwidth = binwidth, color = "black", fill = "#00AFBB") +
            geom_density(aes(y = after_stat(density) * (nrow(L) * binwidth), alpha=.2), color="red", fill="pink", show.legend = F) +
            scale_y_continuous(name = "Counts") +
            ggtitle("Sequences Length") +
            xlab("Length") +
            plotTheme

# Get highest bin value
ggh_build <- ggplot_build(ggh)
counts_df <- ggh_build[["data"]][[1]]
top <- counts_df[counts_df$y==max(counts_df$y),]

# Plot
pdf(output_file, width=10, height=5)
ggh + annotate(geom = "text", x=top$x, y=top$y, label=paste0("Length: ",top$x), hjust = -0.2, color = "red")
dev.off()
