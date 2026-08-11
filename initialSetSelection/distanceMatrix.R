# Heatmap of the distancematrix and dendrogram of alignment

# Libraries
library(seqinr)
library(dplyr)
library(reshape2)
library(ggdendro)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library("optparse")

# Input options
option_list = list(
  make_option(c("-a", "--alignment"), type="character", default="empty", help="Alignment file", metavar="character"),
  make_option(c("-f", "--format"), type="character", default="fasta", help="Alignment format", metavar="character"),
  make_option(c("-m", "--matrix_method"), type="character", default="similarity", help="Method used to compute distance matrix", metavar="character"),
  make_option(c("-c", "--cluster_method"), type="character", default="average", help="Method used to compute clusters", metavar="character"),
  make_option(c("-t", "--taxonomy"), type="character", default = "none", help="TAB delimited file with matching sequences IDs from alignment and their corresponding taxonomy group name. e.g(Seq1  Tax1)", metavar="character"),
  make_option(c("-g", "--group"), type="logical", default = FALSE, help="TAB delimited file with matching sequences IDs from alignment and their corresponding taxonomy group name. e.g(Seq1  Tax1)"),
  make_option(c("-s", "--sequences_name_in_plot"), type="logical", default=FALSE, help="bool; If TRUE includes seqnames in plot"),
  make_option(c("-o", "--output"), type="character", default="heatmap.pdf", help="Output file", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# check if alignment is provided
if (!(opt$alignment == "empty")) {
  if (!file.exists(opt$alignment)) {stop("Alignment file must exist.")}
  alignment <- opt$date
} else {
  ## Set the variables here for local runs (run with working directory set to
  ## initialSetSelection/, e.g. via RStudio "Set As Working Directory")
  group <- "Qa.III"
  
  opt$alignment <- paste0("../alignments/", group, ".trimmed.aln")
  opt$taxonomy <- "./seqsClassification.Qa.III.txt"
  opt$group <- TRUE
  opt$sequences_name_in_plot <- TRUE
  
  stop("Alignment file must be provided and exist. See script usage (--help)")
}

# read alignment  
myseqs <- read.alignment(file = opt$alignment, format = opt$format)
myseqs$nam <- gsub("/.*", "", myseqs$nam)
dm <- as.matrix(dist.alignment(myseqs, matrix = opt$matrix_method))

# sort by clusters
clust <- hclust(as.dist(dm), method = opt$cluster_method)
dm <- dm[clust$order, clust$order]
dmm <- melt(dm)

# DF for taxonomy and groups
taxColor <- data.frame("seqs" = rownames(dm))

if (!(opt$taxonomy=="none") & file.exists(opt$taxonomy)) {
  taxColor <- left_join(x = taxColor, y = read.table(opt$taxonomy, header = F), by = c("seqs" = "V1"))
  colnames(taxColor) <- c("seqs", "tax")
  taxColor[which(is.na(taxColor$tax)),2] <- "Other"
  taxColor$seqs <- factor(taxColor$seqs, levels = taxColor$seqs)
}

if (opt$group) {
  # Add groups
  taxColor$group <- "Other"
  groups <- list.files("./iterations/")[grepl(paste0(group, "\\..*.fasta"), list.files("./iterations/"))]
  for (file_group in groups) {
    file_name <- paste0("./iterations/", file_group)
    group_name <- gsub(".fasta", "", file_group)
    group_lines <- unlist(lapply(readLines(file_name)[grepl(pattern = ">", readLines(file_name))], function(x) strsplit(x, ">")[[1]][2]))
    taxColor <- mutate(taxColor, group = case_when(seqs %in% group_lines ~ group_name,
                                                   TRUE ~ group))
  }
}

# PLOTS
## dendrograms
dhc <- as.dendrogram(clust)
data <- dendro_data(dhc, type = "rectangle")
gg_dendro_top <- ggplot(segment(data)) + 
                    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), linewidth=0.1) +
                    theme_void() +
                    scale_x_continuous(expand=c(0,0)) +
                    theme(plot.margin = margin(10, 5, 2.5, 12))

gg_dendro_right <- ggplot(segment(data)) + 
                      geom_segment(aes(x = xend, y = yend, xend = x, yend = y), linewidth=0.1) +
                      coord_flip(expand = FALSE) +
                      theme_void() +
                      theme(plot.margin = margin(5, 5, 10, 5))

## taxonomy/group heatmap
if (!(opt$taxonomy=="none") & file.exists(opt$taxonomy)) {
  gg_tax_right <- ggplot(taxColor, aes(x=0, y=seqs, fill=tax)) +
                      geom_tile() +
                      scale_fill_manual(values = c("Metazoa.a" = "green", "Metazoa.b" = "darkgreen", "Metazoa.c" = "red",
                                                   "Fungi.a" = "steelblue1", "Fungi.b" = "navy",
                                                   "Viridiplantae.a" = "yellow", "Viridiplantae.b" = "orange",
                                                   "Other" = "black"),
                                        name = "Taxonomy") +
                      #scale_fill_manual(values = c("Metazoa" = "darkgreen", "Fungi" = "blue", "Viridiplantae" = "yellow", "Other" = "black"),
                      #                  name = "Taxonomy") +
                      theme_void() +
                      scale_x_discrete(expand=c(0,0)) +
                      theme(legend.text=element_text(size=40),
                            legend.title=element_text(size=40))

  gg_tax_top <- gg_tax_right + coord_flip() + 
                      theme(legend.position = "none",
                            plot.margin = margin(2.5, 2.5, 2.5, 10))
  
  tax_legend <- get_legend(gg_tax_right)
  gg_tax_right <- gg_tax_right + theme(plot.margin = margin(5, 2.5, 10, 2.5),
                                       legend.position = "none")
                      
}

if (opt$group) {
  gg_group_right <- ggplot(taxColor, aes(x=0, y=seqs, fill=group)) +
                        geom_tile() +
                        #scale_fill_manual(values = c("Other" = "black", palette()),
                        #                  name = "Groups") +
                        theme_void() +
                        scale_x_discrete(expand=c(0,0)) +
                        theme(legend.text=element_text(size=40),
                              legend.title=element_text(size=40))
  
  gg_group_top <- gg_group_right + coord_flip() + 
                        theme(legend.position = "none",
                              plot.margin = margin(2.5, 2.5, 2.5, 10))
  
  group_legend <- get_legend(gg_group_right)
  gg_group_right <- gg_group_right + theme(plot.margin = margin(5, 2.5, 10, 2.5),
                                           legend.position = "none")
                        
}

if (opt$sequences_name_in_plot) {
  gg_seqname_right <- ggplot(taxColor, aes(x=0, y=seqs, label=seqs)) +
                          geom_text() +
                          theme_void() + 
                          theme(plot.margin = margin(5, 2.5, 10, 2.5))
}

## heatmap
gg_heatmap <- ggplot(dmm, aes(Var1, Var2, fill = value)) +
                  geom_raster() +
                  scale_fill_gradient2(low = "black", mid = "red", high = "white", guide = "none") +
                  theme_void() +
                  scale_x_discrete(expand=c(0,0)) +
                  theme(plot.margin = margin(5, 5, 10, 10))

# labels
gg_group <- ggplot() + annotate(geom="text", x=0.05, y=0, hjust=0, label= "Group", size=10) + xlim(0.05, 0.5) +
                       theme(plot.margin = margin(2.5, 0, 2.5, 0)) + theme_void()
gg_taxonomy <- ggplot() + annotate("text", x=0.05, y=0, hjust=0, label= "Taxonomy", size=10) +xlim(0.05, 0.5) +
                       theme(plot.margin = margin(2.5, 0, 2.5, 0)) + theme_void()

# Plot and save
if (!is.na(opt$taxonomy) & file.exists(opt$taxonomy) & opt$group & opt$sequences_name_in_plot) {
  top_row <- plot_grid(gg_dendro_top, tax_legend, group_legend, 
                       nrow = 1, rel_widths = c(7/10, 1.5/10, 1.5/10))
  cluster_rows <- plot_grid(gg_tax_top, gg_taxonomy,
                            gg_group_top, gg_group, 
                            nrow = 2, rel_widths = c(7/10, 3/10))
  bottom_row <- plot_grid(gg_heatmap, gg_group_right, gg_tax_right, gg_seqname_right, gg_dendro_right,
                          nrow = 1, rel_widths = c(70/100, 2/100, 2/100, 8/100, 18/100))
  
  p1 <- plot_grid(top_row, cluster_rows, bottom_row, nrow = 3, rel_heights = c(16/100, 4/100, 70/100))
  
} else {
  p1 <- plot_grid(gg_dendro_top, get_legend(gg_tax),
                      gg_heatmap, gg_dendro_right,
                      nrow = 2, rel_heights = c(1/4, 3/4), rel_widths = c(3/4, 1/4))
}

ggsave(opt$output, p1, width = 30, height = 30)

