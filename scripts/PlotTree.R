library("ape", quietly = TRUE)
library("ggplot2", quietly = TRUE)
library("ggtree", quietly = TRUE)
library("tidytree", quietly = TRUE)
library("stringr", quietly = TRUE)

library("seqinr", quietly = TRUE)
library("dplyr", quietly = TRUE)

# Run with the working directory set to the repo root.
#args = commandArgs(trailingOnly=TRUE)

# Read classification file
group = "kingdom"
shortnameToId <- read.csv("utils/sequenceIDtaxonomyID.csv", header = T)
idToTaxonomy  <- read.csv("utils/taxonomyTable.csv", header = T)
colnames(idToTaxonomy)[1] <- "tracey_taxonomy_id"
df_taxonomies <- merge(shortnameToId, idToTaxonomy, by = "tracey_taxonomy_id", all.x = T)
df_taxonomies <- df_taxonomies[, c('shortname', group)]

treefiles <- dir('trees')
for ( nwk in treefiles[endsWith(treefiles, suffix = "treefile")] ) {
  
  # Input files
  #nwk <- "Ha.II.trimmed.aln.treefile"
  
  plotTitle <- paste0(c(head(strsplit(nwk, "\\.")[[1]], -3)), collapse = ".")
  plotFileName <- paste0("trees/", plotTitle, ".pdf")
  # Skip if pdf already exists
  if (file.exists(plotFileName)) {next}
  print(nwk)

  
  
  # Read tree and Rename tip.labels
  tree <- read.tree(file = paste0("trees/", nwk))
  #split_shortnames <- unlist(strsplit(tree$tip.label, split = "\\|"))
  split_shortnames <- unlist(lapply(strsplit(tree$tip.label, split = "/"), "[", 1))
  tree$tip.label <- unlist(lapply(strsplit(tree$tip.label, split = "/"), "[", 1))
  #if (length(split_shortnames) == length(tree$tip.label)) {
  #  tree$tip.label <- unlist(strsplit(tree$tip.label, split = "\\|"))
  #} else {
  #  tree$tip.label <- unlist(strsplit(tree$tip.label, split = "\\|"))[2*(1:length(tree$tip.label))-1]
  #}
  
  # Get only sequences in tree from classification
  matches <- match(tree$tip.label, df_taxonomies$shortname)
  df_tax <- df_taxonomies[matches,]
  df_tax$group <- paste0(df_tax$shortname, "_", df_tax[,group])
  df_tax <- df_tax[!duplicated(df_tax), ]
  colnames(df_tax) <- c("shortname", "taxonomy", "shortnameTax")
  tree$tip.label <- df_tax$shortname
  
  # Split tree by groups
  treeTax <- as.treedata( full_join(as_tibble(tree), df_tax, by = c("label" = "shortname")) )
  tree_plot <- ggtree(treeTax, color="lightgrey", layout="equal_angle") +
                      geom_tiplab(aes(color=taxonomy, label=label), size=2) +
                      scale_color_manual(values = c("Metazoa" = "#35B779FF", "Viridiplantae" = "#FDE725FF", "Fungi"= "#31688EFF", "-" = "lightgrey"), guide="none") +
                      geom_polygon(aes(fill=taxonomy, x = 0, y = 0)) +
                      scale_fill_discrete(type = c("Metazoa" = "#35B779FF", "Viridiplantae" = "#FDE725FF", "Fungi"= "#31688EFF", "Other" = "lightgrey"),
                                          limits = c("Metazoa", "Viridiplantae", "Fungi", "Other")) +
                      ggtitle(str_to_title(plotTitle)) +
                      theme(plot.title = element_text(hjust = 0.5, size = 24),
                            legend.position = "top",
                            legend.title = element_text(size = 18),
                            legend.text = element_text(size = 12),
                            legend.background = element_rect(fill="lightgrey", linetype="solid"))
  #tree_plot
  ggsave(plotFileName, tree_plot, width = 30, height = 30)
  #pdf(file=plotFileName, width = 30, height = 30)
  #tree_plot
  #dev.off()

}
