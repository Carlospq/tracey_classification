library("ggplot2")
library("ggtree", quietly = TRUE)
library("dplyr", quietly = TRUE)
library("stringr", quietly = TRUE)
library("tidytree", quietly = TRUE)
library("gridExtra")

# Run this script with the working directory set to initialSetSelection/
# (e.g. RStudio "Set As Working Directory", or setwd("initialSetSelection")
# from the repo root) since all paths below are relative to it.

group = "SNAP.b"
treeGroup = "SNAP"

# Plot eval distribution of all sequence from treeGroup
for (i in 1){
  lines <- readLines(paste0('iterations/', group, '.results'))
  filteredLines <- lines[!grepl("#", lines)]
  df_distribution <- NULL;
  for (l in filteredLines){
    #Some code that generates new row
    df_distribution <- rbind(df_distribution, strsplit(l, "\\s+")[[1]])
  }
  df_distribution <- as.data.frame(df_distribution)
  df_distribution[,8] <- as.numeric(df_distribution[,8])
  df_distribution[,1] <- unlist(lapply(df_distribution[,1], function(x) unlist(strsplit(x, split='/'))[[1]]))
  #df_distribution <- df_distribution[!(df_distribution[dim(df_distribution)[2]] == "dupelabel1"), ]
  colnames(df_distribution)[1] <- "label"
  colnames(df_distribution)[3] <- "group"
  colnames(df_distribution)[4] <- "group2"
  colnames(df_distribution)[8] <- "evals"

  df_distribution <- df_distribution %>% 
            group_by(label) %>%
            slice_max(-log2(evals), n = 1) %>%
            ungroup()  
  
  p1 <- ggplot(df_distribution, aes(x=-log2(evals), fill=group)) +
            geom_density(alpha=.2) + 
            geom_vline(xintercept = mean(-log2(df_distribution$evals)), linetype="dashed", color="red", linewidth=1) +
            theme_minimal()
  
  d <- density(-log2(df_distribution$evals))
  p1color <- ggplot(data.frame(x = d$x, y = d$y), aes(x, y)) + geom_line(linewidth=2) + 
                geom_segment(aes(xend = x, yend = 0, colour = x), linewidth=2) + 
                scale_color_gradient2(name = "-log2(eval)", low = "red", mid = "yellow", high = "darkgreen", midpoint = mean(-log2(df_distribution$evals))) +
                geom_vline(xintercept = 70, linetype="dotted", color="red", linewidth=1)
}
p1color


# Tree colored by evalues according to HMM of "group"
for (i in 1) {
  
  plotTitle <- group
  plotFileName <- paste0("../trees/", treeGroup, ".trimmed.aln.treefile")
  tree <- read.tree(file = plotFileName)
  tree$tip.label <- unlist(lapply(strsplit(tree$tip.label, split = "/"), "[", 1))
  
  treeDfFilled <- left_join(as_tibble(tree), df_distribution[,c(1,3,8)], by = "label", keep = FALSE) %>% 
                        mutate(evals = ifelse(is.na(evals), 0, evals),
                               group = ifelse(is.na(group), "Other", group))
  
  midpoint = mean(-log2(df_distribution$evals))
  treeTax <- as.treedata( treeDfFilled )

  g1 <- ggtree(treeTax, color="lightgrey", layout="equal_angle") +
            geom_tiplab(aes(color=-log2(evals), label=label), size=2) +
            scale_color_gradient2(name = "-log2(eval)", low = "red", mid = "yellow", high = "darkgreen", midpoint = midpoint, guide = "none") +
            ggtitle(str_to_title(plotTitle)) +
            theme(plot.title = element_text(hjust = 0.5, size = 24),
                  legend.position = "top",
                  legend.title = element_text(size = 18),
                  legend.text = element_text(size = 12),
                  legend.background = element_rect(fill="lightgrey", linetype="solid"))
}
grid.arrange(p1color + theme(plot.margin = unit(c(1,0,0,0), "cm")),
             g1,
             ncol = 1, nrow = 2,
             heights = c(1, 5))


# Tree colored for selected sequences of "group"
for (i in 1) {
  
  selectedSeqs <- read.csv(paste0("iterations/", group, ".fasta"), sep=' ', header=FALSE)[,1]
  selectedSeqs <- data.frame(shortname = unlist(lapply(grep("^>", selectedSeqs, value=TRUE), function(x) unlist(strsplit(x, split='>'))[2])), selected = "Selected")
  
  plotTitle <- group
  plotFileName <- paste0("../trees/", treeGroup, ".trimmed.aln.treefile")
  tree <- read.tree(file = plotFileName)
  tree$tip.label <- unlist(lapply(strsplit(tree$tip.label, split = "/"), "[", 1))
  
  treeDf <- left_join(as_tibble(tree), selectedSeqs, by = c("label" = "shortname"))
  treeDfFilled <- left_join(treeDf, df_distribution[,c(1,3,8)], by = "label", keep = FALSE) %>% 
    mutate(selected = ifelse(is.na(selected), "Not selected", "Selected"),
           evals = ifelse(is.na(evals), 1, evals))
    
  treeTax <- as.treedata( treeDfFilled )
  g2 <- ggtree(treeTax, color="lightgrey", layout="equal_angle") +
          geom_tiplab(aes(color=selected, label=label), size=2) +
          scale_color_manual(values = c("Selected" = "#35B779FF", "Not selected" = "#31688EFF"), guide="none") +
          #geom_tiplab(aes(color=-log2(evals), label=label), size=2) +
          #scale_color_viridis_c() +#(low="black", high="red") +
          ggtitle(str_to_title(plotTitle)) +
          theme(plot.title = element_text(hjust = 0.5, size = 24),
legend.position = "top",
legend.title = element_text(size = 18),
legend.text = element_text(size = 12),
legend.background = element_rect(fill="lightgrey", linetype="solid"))
}
g2


# Tree colored by selected groups
for (i in 1) {
  
  groupsList <- list.files(path = 'initialSets/')
  groups <- groupsList[grep(paste0(treeGroup, "\\..*\\.initialSet\\.txt"), groupsList)]
  
  groups <- unique(unlist(lapply(groups, function(x) gsub(paste0(treeGroup,'\\.(.*)\\.initialSet\\.txt'), '\\1', x))))
  
  #groups <- unique(unlist(lapply(groups, function(x) unlist(strsplit(x, split='\\.'))[[3]])))
  #groups <- c("Fungi", "Metazoa", "SAR", "Viridiplantae")

  initial_df <- data.frame()
  for (group in groups) {
    selectedSeqs <- read.csv(paste0("initialSets/", treeGroup, ".", group, ".initialSet.txt"), sep=' ', header=FALSE)[,1]
    selectedSeqs <- data.frame(shortname = selectedSeqs, selected = paste0("Initial_", group))
    initial_df <- rbind(initial_df, selectedSeqs)
  }
  
  selected_df <- data.frame()
  for (group in groups) {
    selectedSeqs <- read.csv(paste0("iterations/", treeGroup, ".", group, ".fasta"), sep=' ', header=FALSE)[,1]
    selectedSeqs <- data.frame(shortname = unlist(lapply(grep("^>", selectedSeqs, value=TRUE), function(x) unlist(strsplit(x, split='>'))[2])), selected = group)
    selectedSeqs <- selectedSeqs[!(selectedSeqs$shortname %in% initial_df$shortname), ]
    selected_df <- rbind(selected_df, selectedSeqs)
  }

  df <- rbind(initial_df, selected_df)
  n_occur <- data.frame(table(df$shortname))
  n_occur_names <- n_occur[n_occur$Freq > 1, 1]
  
  df <- df %>% mutate(selected = ifelse(shortname %in% n_occur_names, "Both", selected))
  df <- df[!duplicated(df), ]
  
  plotTitle <- group
  plotFileName <- paste0("../trees/", treeGroup, ".trimmed.aln.treefile")
  tree <- read.tree(file = plotFileName)
  tree$tip.label <- unlist(lapply(strsplit(tree$tip.label, split = "/"), "[", 1))
  
  treeDf <- left_join(as_tibble(tree), df, by = c("label" = "shortname")) %>% 
        mutate(selected = ifelse(is.na(selected), "Not selected", selected))
  
  treeTax <- as.treedata( treeDf )
  g3 <- ggtree(treeTax, color="lightgrey", layout="equal_angle", mapping = aes(color=selected)) +
    geom_tiplab(aes(color=selected, label=label), size=2, show.legend = FALSE) +
    #scale_color_viridis_d() +
    scale_color_discrete() +
    geom_point(size = 0, shape = 15) +
    ggtitle(str_to_title(plotTitle)) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme(plot.title = element_text(hjust = 0.5, size = 24),
          legend.position = "top",
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 12),
          legend.background = element_rect(fill="lightgrey", linetype="solid"))
  
}
g3
# 


##### Code for df duplicated values in column #####
n_occur <- data.frame(table(df_distribution$label))
df_distribution[df_distribution$label %in% n_occur$Var1[n_occur$Freq > 1],]

n_occur <- data.frame(table(as_tibble(tree)$label))
as_tibble(tree)[as_tibble(tree)$label %in% n_occur$Var1[n_occur$Freq > 1],]
###################################################
  
  






