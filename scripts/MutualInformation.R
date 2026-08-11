########## plotting AMI calculated with python ##########
#########################################################

########## Libraries ##########
library(ggplot2)
library(tidyr)
library(tibble)
library(ggdendro)
library(grid)
library(seqinr)
library(ggseqlogo)
library(dplyr)
library(cowplot)

########## Functions ##########
savePlot <- function(ggObject, plotname, path=MIpath, width=10, height=7){
  pdf(file = file.path(path, plotname), width = width, height = height)
  print(ggObject)
  dev.off()
}

########## Variables ##########
args = commandArgs(trailingOnly=TRUE)
AlignmentFile <- toString(args[1])
outDir <- toString(args[2])                 # e.g. coevolution/MI/
AMIcutOff <- toString(args[3])
MIpath <- file.path(outDir)
theme <- theme_bw() +
         theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
               axis.title = element_text(hjust = 0.5, size = 16, face = "bold"),
               axis.text = element_text(size = 16))

##########    Code   ##########
### Occupancy plot
occupancy <- read.csv(file.path(MIpath, "occupancy.txt"), header = FALSE)
colnames(occupancy) <- "Ocuppancy"
occupancy$id <- seq(1,length(occupancy[,1]))
plot_occupancy <- ggplot(occupancy, aes(x=id, y=Ocuppancy, fill=Ocuppancy)) +
                    geom_bar(stat="identity", show.legend = FALSE) +
                    scale_fill_viridis_b() +
                    ggtitle("Ocuppancy per residue in alignment") +
                    xlab("Residues") +
                    ylab("Ocuppancy") +
                    theme
savePlot(plot_occupancy, "Occupancy.pdf")


### Entropy plot
entropy <- read.csv(file.path(MIpath, "entropy.txt"), header = FALSE)
colnames(entropy) <- "Entropy"
entropy$id <- seq(1,length(entropy[,1]))
plot_entropy <- ggplot(entropy, aes(x=id, y=Entropy, fill=Entropy)) +
  geom_bar(stat="identity") +
  scale_fill_viridis_b() +
  ggtitle("Entropy per residue in alignment") +
  xlab("Residues") +
  ylab("Entropy") +
  theme
savePlot(plot_entropy, "Entropy.pdf")


### Conservation plot
conservation <- read.csv(file.path(MIpath, "conservation.txt"), skip = 1, sep="\t", header = T)
colnames(conservation) <- c("column_number", "Conservation", "column")
conservation$id <- seq(1,length(conservation[,1]))
plot_conservation <- ggplot(conservation, aes(x=id, y=Conservation, fill=Conservation)) +
  geom_bar(stat="identity") +
  scale_fill_viridis_b() +
  ggtitle("Conservation per residue in alignment") +
  xlab("Residues") +
  ylab("Conservation") +
  theme
savePlot(plot_conservation, "Conservation.pdf")


### AMI
AMI <- read.csv(file.path(MIpath, "AMI.txt"), header = FALSE, sep = ",")
#AMIjaccard <- read.csv(file.path(MIpath, "AMIjaccard.txt"), header = FALSE, sep = ",")

n_row <- length(AMI[,1])
n_col <- length(AMI[1,])

colnames(AMI) <- seq(1,n_col)
rownames(AMI) <- seq(1,n_row)

AMI_melt <- AMI %>%
  as.data.frame() %>%
  rownames_to_column("f_id") %>%
  pivot_longer(-c(f_id), names_to = "samples", values_to = "AMI")

# Sort residues by order in protein sequence
AMI_melt_res <- AMI_melt
AMI_melt_res$samples <- factor(AMI_melt$samples, levels=order(levels(factor(AMI_melt$samples)), method = "shell"))
AMI_melt_res$f_id    <- factor(AMI_melt$f_id, levels=order(levels(factor(AMI_melt$f_id)), method = "shell"))
AMI_plot <- ggplot(AMI_melt_res, aes(x=f_id, y=samples, fill=AMI)) +
              geom_tile() +
              scale_fill_viridis_c() +
              ggtitle("Adjusted Mutual information (AMI; Sequence residue order)") +
              xlab("Residues") +
              ylab("") +
              theme + theme(axis.text = element_text(size = 8))
savePlot(AMI_plot, "AMI_resOrder.pdf")

# Sort residues by dendrogram clustering - Build distance matrix between residues in AMI to obtain dendrogram
dendro <- as.dendrogram(hclust(d = dist(x = AMI)))
dendro_order <- order.dendrogram(dendro)

AMI_melt_dendro <- AMI_melt_res
AMI_melt_dendro$f_id    <- factor(x = AMI_melt_dendro$f_id, levels = levels(factor(AMI_melt_dendro$f_id))[rev(dendro_order)], ordered = TRUE)
AMI_melt_dendro$samples <- factor(x = AMI_melt_dendro$samples, levels = levels(factor(AMI_melt_dendro$samples))[dendro_order], ordered = TRUE)

# Ploting
dendro_plot  <- ggdendrogram(data = dendro, rotate = TRUE) + theme(axis.text.y = element_text(size = 6))
heatmap_plot <- ggplot(data = AMI_melt_dendro, aes(x = f_id, y = samples)) +
                  geom_tile(aes(fill = AMI)) +
                  #scale_fill_gradient2() +
                  scale_fill_viridis_c() +
                  scale_y_discrete(position = "right") +
                  ggtitle("Adjusted Mutual Information (AMI; Clustered order)") +
                  xlab("Residues") + theme +
                  theme(axis.title.y = element_text(size = 0),
                        axis.text.y = element_blank(),
                        axis.text.x = element_text(size = 6),
                        legend.position = "top")

# Generating PDF with plot
pdf(file = file.path(MIpath, "AMI_dendroOrder.pdf"), width = 10)
grid.newpage()
print(heatmap_plot,
      vp = viewport(x = 0.41, y = 0.5, width = 0.80, height = 1.0))
print(dendro_plot,
      vp = viewport(x = 0.90, y = 0.425, width = 0.2, height = 0.91))
dev.off()



### Weblogo
# Read aligned sequences
aln <- read.fasta(AlignmentFile, as.string = TRUE, forceDNAtolower = FALSE)
seqs <- sapply(aln,"[[",1)

# Read significant interactions
sigRes <- read.table(file.path(MIpath, "significant_interactions.txt"), colClasses = "numeric", sep = "\t", header = TRUE, row.names = NULL)[,c(1,2,3)]
colnames(sigRes) <- c('residue1', 'residue2', 'AMI')
sigRes$residue1 <- sigRes$residue1
sigRes$residue2 <- sigRes$residue2
sigRes <- sigRes[order(sigRes$residue2),]
rownames(sigRes) <- NULL
# Get only significant residue pairs with entropy != 0
sigRes <- sigRes[!(entropy[sigRes$residue1,]$Entropy == 0) & !(entropy[sigRes$residue2,]$Entropy == 0), ]

#hmmLength <- 82
#sigRes <- sigRes[((sigRes$residue1 <= hmmLength & sigRes$residue2 >= hmmLength) | (sigRes$residue2 <= hmmLength & sigRes$residue1 >= hmmLength)) & sigRes$AMI>0.8,]


# Weblogo plot
plotTheme <- theme(legend.position = 'right',
                   plot.margin = unit(c(0, 0, 0, 1), "cm"),
                   axis.text.x = element_text(size = 10, angle = 45),
                   axis.line.x.bottom = element_line(),
                   axis.line.y.left = element_line(),
                   axis.title = element_text(face = "bold"))

cs1 = make_col_scheme(chars = LETTERS, groups = c(rep("gr1",length(LETTERS))), 
                      cols = c(rep("white",length(LETTERS))))

p1 <- ggseqlogo(data = seqs,
                col_scheme = 'chemistry') +                    # chemistry, hydrophobicity, clustalx, taylor
                plotTheme +
                theme(axis.text.x = element_blank()) 

p2 <- ggplot(conservation, aes(id, Conservation, fill=Conservation)) +
            geom_bar(stat='identity') + 
            scale_fill_viridis_c() +
            #scale_x_continuous("Residues", breaks = entropy$id) + 
            theme_logo() + 
            plotTheme +
            theme(axis.text.x = element_blank()) +
            theme(axis.title.x = element_blank()) +
            theme(plot.margin = margin(0, 0, 1, 0, unit = "cm"))


p3 <- ggplot(entropy, aes(id, Entropy, fill=Entropy)) +
            geom_bar(stat='identity', show.legend = FALSE) + 
            scale_fill_viridis_c() +
            theme_logo() + 
            plotTheme + 
            theme(axis.text.x = element_blank()) +
            theme(axis.title.x = element_blank()) +
            theme(plot.margin = margin(0, 0, 1, 0, unit = "cm"))

p4 <- ggplot(occupancy, aes(x=id, y=Ocuppancy, fill=Ocuppancy)) +
            geom_bar(stat="identity", show.legend = FALSE) +
            scale_fill_viridis_b() +
            scale_x_continuous("Residues", breaks = occupancy$id) + 
            theme_logo() + 
            plotTheme
            
# Highlights and annotation
n <- 0
highlights <- list()
highlighted <- c()
annotations <- list()

p1High <- max(p1$layers[[1]]$data$y + 0.2)
for (i in 1:nrow(sigRes)){
  
  # Read significant residues
  r1 <- as.numeric(sigRes[i,1])
  r2 <- as.numeric(sigRes[i,2])
  v <- as.numeric(sigRes[i,3])
  
  col <- case_when(
          v <= 0.75 ~ "grey",
          v >= 0.9  ~ "red",
          v >  0.75 ~ "green"
  )
  
  # Add highlighting for residues
  for (r in c(r1, r2)){
    if (r %in% highlighted){
      next
    } else {
      highlights <- append(highlights,
                          annotate('rect', xmin = r-0.5, xmax = r+0.5, ymin = -0.05, ymax = p1High, size = 0.05, alpha = .1, col='grey', fill='yellow'))
      highlighted <- c(highlighted, r)
    }
  } 
  
  # Add segment for residues pairs
  a <- annotate('segment', x = r1, xend=r2, y=n, yend=n, size=0.5, col = col,
                arrow = arrow(angle = 40, ends='both', length = unit(0.2,"cm")), lineend = 'butt', linejoin = 'bevel')
  annotations <- append(annotations, a)
  n <- n+0.05

}

# Annotation legend
#annotations <- append(annotations,
#                      annotate('text', x=1, y=n-0.05, hjust = 0, label=paste('AMI:', AMIcutOff,'(grey) < 0.75 (green) < 0.9 (red)') ))

p0 <- ggseqlogo(data = seqs[1], col_scheme = cs1)  +
            labs(title = "Weblogo - Higlighted residues with significant AMI",
                 subtitle = paste('AMI:', AMIcutOff,'(grey) < 0.75 (green) < 0.9 (red)')) +
            scale_color_manual(values = "white") +
            scale_x_continuous("", breaks = sort(unique(c(as.numeric(sigRes$residue1), as.numeric(sigRes$residue2))))) +          
            theme_minimal() + 
            annotations +
            ylim(-0.05, 0.05 * (length(annotations) - 1 ) ) +
            ylab(label = "Significant\npairs") +
            theme(legend.title = element_text(color = "white"),
                  legend.text = element_text(color = "white"),
                  axis.line = element_blank(),
                  axis.text.y = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(face = "bold"),
                  panel.grid = element_blank(),
                  plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

p1 <- p1 + highlights
plot_Weblogo <- plot_grid(p0, p1, p2, p3, p4,  ncol = 1, align = 'v', axis = 'rl', rel_heights = c(.2, .4, .2, .1, .1))

w <- length(conservation$id)*2/10
if (w < 15){
  W <- 15
}else if (w > 200){
  W <- 200
}else{
  W <- w
}
savePlot(plot_Weblogo, "Weblogo.pdf", width = W, height = 15)
