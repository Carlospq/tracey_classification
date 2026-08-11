library(ggplot2)
# Run this script with the working directory set to the repo root
# (e.g. RStudio "Set As Working Directory").
mainGroup = "QIII"
#output_file = "testEvaluesDistribution.pdf"
###############################################################

groups <- list("QI"   = c("Qa.I", "Qb.I", "Qc.I", "R.I", "Ha.I", "Hb.I", "Hc.I",
                          "Qa.I.Syx18.I", "Ha.I.Syx18.I", "Qa.I.Ufe1.I", "Ha.I.Ufe1.I"),
               "QII"  = c("Qa.II", "Qb.II", "Qc.II", "R.II", "Ha.II", "Hb.II",
                          "Qa.II.Sed5.II", "Qa.II.Syx5.II", "Qb.II.Gos1.II", "Qb.II.Bos1.II", "Qb.II.Membrin.II",
                          "Hb.II.Gos1.II", "Hb.II.Bos1.II", "Hb.II.Membrin.II"),
               "QIII" = c("Qa.III", "Qb.III", "Qc.III", "R.III", "Ha.III", "Hb.III", "Hc.III",
                          "Ha.III.Ha.III.b", "Ha.III.Ha.III.a",    "Qa.III.Qa.III.b", "Qa.III.Qa.III.a",
                          "Hb.III.Hb.III.b", "Hb.III.Hb.III.d",    "Qb.III.Qb.III.b", "Qb.III.Qb.III.d",
                          "Hc.III.Hc.III.b", "Hc.III.Hc.III.c",    "Qc.III.Qc.III.b", "Qc.III.Qc.III.c"),
               "QIV"  = c("Ha.IV", "Qa.IV",  "R.IV", "R.Reg",
                          "Ha.IV.Sso.IV",    "Ha.IV.Syx.IV",       "Qa.IV.Sso.IV",    "Qa.IV.Syx.IV"),
               "SNAP" = c("SNAP", "SNAPb", "SNAPc")
)

df_colors <- data.frame(motif = c("Ha.I",   "Hb.I",   "Hc.I",   "Qa.I",   "Qb.I",   "Qc.I",   "R.I",
                                  "Ha.II",  "Hb.II",  "Hc.II",  "Qa.II",  "Qb.II",  "Qc.II",  "R.II",
                                  "Ha.III", "Hb.III", "Hc.III", "Qa.III", "Qb.III", "Qc.III", "R.III",
                                  "Ha.IV",  "Qa.IV",  "R.IV",   "R.Reg",
                                  "Ha.I.Syx18.I",    "Ha.I.Ufe1.I",      "Qa.I.Syx18.I", "Qa.I.Ufe1.I",
                                  "Qa.II.Syx5.II",   "Qa.II.Sed5.II",
                                  "Hb.II.Gos1.II",   "Hb.II.Membrin.II", "Hb.II.Bos1.II",
                                  "Qb.II.Gos1.II",   "Qb.II.Membrin.II", "Qb.II.Bos1.II",
                                  "Ha.III.Ha.III.b", "Ha.III.Ha.III.a",  "Qa.III.Qa.III.b", "Qa.III.Qa.III.a",
                                  "Hb.III.Hb.III.b", "Hb.III.Hb.III.d",  "Qb.III.Qb.III.b", "Qb.III.Qb.III.d",
                                  "Hc.III.Hc.III.b", "Hc.III.Hc.III.c",  "Qc.III.Qc.III.b", "Qc.III.Qc.III.c",
                                  "Ha.IV.Sso.IV",    "Ha.IV.Syx.IV",     "Qa.IV.Sso.IV",    "Qa.IV.Syx.IV",
                                  "SNAP",            "SNAPb",            "SNAPc"),
                        color = c("peachpuff1", "yellow3",    "#D0BFD4", "peachpuff1", "yellow3", "#D0BFD4","#D89B64",
                                  "peachpuff1", "yellow3",    "#D0BFD4", "peachpuff1", "yellow3", "#D0BFD4","#D89B64",
                                  "peachpuff1", "yellow3",    "#D0BFD4", "peachpuff1", "yellow3", "#D0BFD4","#D89B64",
                                  "peachpuff1", "peachpuff1", "#D89B64", "#D89B64",
                                  "#FEABD5",    "#F9B29E",    "#FEABD5", "#F9B29E",
                                  "#FEABD5",    "#F9B29E",
                                  "#E1ECA6",    "#C9E5CF",    "#BED5A4",
                                  "#E1ECA6",    "#C9E5CF",    "#BED5A4",
                                  "#FEABD5",    "#F9B29E",    "#FEABD5", "#F9B29E",
                                  "#E1ECA6",    "#C9E5CF",    "#E1ECA6", "#C9E5CF",
                                  "#A6BAD8",    "#BDB2D0",    "#A6BAD8", "#BDB2D0",
                                  "#FEABD5",    "#F9B29E",    "#FEABD5", "#F9B29E",
                                  "#A6BAD8", "#BDB2D0","#A6BAD8"
                                  )
                        )

# Get command line arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plotEvaluesDistributionHMMs.R <main_Group[QI, QII, QIII, QIV]> <output_file>")
}
mainGroup = args[1]
output_file = args[2]
if (length(args) <= 3) {
  ylim <- args[3]
}
if (length(args) <= 4) {
  width <- args[4]
}
if (length(args) <= 5) {
  height <- args[5]
}

plotTheme <- theme_bw() + theme(plot.title = element_text(hjust = 0.5, size=20),
                                axis.title = element_text(size=18),
                                axis.text = element_text(size=12),
                                plot.subtitle = element_text(hjust=0.5)
                                )

distributionPlot <- function(df, output_file, plotTitle, ylim=0.5, width=width, height=height) {
  pdf(file=output_file, width=width, height=height)
  ggplot(dfValues, aes(x=log(eval), color=group, fill=group)) +
      ggtitle(plotTitle) +
      labs(subtitle = "Rows=Sequences groups; Columns=HMM models") +
      ylab("Density") + 
      xlab("Log(e-values)") +
      geom_density() +
      ylim(0, 0.05) +
      geom_vline(xintercept = log(1e-10), color="red") +
      annotate(geom="text", x=log(1e-10), y=0, vjust=1.5, hjust=1.1, label="1e-10", color="red") +
      facet_grid(group ~ motif , scales = "free_y") +
      scale_color_manual(values = alpha(unique(df_colors$color), 1)) +
      scale_fill_manual(values  = alpha(unique(df_colors$color), 0.7)) +
      plotTheme
  dev.off()
}
###############################################################

# Read taxonomy data
shortnameToId <- read.csv("utils/sequenceIDtaxonomyID.csv", header = T)
idToTaxonomy  <- read.csv("utils/taxonomyTable.csv", header = T)
colnames(idToTaxonomy)[1] <- "tracey_taxonomy_id"
df_taxonomies <- merge(shortnameToId, idToTaxonomy, by = "tracey_taxonomy_id", all.x = T)
df_taxonomies <- df_taxonomies[, c('shortname', 'kingdom')]

# Read evalues for each group
dfValues <- data.frame()
for (group in df_colors$motif) {
#for (group in groups[mainGroup][[1]]) {
  valuesFile <- paste0('regressionModel/', group, '.domtblout')
  if (!file.exists(valuesFile)) {next}
  print(group)
  dfTmp <- read.csv(valuesFile, header = F, sep = "")
  dfTmp$group <- group
  dfTmp <- dfTmp[,c("V1", "group", "V4", "V12")]
  colnames(dfTmp) <- c("shortname", "group", "motif", "eval")
  dfTmp$shortname <- sapply(dfTmp$shortname, function(x){ strsplit(x, "\\|")[[1]][1] })
  dfValues <- rbind(dfValues, dfTmp)
}
dfValues <- dfValues[!(duplicated(dfValues)),]
colnames(dfValues) <- c("shortname", "group", "motif", "eval")
dfValues$eval <- as.numeric(dfValues$eval)

# Merge evalues df with taxonomy and colors
dfValuesMerged <- merge(dfValues, df_taxonomies, by = "shortname")
dfValuesMerged <- merge(dfValuesMerged, df_colors, by= "motif")
dfValuesMerged <- dfValuesMerged[dfValuesMerged$motif %in% df_colors$motif,]
dfValuesMerged <- dfValuesMerged[dfValuesMerged$kingdom!="Bamfordvirae",]

# Plot distributions
distributionPlot(dfValuesMerged, output_file, "e-values distribution of Habc HMMs")


keepGroups <- c("Hb.III.Hb.III.b", "Hb.III.Hb.III.d")
df_ <- dfValuesMerged[dfValuesMerged$group %in% keepGroups & dfValuesMerged$motif %in% keepGroups,]

#df_$prot <- "none"
#df_[grepl("Syx16", df_$shortname), "prot"] <- "Syx16"
#pdf('regressionModel/QI.evalsDistribution.groupsColor.pdf', width=15, height=10)
ggplot(df_, aes(x=log(eval), color=kingdom, fill=kingdom)) +
  ggtitle("E-values distribution") +
  labs(subtitle = "Rows=Sequences groups; Columns=HMM models") +
  ylab("Density") + 
  xlab("Log(e-values)") +
  geom_density(linewidth=0.5) +
  #ylim(0, 0.2) +
  geom_vline(xintercept = log(1e-10), color="red") +
  annotate(geom="text", x=log(1e-10), y=0, vjust=1.5, hjust=1.1, label="1e-10", color="red") +
  facet_grid(group ~ motif , scales = "free_y") +
  scale_color_manual(guide="none", values = alpha(unique(df_colors$color), 1)) +
  scale_fill_manual(name="Group", values  = alpha(unique(df_colors$color), 0.7)) +
  plotTheme
#dev.off()



z <- df_[df_$group=="Qa.III" & df_$motif!="SNAPb",]
x <- z[z$motif=="Qa.III",]
y <- z[z$motif=="Qa.IV",]
z2 <- merge(x, y, by = "shortname")

ggplot(z2, aes(x=log(eval.x), y=log(eval.y), color=kingdom.x)) + geom_point(size=3) + geom_abline() + xlim(-150, 10) + ylim(-150, 10)


