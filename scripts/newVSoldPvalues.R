# Run this app with the working directory set to the repo root
# (e.g. RStudio "Set As Working Directory").
library("shiny")
library(bslib)

library("ggplot2")
library("dplyr")
library("reshape2")
library("ggrepel")
library("patchwork")
library("ggpubr")

###### SHINY APP ######
results <- read.csv("tmp/results.csv", header = TRUE, stringsAsFactors = FALSE)
results_wide <- read.csv("tmp/results_wide.csv", header = TRUE, stringsAsFactors = FALSE)
plotTheme <- theme(plot.title = element_text(size = 16, face = "bold"),
                   plot.subtitle = element_text(size= 14, face = "italic"),
                   axis.title = element_text(size = 14),
                   axis.text = element_text(size = 12))
                   

ui <-  page_fluid(
  
          tags$head(tags$style(HTML("
            .bslib-card, .tab-content, .tab-pane, .card-body {
              overflow: visible !important;
            }
          "))),
  
          titlePanel("HMM p-values for SNARE motifs!"),
        
          sidebarLayout(
            
            sidebarPanel(
              
              card(
                card_header("Select Fasta file"),
                selectInput("fasta", label = "Fasta", choices = sort(unique(results$fasta)), selected = "Qa.I"),
              ),
              
              card(
                card_header("Select Parameters to compare"),
                selectInput("value", label = "Value", choices = c("iEval", "cEval", "Score", "Bias"), selected = "iEval"),
                checkboxInput("showNames", label = "Show names", value = TRUE),
              ),
              
              card(
                card_header("Select HMM models"),
                selectInput("newhmm", label = "New HMM", choices = sort(unique(results[results$fasta == results$motif & results$newiEval > -25,]$motif)), selected = "Qa.I"),
                selectInput("oldhmm", label = "Old HMM", choices = sort(unique(results[results$fasta == results$motif & results$oldiEval > -25,]$motif)), selected = "Qa.I"),
              ),
              
              # card(
              #   card_header("Colored by ..."),
              #   selectInput("colorBy", label = "Var used for color", choices = c("Bias", "Score", "iEval", "cEval", "Length"), selected = "Bias"),
              # ),
              
              width = 3, height = "85vh"#"20vw", height = "80vh"
            ),
            
            mainPanel(
              
              navset_card_underline(
                
                nav_panel("Eval",           plotOutput("plot",  height = "85vh"), height = "100%"),
                nav_panel("Score and bias", plotOutput("plot2", height = "85vh"), height = "100%")
                
              ),
              
              width = 9, height = "100%" #, width = "80vw", height = "80vh")
              
            )
          
          )
        
        )
        

## SERVER
server <- function(input, output, session) {
  
  pvalPlot <- function(fasta, hmm1, hmm2, value, showNames=TRUE) {#, colorBy) {
    
    v1 <- paste0("new", value, "_", hmm1)
    v2 <- paste0("old", value, "_", hmm2)
    # c1 <- paste0("new", colorBy, "_", hmm1)
    # c2 <- paste0("old", colorBy, "_", hmm2)
    l1 <- paste0("newLength_", hmm1)
    l2 <- paste0("oldLength_", hmm2)
    
    #### Plots from results_wide
    subset1 <- results_wide[results_wide$fasta == fasta, c('id', v1, v2, l1, l2) ]
    subset1$Ldiff <- (abs(subset1[, l1] - subset1[, l2]) + mean(c(subset1[, l1], subset1[, l2]))) - 53
    subset1 <- subset1 %>% arrange(desc(Ldiff))
    colnames(subset1) <- c("id", "newPval", "oldPval", "newLength", "oldLength", "Ldiff")
    maxX <- max(subset1[subset1$oldPval > subset1$newPval, "newPval"])
    subset1 <- mutate(subset1, 
                      Lshape = case_when(newLength == oldLength ~ "0",
                                         newLength > oldLength ~ "1",
                                         newLength < oldLength ~ "-1"
                                         )
    )
    
    ## Scatterplot
    p1 <- ggplot(subset1, aes(x = newPval, y = oldPval)) +
              geom_point(aes(color = Ldiff, shape=Lshape), size = 3, alpha = 0.5) +
              annotate("line", x = c(-25, max(subset1$newPval, subset1$oldPval) + 10), y = c(-25, max(subset1$newPval, subset1$oldPval) + 10), linetype = "dashed", color = "grey") +
              annotate("line", x = c(maxX, maxX), y = c(-25, max(subset1$newPval, subset1$oldPval) + 10), linetype = "dashed", color = "red") +
              annotate("text", x = maxX + 2, y = -20, label = paste0(round(maxX, 2)), color = "red", size=4, hjust=0) +
              ylim(c(-25, max(subset1$newPval, subset1$oldPval) + 10)) +
              annotate("text", x = c(-20, -20, -20), y = c(70, 64, 58), 
                       label = c(paste0("n = ", nrow(subset1)),
                                 paste0("n = ", nrow(subset1[subset1$newPval > subset1$oldPval, ]), "(new > old; ", round(nrow(subset1[subset1$newPval > subset1$oldPval, ]) / nrow(subset1) * 100, 2), "%)"),
                                 paste0("n = ", nrow(subset1[subset1$newPval < subset1$oldPval, ]), "(new < old; ", round(nrow(subset1[subset1$newPval < subset1$oldPval, ]) / nrow(subset1) * 100, 2), "%)")),
                       color = c("black", "#4DAF4A", "darkred"), size=6, hjust=0) +
              scale_color_viridis_c(name = paste0("Length diff (abs(L1-L2)+mean(L1,L2)-53)"), begin=0.30) +
              scale_shape_discrete(name = paste0("HMM match Length"), labels = c("new < old", "new = old", "new > old")) +
              labs(title = paste0("HMM Score distribution for ", fasta, " motifs"),
                   subtitle = paste0("score = ", v1, " vs ", v2),
                   x = "new Score",
                   y = "old Score") +
              theme_minimal() +
              plotTheme

    if (showNames) {
      p1 <- p1 + geom_text_repel(aes(label = ifelse(subset1$newPval < subset1$oldPval, subset1$id, "")),
                                 size = 5, show.legend = FALSE, nudge_x = 0.5, nudge_y = 0.5, segment.color = "grey")
    }

    ## Boxplot
    subset2 <- results_wide[results_wide$fasta == fasta, c('id', v1, v2) ]
    msubset <- melt(subset2)
    my_comparisons <- list( c(v1, v2) )
    p2 <- ggplot(msubset, aes(x = variable, y = value, fill = variable)) +
              geom_boxplot(outlier.shape = NA, alpha = 0.5) +
              ylim(c(-25, max(msubset$value) + 10)) +
              stat_compare_means(label = "p.signif", method = "t.test", comparison = my_comparisons) +
              labs(subtitle = paste0("P-value distribution (new vs old) for ", fasta, " motifs"),
                   x = "iEvalue Type",
                   y = "iEvalue (-log scale)") +
              theme_minimal() +
              # scale_fill_manual(values = c("#4DAF4A", "#377EB8"), name = "P-value Type") +
              scale_fill_viridis_d() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              plotTheme
    
    
    #### Boxplot for all results
    # 0. Melt results
    mresults <- melt(results, id.vars = c("id", "fasta", "motif"), measure.vars = c(paste0("new", value), paste0("old", value)), variable.name = "Type", value.name = "Score")
    
    # 1. Calculate group means
    group_means <- mresults %>%
      group_by(fasta, motif, Type) %>%
      summarise(mean_value = median(Score))
    
    # 2. Filter the data
    fasta_filter <- fasta
    low_means <- group_means %>%
      filter(fasta==fasta_filter, mean_value < -25) %>%
      pull(motif) %>% unique()
    
    filtered_data <- mresults %>%
      filter(fasta==fasta_filter, !motif %in% low_means)
    
    # 3. Create the boxplot
    p3 <- ggplot(filtered_data, aes(x = motif, y = Score, fill = Type)) +
              geom_boxplot(outlier.shape = NA, alpha = 0.5) +
              # scale_fill_manual(values = c("#4DAF4A", "#377EB8"), name = "Score") +
              scale_fill_viridis_d() +
              labs(title = paste0("Score distribution for ", fasta, " motifs - if both motifs in new and old"),
                   x = "Motif",
                   y = "Score",
                   color = "Type") +
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              plotTheme
    
    ## Combine plots
    return(p1 + p2 + p3 +
             plot_layout(guides = "collect",
                         design = "11112
                                   11112
                                   33333") & 
             theme(legend.position = "bottom")
    )
        
  }
  scoreBiasPlot <- function(fasta, hmm1, hmm2, value) {#, colorBy) {
    
    v1 <- paste0("new", value, "_", hmm1)
    v2 <- paste0("old", value, "_", hmm2)
    s1 <- paste0("newScore_", hmm1)
    s2 <- paste0("oldScore_", hmm2)
    # c1 <- paste0("new", colorBy, "_", hmm1)
    # c2 <- paste0("old", colorBy, "_", hmm2)
    b1 <- paste0("newBias_", hmm1)
    b2 <- paste0("oldBias_", hmm2)
    l1 <- paste0("newLength_", hmm1)
    l2 <- paste0("oldLength_", hmm2)
    
    #### Plots from results_wide
    subset3 <- results_wide[results_wide$fasta == fasta, c('id', v1, v2, s1, s2, b1, b2, l1, l2) ]
    subset3$valueDiff <- subset3[, v1] - subset3[, v2] 
    subset3$Ldiff <- subset3[, l1] - subset3[, l2]
    subset3$ScoreDiff <- subset3[, s1] - subset3[, s2]
    subset3$BiasDiff <- subset3[, b1] - subset3[, b2]
    
    subset3 <- subset3 %>% arrange(-desc(Ldiff))
    colnames(subset3) <- c("id", "newPval", "oldPval", "newScore", "oldScore", "newBias", "oldBias", "newLength", "oldLength",
                           "valueDiff", "Ldiff", "ScoreDiff", "BiasDiff")
    
    ## Scatterplot
    p4 <- ggplot(subset3, aes(x = valueDiff , y = ScoreDiff, color = Ldiff)) +
              geom_point(size = 3, alpha = 0.8) +
              annotate("line", x = c(-50, 50), y = c(-50, 50), linetype = "dashed", color = "grey") +
              scale_color_viridis_c() +
              xlim(c(-55, 55)) +
              ylim(c(-55, 55)) +
              labs(title = paste0("delta ", v1, " vs delta ", v2),
                   subtitle = paste0("Delta score vs Delta ", value, " (new vs old)"),
                   x = paste0("Delta ", value),
                   y = "Delta score") +
              theme_minimal() +
              plotTheme
    
    return(p4)
    
  }
  
  output$plot <- renderPlot(pvalPlot(input$fasta,
                                     input$newhmm,
                                     input$oldhmm,
                                     input$value,
                                     input$showNames
                                     # input$colorBy
                  ))
  
  output$plot2 <- renderPlot(scoreBiasPlot(input$fasta,
                                           input$newhmm,
                                           input$oldhmm,
                                           input$value
                                           # input$colorBy
                  ))
}
        
## RUN APP                           
shinyApp(ui, server)


