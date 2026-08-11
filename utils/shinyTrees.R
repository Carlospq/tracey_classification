# ===== Section 1: Libraries =====
library(shiny)
library(ggplot2)
library(ggtree)
library(dplyr)
library(stringr)
library(tidytree)
library(DT)

# ===== Section 2: Global setup =====
# Run this app with the working directory set to the repo root
# (e.g. RStudio "Set As Working Directory").
ISEL_DIR  <- "initialSetSelection"
TREES_DIR <- "trees"

ALL_TREE_GROUPS <- sort(gsub("\\.trimmed\\.aln\\.treefile$", "",
                             list.files(TREES_DIR, pattern = "\\.trimmed\\.aln\\.treefile$")))
ALL_TREE_GROUPS <- ALL_TREE_GROUPS[!startsWith(ALL_TREE_GROUPS, "SNARE")]
ALL_TREE_GROUPS <- ALL_TREE_GROUPS[!startsWith(ALL_TREE_GROUPS, "H")]

# ===== Section 3: Helper functions =====

read_results <- function(group) {
  path <- file.path(ISEL_DIR, "iterations", paste0(group, ".results"))
  lines <- readLines(path)
  lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
  df <- do.call(rbind, lapply(lines, function(l) strsplit(trimws(l), "\\s+")[[1]]))
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df$label <- unlist(lapply(df[, 1], function(x) strsplit(x, "/")[[1]][1]))
  df$evals  <- as.numeric(df[, 8])
  # Replace zero e-values to avoid Inf in -log2()
  df$evals[df$evals == 0] <- 1e-300
  df %>%
    group_by(label) %>%
    slice_max(-log2(evals), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(label, evals)
}

read_fasta_ids <- function(filepath) {
  lines   <- readLines(filepath)
  headers <- grep("^>", lines, value = TRUE)
  ids     <- sub("^>", "", headers)
  unlist(lapply(strsplit(ids, "/"), `[`, 1))
}

read_initial_set <- function(filepath) {
  lines <- readLines(filepath)
  lines[nchar(trimws(lines)) > 0]
}

get_groups_for_tree <- function(treeGroup, all_groups) {
  all_groups[all_groups == treeGroup |
               startsWith(all_groups, paste0(treeGroup, "."))]
}

get_initial_set_labels <- function(treeGroup) {
  escaped <- gsub("\\.", "\\\\.", treeGroup)
  files <- list.files(file.path(ISEL_DIR, "initialSets"),
                      pattern = paste0("^", escaped, "(\\..*)?", "\\.initialSet\\.txt$"))
  if (length(files) == 0) return(character(0))
  inner_labels <- gsub(paste0("^", escaped, "(\\.(.+))?\\.initialSet\\.txt$"),
                       "\\2", files)
  inner_labels[inner_labels == ""] <- treeGroup
  setNames(file.path(ISEL_DIR, "initialSets", files), inner_labels)
}

# ===== Section 4: UI =====

ui <- fluidPage(
  titlePanel("Phylogenetic Tree Explorer — SNARE Classification"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      selectInput("treeGroup",
                  label = "Tree Group",
                  choices = ALL_TREE_GROUPS,
                  selected = ALL_TREE_GROUPS[1]),

      hr(),
      h5("Overlay: Iteration Results"),
      uiOutput("groupUI"),

      hr(),
      h5("Overlay: Initial Sets"),
      uiOutput("initialSetSubgroupsUI"),

      hr(),
      h5("Display Options"),
      sliderInput("tipLabelSize", "Tip label size",
                  min = 0.5, max = 5, value = 2, step = 0.25),
      radioButtons("treeLayout", "Tree layout",
                   choices = c("Equal angle" = "equal_angle",
                               "Fan"         = "fan",
                               "Rectangular" = "rectangular"),
                   selected = "equal_angle"),

      hr(),
      actionButton("renderBtn", "Render", class = "btn-primary btn-block")
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        id = "tabs",
        tabPanel("E-value Distribution",
                 br(),
                 plotOutput("plotEvalDensity", height = "350px"),
                 br(),
                 fluidRow(
                   column(3, numericInput("pdfWidthEvalDensity",  "Width (in)",  value = 8,  min = 1)),
                   column(3, numericInput("pdfHeightEvalDensity", "Height (in)", value = 5,  min = 1)),
                   column(3, br(), downloadButton("dlEvalDensity", "Download PDF"))
                 )),
        tabPanel("Tree: E-values",
                 br(),
                 plotOutput("plotTreeEval", height = "750px"),
                 br(),
                 fluidRow(
                   column(3, numericInput("pdfWidthTreeEval",  "Width (in)",  value = 12, min = 1)),
                   column(3, numericInput("pdfHeightTreeEval", "Height (in)", value = 12, min = 1)),
                   column(3, br(), downloadButton("dlTreeEval", "Download PDF"))
                 )),
        tabPanel("Tree: Selected Sequences",
                 br(),
                 plotOutput("plotTreeSelected", height = "750px"),
                 br(),
                 fluidRow(
                   column(3, numericInput("pdfWidthTreeSelected",  "Width (in)",  value = 12, min = 1)),
                   column(3, numericInput("pdfHeightTreeSelected", "Height (in)", value = 12, min = 1)),
                   column(3, br(), downloadButton("dlTreeSelected", "Download PDF"))
                 )),
        tabPanel("Tree: Initial Sets",
                 br(),
                 plotOutput("plotTreeInitialSets", height = "750px"),
                 br(),
                 fluidRow(
                   column(3, numericInput("pdfWidthTreeInitialSets",  "Width (in)",  value = 12, min = 1)),
                   column(3, numericInput("pdfHeightTreeInitialSets", "Height (in)", value = 12, min = 1)),
                   column(3, br(), downloadButton("dlTreeInitialSets", "Download PDF"))
                 )),
        tabPanel("Data Table",
                 br(),
                 DT::dataTableOutput("tableResults"))
      )
    )
  )
)

# ===== Section 5: Server =====

server <- function(input, output, session) {

  # --- Data discovery: groups available for the selected treeGroup ---
  available_groups <- reactive({
    req(input$treeGroup)
    all_results <- gsub("\\.results$", "",
                        list.files(file.path(ISEL_DIR, "iterations"), pattern = "\\.results$"))
    grps <- get_groups_for_tree(input$treeGroup, all_results)
    sort(grps)
  })

  # --- Data discovery: initial set subgroup files ---
  initial_set_files <- reactive({
    req(input$treeGroup)
    get_initial_set_labels(input$treeGroup)
  })

  # --- renderUI: group dropdown (rebuilt whenever treeGroup changes) ---
  output$groupUI <- renderUI({
    grps <- available_groups()
    if (length(grps) == 0) {
      selectInput("group", "Group (HMM iteration)",
                  choices = c("(none available)" = ""))
    } else {
      selectInput("group", "Group (HMM iteration)",
                  choices = grps, selected = grps[1])
    }
  })

  # --- renderUI: initial set subgroup checkboxes ---
  output$initialSetSubgroupsUI <- renderUI({
    isets <- initial_set_files()
    if (length(isets) == 0) {
      checkboxGroupInput("initialSetSubgroups", "Subgroups to color",
                         choices = character(0))
    } else {
      checkboxGroupInput("initialSetSubgroups", "Subgroups to color",
                         choices  = names(isets),
                         selected = names(isets))
    }
  })

  # --- Button-gated data loading ---

  results_data <- reactive({
    req(input$group, nchar(input$group) > 0)
    path <- file.path(ISEL_DIR, "iterations", paste0(input$group, ".results"))
    validate(need(file.exists(path),
                  paste("No results file found for group:", input$group)))
    read_results(input$group)
  }) |> bindEvent(input$renderBtn, ignoreNULL = FALSE)

  tree_data <- reactive({
    req(input$treeGroup)
    path <- file.path(TREES_DIR,
                      paste0(input$treeGroup, ".trimmed.aln.treefile"))
    print(path)
    validate(need(file.exists(path),
                  paste("Tree file not found for:", input$treeGroup)))
    tree <- read.tree(file = path)
    tree$tip.label <- unlist(lapply(strsplit(tree$tip.label, "/"), `[`, 1))
    tree
  }) |> bindEvent(input$renderBtn, ignoreNULL = FALSE)

  # --- Joined data reactives ---

  tree_eval_data <- reactive({
    tree <- tree_data()
    df   <- results_data()
    left_join(as_tibble(tree), df, by = "label") %>%
      mutate(evals = ifelse(is.na(evals), 1e-300, evals)) %>%
      distinct(label, .keep_all = TRUE)
  })

  tree_selected_data <- reactive({
    tree  <- tree_data()
    df    <- results_data()
    grp   <- input$group
    fasta_path <- file.path(ISEL_DIR, "iterations", paste0(grp, ".fasta"))
    validate(need(file.exists(fasta_path),
                  paste("No fasta file for group:", grp)))
    selected_ids <- read_fasta_ids(fasta_path)
    selected_df  <- data.frame(label    = selected_ids,
                                selected = "Selected",
                                stringsAsFactors = FALSE)
    left_join(as_tibble(tree), selected_df, by = "label") %>%
      left_join(df, by = "label") %>%
      mutate(selected = ifelse(is.na(selected), "Not selected", selected),
             evals    = ifelse(is.na(evals), 1, evals)) %>%
      distinct(label, .keep_all = TRUE)
  })

  tree_initialsets_data <- reactive({
    tree         <- tree_data()
    tg           <- input$treeGroup
    sel_subgroups <- input$initialSetSubgroups
    req(length(sel_subgroups) > 0)

    # Build initial_df from initialSets/
    initial_df <- data.frame()
    for (sg in sel_subgroups) {
      fp <- file.path(ISEL_DIR, "initialSets",
                      paste0(tg, ".", sg, ".initialSet.txt"))
      if (!file.exists(fp)) next
      ids <- read_initial_set(fp)
      initial_df <- rbind(initial_df,
                          data.frame(label    = ids,
                                     category = paste0("Initial_", sg),
                                     stringsAsFactors = FALSE))
    }

    # Build selected_df from iterations/ fasta (sequences NOT in initial set)
    selected_df <- data.frame()
    for (sg in sel_subgroups) {
      fp <- file.path(ISEL_DIR, "iterations", paste0(tg, ".", sg, ".fasta"))
      if (!file.exists(fp)) next
      ids <- read_fasta_ids(fp)
      ids <- ids[!(ids %in% initial_df$label)]
      if (length(ids) > 0) {
        selected_df <- rbind(selected_df,
                             data.frame(label    = ids,
                                        category = sg,
                                        stringsAsFactors = FALSE))
      }
    }

    combined_df <- rbind(initial_df, selected_df)

    if (nrow(combined_df) > 0) {
      n_occur <- data.frame(table(combined_df$label))
      combined_df <- combined_df %>%
        mutate(category = ifelse(label %in% n_occur$Var1[n_occur$Freq > 1],
                                 "Both", category)) %>%
        distinct()
    }

    left_join(as_tibble(tree), combined_df, by = "label") %>%
      mutate(category = ifelse(is.na(category), "Not selected", category)) %>%
      distinct(label, .keep_all = TRUE)
  })

  # --- Plot reactives (reused by renderPlot and downloadHandler) ---

  plot_eval_density <- reactive({
    df <- results_data()
    req(nrow(df) > 0)
    log_evals <- -log2(df$evals)
    midpoint  <- mean(log_evals)
    d <- density(log_evals)
    ggplot(data.frame(x = d$x, y = d$y), aes(x, y)) +
      geom_line(linewidth = 1.5) +
      geom_segment(aes(xend = x, yend = 0, colour = x), linewidth = 1.5) +
      scale_color_gradient2(name = "-log2(evals)",
                            low = "red", mid = "yellow", high = "darkgreen",
                            midpoint = midpoint) +
      geom_vline(xintercept = midpoint, linetype = "dashed",
                 color = "red", linewidth = 1) +
      labs(title = paste("E-value distribution —", input$group),
           x = "-log2(E-value)", y = "Density") +
      theme_minimal() +
      theme(plot.title  = element_text(hjust = 0.5, size = 60),
            axis.text   = element_text(size = 30),
            axis.title  = element_text(size = 40),
            legend.text = element_text(size = 30),
            legend.title = element_text(size = 40))
  })

  plot_tree_eval <- reactive({
    df_tree  <- tree_eval_data()
    req(nrow(df_tree) > 0)
    real_evals <- results_data()$evals
    log_evals  <- -log2(real_evals[!is.na(real_evals) & real_evals > 0])
    midpoint   <- if (length(log_evals) > 0) mean(log_evals) else 50
    ggtree(tree_data(), color = "lightgrey", layout = input$treeLayout) %<+%
      results_data() +
      geom_tiplab(aes(color = -log2(evals)), size = input$tipLabelSize, show.legend = TRUE) +
      scale_color_gradient2(name = "-log2(E-value)",
                            low = "red", mid = "yellow", high = "darkgreen",
                            midpoint = midpoint) +
      ggtitle(isolate(input$group)) +
      theme(plot.title  = element_text(hjust = 0.5, size = 60),
            axis.text   = element_text(size = 12),
            axis.title  = element_text(size = 14),
            legend.position = "top",
            legend.title = element_text(size = 20),
            legend.text  = element_text(size = 16))
  })

  plot_tree_selected <- reactive({
    df_tree <- tree_selected_data()
    req(nrow(df_tree) > 0)
    df_meta <- df_tree %>% select(label, selected) %>% filter(!is.na(label))
    ggtree(tree_data(), color = "lightgrey", layout = input$treeLayout) %<+%
      df_meta +
      geom_tiplab(aes(color = selected, label = label),
                  size = input$tipLabelSize, show.legend = TRUE) +
      scale_color_manual(name = "Status",
                         values = c("Selected"     = "#35B779FF",
                                    "Not selected" = "#31688EFF")) +
      guides(color = guide_legend(override.aes = list(shape = 15, size = 5))) +
      ggtitle(isolate(input$group)) +
      theme(plot.title  = element_text(hjust = 0.5, size = 60),
            axis.text   = element_text(size = 12),
            axis.title  = element_text(size = 14),
            legend.position = "top",
            legend.title = element_text(size = 20),
            legend.text  = element_text(size = 16))
  })

  plot_tree_initialsets <- reactive({
    df_tree <- tree_initialsets_data()
    req(nrow(df_tree) > 0)
    df_meta <- df_tree %>% select(label, category) %>% filter(!is.na(label))
    ggtree(tree_data(), color = "lightgrey", layout = input$treeLayout) %<+%
      df_meta +
      geom_tiplab(aes(color = category, label = label),
                  size = input$tipLabelSize, show.legend = TRUE) +
      scale_color_discrete(name = "Category") +
      geom_point(size = 0, shape = 15) +
      guides(color = guide_legend(override.aes = list(shape = 15, size = 5))) +
      ggtitle(paste("Initial sets —", isolate(input$treeGroup))) +
      theme(plot.title  = element_text(hjust = 0.5, size = 60),
            axis.text   = element_text(size = 12),
            axis.title  = element_text(size = 14),
            legend.position = "top",
            legend.title = element_text(size = 20),
            legend.text  = element_text(size = 16))
  })

  # --- Render functions ---

  output$plotEvalDensity     <- renderPlot({ plot_eval_density() })
  output$plotTreeEval        <- renderPlot({ plot_tree_eval() })
  output$plotTreeSelected    <- renderPlot({ plot_tree_selected() })
  output$plotTreeInitialSets <- renderPlot({ plot_tree_initialsets() })

  # --- Download handlers ---

  output$dlEvalDensity <- downloadHandler(
    filename = function() paste0(input$group, "_evalue_dist.pdf"),
    content  = function(file) {
      pdf(file, width = input$pdfWidthEvalDensity, height = input$pdfHeightEvalDensity)
      print(plot_eval_density())
      dev.off()
    }
  )

  output$dlTreeEval <- downloadHandler(
    filename = function() paste0(input$treeGroup, "_", input$group, "_tree_eval.pdf"),
    content  = function(file) {
      pdf(file, width = input$pdfWidthTreeEval, height = input$pdfHeightTreeEval)
      print(plot_tree_eval())
      dev.off()
    }
  )

  output$dlTreeSelected <- downloadHandler(
    filename = function() paste0(input$treeGroup, "_", input$group, "_tree_selected.pdf"),
    content  = function(file) {
      pdf(file, width = input$pdfWidthTreeSelected, height = input$pdfHeightTreeSelected)
      print(plot_tree_selected())
      dev.off()
    }
  )

  output$dlTreeInitialSets <- downloadHandler(
    filename = function() paste0(input$treeGroup, "_tree_initialsets.pdf"),
    content  = function(file) {
      pdf(file, width = input$pdfWidthTreeInitialSets, height = input$pdfHeightTreeInitialSets)
      print(plot_tree_initialsets())
      dev.off()
    }
  )

  output$tableResults <- DT::renderDataTable({
    results_data()
  }, options = list(pageLength = 20, scrollX = TRUE))
}

# ===== Section 6: Run =====
shinyApp(ui, server)
