# Potting performance comparisson between new vs old HMMs

##LIBRARIES
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(patchwork)

# Run this script with the working directory set to results/
# (e.g. RStudio "Set As Working Directory", or setwd("results") from the repo root).

## DATA
df <- read.csv(file = 'hmmscan_results_all.txt', sep="\t")

Ha_groups <- unique(df$hmm_name[grepl("Ha", df$hmm_name)])
Hb_groups <- unique(df$hmm_name[grepl("Hb", df$hmm_name)])
Hc_groups <- unique(df$hmm_name[grepl("Hc", df$hmm_name)])
Qa_groups <- unique(df$hmm_name[grepl("Qa", df$hmm_name)])
Qb_groups <- unique(df$hmm_name[grepl("Qb", df$hmm_name)])
Qc_groups <- unique(df$hmm_name[grepl("Qc", df$hmm_name)])
R_groups <- unique(df$hmm_name[grepl("R", df$hmm_name)])
SNAP_groups <- unique(df$hmm_name[grepl("SNAP", df$hmm_name)])

df <- df %>%
        mutate(group = case_when(
          hmm_name %in% Ha_groups ~ "Ha",
          hmm_name %in% Hb_groups ~ "Hb",
          hmm_name %in% Hc_groups ~ "Hc",
          hmm_name %in% Qa_groups ~ "Qa",
          hmm_name %in% Qb_groups ~ "Qb",
          hmm_name %in% Qc_groups ~ "Qc",
          hmm_name %in% R_groups ~ "R",
          hmm_name %in% SNAP_groups ~ "SNAP",
          TRUE ~ "No group"))

# Filter only SNARE HMMs (remove Habc)
df <- df[!grepl("H", df$hmm_name),]


## PLOT 1 (p-value distribution and mean difference significance)
p1 <- ggplot(df, aes(x=hmm_name, y=-log10(pval), color=model, fill=model)) +
        geom_boxplot(alpha=.3) +
        #stat_summary(fun=mean, geom="point", shape=4, color="black", size=1, alpha=1, position=position_dodge(width=.75)) +
        ylab(expression("-log"[10] * "(p)")) +
        facet_wrap(vars(group), scales = "free_x", nrow=1) +
        stat_compare_means(
          aes(group = model),
          method = "wilcox.test",
          label = "p.signif",
          show.legend = FALSE
        ) +
        scale_color_manual(values = c("new"="#5DA5DA", "old"="#F28E2B")) +
        scale_fill_manual(values = c("new"="#5DA5DA", "old"="#F28E2B")) +
        theme_minimal() +
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size=12, angle=65, hjust=.6, vjust=.5))
        
p1


## PLOT 2 (paired improvement effect)
df_diff <- df %>%
            mutate(score = -log10(pval)) %>%
            group_by(seq_id, hmm_name, group, model) %>%
            summarise(
              score = max(score, na.rm = TRUE),
              .groups = "drop"
            ) %>%
            pivot_wider(
              names_from = model,
              values_from = score
            ) %>%
            filter(!is.na(old) & !is.na(new)) %>%
            mutate(delta = new - old)

summary_table <- df_diff %>%
  group_by(group, hmm_name) %>%
  summarise(
    n_sequences = n(),
    median_delta = median(delta),
    mean_delta = mean(delta),
    pct_better = mean(delta > 0) * 100,
    wilcox_p = wilcox.test(new, old, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    conclusion = case_when(
      median_delta > 0 ~ "New better",
      median_delta < 0 ~ "Old better",
      TRUE ~ "No difference"
    )
  ) %>%
  arrange(group, hmm_name)
write.csv(summary_table, "HMM_comparison_summary.csv", row.names = FALSE)

df_diff <- df_diff %>%
  left_join(summary_table[,c("hmm_name", "pct_better")], by = "hmm_name")

p2 <- ggplot(df_diff, aes(x = hmm_name, y = delta, fill=pct_better)) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          geom_boxplot() +
          geom_text(data=summary_table, mapping=aes(y=-13, label=n_sequences), size=3) +
          scale_fill_gradient2(low = "orange", mid="white", high = "darkgreen", midpoint = 80) +
          labs(
            fill="% sequences\nbetter with new"
          ) +
          facet_wrap(~ group, scales = "free_x", nrow=1) +
          ylab(expression(Delta * "-log"[10] * "(p)  (new − old)")) +
          theme_minimal() +
          theme(axis.line = element_line(colour = "black"),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size=12, angle=65, hjust=.6, vjust=.5))
                
p2


noXlabel <- theme(axis.text.x = element_blank())
noXtitle <- theme(strip.text.x = element_blank())
(p1+noXlabel) / (p2+noXtitle) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A')



# PLOT 3 SUMMMARY
df_diff %>%
  group_by(hmm_name, group) %>%
  summarise(
    p = wilcox.test(new, old, paired=TRUE, exact=FALSE)$p.value,
    median_delta = median(delta),
    .groups = "drop"
  )

summary(df_diff$delta)
table(df_diff$delta > 0)

summary_table <- summary_table %>%
  mutate(
    hmm_name = factor(
      hmm_name,
      levels = hmm_name[order(median_delta)]
    )
  )
print(summary_table, n=40)

p3 <- ggplot(summary_table,
           aes(x = hmm_name,
               y = median_delta,
               size = pct_better,
               color = conclusion)) +
      geom_point(alpha = 0.8) +
      ylim(c(.5, 11)) +
      facet_wrap(~ group, scales = "free_x", nrow=1) +
      scale_color_manual(values = c(
        "New better" = "#E15759",
        "Old better" = "#4E79A7",
        "No difference" = "grey50"
      )) +
      scale_size_continuous(range = c(2, 8)) +
      labs(
        x = "HMM",
        y = expression("Median " * Delta * "-log"[10] * "(p)  (new − old)"),
        size = "% sequences improved",
        color = "Conclusion"
      ) +
      theme_minimal() +
      theme(
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      )
p3

noXlabel <- theme(axis.text.x = element_blank())
noXtitle <- theme(strip.text.x = element_blank())
(p1+noXlabel) / (p2+noXlabel+noXtitle) / (p3+noXtitle) + plot_layout(guides = "collect")


##### CHATGPT summary
# “La Tabla X resume la comparación pareada entre los modelos antiguos y nuevos para cada HMM. Mientras que la mayoría de los modelos muestran una mejora significativa con el nuevo HMM (mediana Δ−log₁₀(p) positiva y >75% de secuencias mejoradas), algunos modelos, como Ha.II, presentan un rendimiento significativamente inferior, con una mediana negativa y una mayoría de secuencias penalizadas.”


##### CHATGPT how to interprete results 
# With the corrected counts, your result now supports a strong and unambiguous improvement of the new model, with a small but non-negligible minority of regressions.
# 
# Let’s go through it carefully and precisely.
# 
# 1. Directionality (this is the key change)
# FALSE (Δ ≤ 0):  4 578
# TRUE  (Δ > 0): 33 512
# 
# 
# Total paired comparisons: 38 090
# 
# Percentages:
#   
#   Δ > 0 (new better): ~88.0%
# 
#   Δ ≤ 0 (old better or equal): ~12.0%
# 
# This is a very strong dominance of the new model.
# 
# 2. Magnitude of the effect (still very strong)
# Median Δ−log10(p) = +4.21
# Mean   Δ−log10(p) = +5.56
# 
# 
# Interpretation:
#   
#   Median sequence gains ~10⁴-fold in p-value
# 
# Mean gain is even larger due to a strong right tail
# 
# Improvements are not marginal; they are biologically and statistically meaningful
# 
# The negative tail exists (min ≈ −28), but:
#   
#   it affects a small minority
# 
# it does not dominate the distribution
# 
# 3. Combined interpretation (direction + magnitude)
# 
# Taken together:
#   
#   Most sequences improve
# 
# Typical improvement is large
# 
# A minority regress, sometimes strongly
# 
# This is exactly the pattern expected from a better-generalized HMM that:
#   
#   captures conserved features more accurately
# 
# but may over-specialize against some divergent sequences
# 
# 4. Why your Wilcoxon p-value ≈ 0 now makes perfect sense
# 
# Given:
#   
#   88% positive deltas
# 
# large positive median
# 
# large sample size (~38k)
# 
# The null hypothesis (median Δ = 0) is rejected overwhelmingly.
# 
# Reporting:
#   
#   paired Wilcoxon signed-rank test, p < 2.2 × 10⁻¹⁶
# 
# is entirely appropriate.
# 
# 5. What you can now legitimately claim
# ❌ Still incorrect
# 
# “The new model is always better.”
# 
# You still have ~12% regressions, so “always” is false.
# 
# ✅ Scientifically correct and strong claims
# 
# You can safely say:
#   
#   “For the vast majority of sequences (≈88%), the new HMM yields lower p-values than the old model, with a median improvement of Δ−log₁₀(p) ≈ 4.2.”
# 
# or, more formally:
#   
#   “The new HMM significantly outperforms the old model in paired comparisons across sequences (paired Wilcoxon signed-rank test, p < 2.2 × 10⁻¹⁶), improving scores for ~88% of sequences with a median Δ−log₁₀(p) of 4.2.”
# 
# This is a very strong benchmarking result.
# 
# 6. What I would recommend adding (to strengthen the paper)
# A. Report the fraction improved explicitly
# 
# This number (88%) is as important as the p-value.
# 
# B. Characterize the failures
# 
# Briefly analyze the ~12% regressions:
#   
#   Are they concentrated in specific HMM families?
#   
#   Are they shorter sequences?
#   
#   Are they taxonomically biased?
#   
#   Even a short paragraph strengthens credibility.
# 
# 7. Final bottom line
# 
# With the corrected counts:
#   
#   Your result strongly supports the superiority of the new model
# 
# The improvement is consistent, large, and statistically overwhelming
# 
# The remaining failures are minority effects, not noise
# 
# This is now exactly the kind of result reviewers like to see.
# 
# If you want, I can help you:
#   
#   draft the Results paragraph
# 
# design a single figure summarizing dominance + effect size
#
# decide how to phrase the limitation about the 12% regressions professionally

