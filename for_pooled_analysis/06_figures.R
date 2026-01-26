# ==============================================================================
# 06_figures.R
# Publication-quality figures for pooled analyses
# ==============================================================================

# setup ------------------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(scales)
library(data.table)
library(UpSetR)

# Ensure output directories exist
if (!dir.exists(here("output", "figures"))) {
  dir.create(here("output", "figures"), recursive = TRUE)
}

# ==============================================================================
# THEME AND COLOR PALETTE
# ==============================================================================

# Colorblind-friendly palette (blue/orange from ColorBrewer)
pal_cancer = c("Non-cancer" = "#0072B2", "Cancer" = "#D55E00")
pal_type   = c("Solid" = "#009E73", "Hematologic" = "#CC79A7")
pal_sites  = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

theme_paper = theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank(),
    strip.background   = element_rect(fill = "gray95", color = NA),
    strip.text         = element_text(face = "bold", size = 10),
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(color = "gray40", size = 9),
    axis.title         = element_text(size = 10),
    axis.text          = element_text(size = 9)
  )

# ==============================================================================
# FIGURE 1: Score Distribution Histograms (Max Score per Encounter)
# ==============================================================================

message("\n== Creating Figure 1: Score Distributions ==")

# Expand maxscores for histogram (using counts)
fig1_data = maxscores_ca_raw[analysis == "main" & !is.na(max_value)]
fig1_data[, ca_lab := fifelse(ca_01 == 1, "Cancer", "Non-cancer")]
fig1_data[, score_lab := factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)]

fig1 = ggplot(fig1_data, aes(x = max_value, weight = n, fill = ca_lab)) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    alpha    = 0.6,
    binwidth = 1,
    color    = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(values = pal_cancer, name = NULL) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_continuous(labels = percent_format()) +
  facet_wrap(~ score_lab, scales = "free", ncol = 3) +
  labs(
    x        = "Maximum Score During Ward Stay",
    y        = "Density",
    title    = "Distribution of Maximum Early Warning Scores",
    subtitle = "Per encounter, stratified by cancer status"
  ) +
  theme_paper +
  theme(legend.position = "top")

ggsave(here("output", "figures", "fig1_score_distributions.pdf"), fig1,
       width = 10, height = 6, dpi = 300)
ggsave(here("output", "figures", "fig1_score_distributions.png"), fig1,
       width = 10, height = 6, dpi = 300)

message("  Figure 1 saved")

# ==============================================================================
# FIGURE 2: AUROC Forest Plot (Main Analysis)
# ==============================================================================

message("\n== Creating Figure 2: AUROC Forest Plot ==")

# Use site-level meta results if available, otherwise weighted
if (exists("auroc_results_final") && nrow(auroc_results_final) > 0) {
  fig2_data = auroc_results_final[analysis == "main"]
} else if (exists("auroc_weighted_final") && nrow(auroc_weighted_final) > 0) {
  fig2_data = auroc_weighted_final[analysis == "main"]
} else if (exists("auroc_meta_final") && nrow(auroc_meta_final) > 0) {
  fig2_data = auroc_meta_final[analysis == "main"]
} else {
  message("  WARNING: No AUROC data available for Figure 2")
  fig2_data = data.table()
}

if (nrow(fig2_data) > 0) {
  
  fig2_data[, ca_lab := fifelse(ca_01 == 1, "Cancer", "Non-cancer")]
  fig2_data[, score_lab := factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)]
  
  fig2 = ggplot(fig2_data, aes(x = auroc, y = score_lab, color = ca_lab)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray60", linewidth = 0.5) +
    geom_vline(xintercept = 0.7, linetype = "dotted", color = "gray80", linewidth = 0.5) +
    geom_pointrange(
      aes(xmin = ci_lower, xmax = ci_upper),
      position = position_dodge(width = 0.6),
      size     = 0.7,
      fatten   = 3
    ) +
    scale_color_manual(values = pal_cancer, name = NULL) +
    scale_x_continuous(
      limits = c(0.5, 0.85),
      breaks = seq(0.5, 0.85, 0.05),
      expand = c(0.02, 0)
    ) +
    facet_wrap(~ metric_lab, ncol = 2) +
    labs(
      x        = "AUROC (95% CI)",
      y        = NULL,
      title    = "Discrimination of Early Warning Scores",
      subtitle = "Random-effects meta-analysis of site-level AUROCs"
    ) +
    theme_paper +
    theme(
      legend.position    = "top",
      panel.grid.major.y = element_blank()
    )
  
  ggsave(here("output", "figures", "fig2_auroc_forest.pdf"), fig2,
         width = 10, height = 5, dpi = 300)
  ggsave(here("output", "figures", "fig2_auroc_forest.png"), fig2,
         width = 10, height = 5, dpi = 300)
  
  message("  Figure 2 saved")
} else {
  message("  Figure 2 skipped - no data")
}

# ==============================================================================
# FIGURE 3: Sensitivity Analysis Forest Plot (Faceted by Score)
# ==============================================================================

message("\n== Creating Figure 3: Sensitivity Analysis Forest Plot ==")

# Get all analysis variants for encounter-level
if (exists("auroc_results_final") && nrow(auroc_results_final) > 0) {
  fig3_data = auroc_results_final[metric == "Encounter max" & !is.na(auroc)]
} else {
  fig3_data = data.table()
}

if (nrow(fig3_data) > 0) {
  
  fig3_data[, ca_lab := fifelse(ca_01 == 1, "Cancer", "Non-cancer")]
  fig3_data[, score_lab := factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)]
  fig3_data[, analysis_lab := factor(analysis, 
                                      levels = c("main", "fullcode_only", "no_ed_req", "win0_96h", "one_enc_per_pt"),
                                      labels = c("Main", "Full-code only", "No ED req.", "0-96h window", "1 enc/patient"))]
  
  fig3 = ggplot(fig3_data, aes(x = auroc, y = analysis_lab, color = ca_lab)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray60", linewidth = 0.4) +
    geom_pointrange(
      aes(xmin = ci_lower, xmax = ci_upper),
      position = position_dodge(width = 0.6),
      size     = 0.5,
      fatten   = 2.5
    ) +
    scale_color_manual(values = pal_cancer, name = NULL) +
    scale_x_continuous(
      limits = c(0.5, 0.85),
      breaks = seq(0.5, 0.85, 0.1)
    ) +
    facet_wrap(~ score_lab, ncol = 3) +
    labs(
      x        = "AUROC (95% CI)",
      y        = NULL,
      title    = "Sensitivity Analyses",
      subtitle = "Encounter-level discrimination across analysis variants"
    ) +
    theme_paper +
    theme(
      legend.position    = "top",
      panel.grid.major.y = element_blank(),
      strip.text         = element_text(size = 9)
    )
  
  ggsave(here("output", "figures", "fig3_sensitivity_forest.pdf"), fig3,
         width = 11, height = 7, dpi = 300)
  ggsave(here("output", "figures", "fig3_sensitivity_forest.png"), fig3,
         width = 11, height = 7, dpi = 300)
  
  message("  Figure 3 saved")
} else {
  message("  Figure 3 skipped - no data")
}

# ==============================================================================
# FIGURE 4: Meta-Analysis Forest Plot (Interaction Terms)
# ==============================================================================

message("\n== Creating Figure 4: Meta-Analysis Forest Plot ==")

if (exists("forest_data_final") && nrow(forest_data_final) > 0) {
  
  fig4_data = forest_data_final[!is.na(or_int)]
  fig4_data[, site_label := fifelse(is_pooled == TRUE, "Pooled", site)]
  
  fig4 = ggplot(fig4_data, aes(x = or_int, y = reorder(site_label, -is_pooled))) +
    geom_vline(xintercept = 1, linetype = "solid", color = "gray40", linewidth = 0.5) +
    geom_vline(xintercept = c(0.98, 1.02), linetype = "dashed", color = "#D55E00", 
               linewidth = 0.4, alpha = 0.7) +
    geom_pointrange(
      aes(xmin = or_int_lower, xmax = or_int_upper,
          color = is_pooled, size = is_pooled),
      fatten = 2
    ) +
    scale_color_manual(values = c("TRUE" = "#0072B2", "FALSE" = "gray50"), guide = "none") +
    scale_size_manual(values = c("TRUE" = 0.8, "FALSE" = 0.5), guide = "none") +
    scale_x_continuous(
      trans  = "log",
      breaks = c(0.9, 0.95, 1, 1.05, 1.1),
      limits = c(0.85, 1.15)
    ) +
    facet_wrap(~ score_lab, ncol = 3) +
    labs(
      x        = "Interaction OR per Score Point (95% CI)",
      y        = NULL,
      title    = "Score × Cancer Interaction",
      subtitle = "Random-effects meta-analysis; dashed lines = equivalence bounds (0.98-1.02)"
    ) +
    theme_paper +
    theme(panel.grid.major.y = element_blank())
  
  ggsave(here("output", "figures", "fig4_meta_forest.pdf"), fig4,
         width = 11, height = 7, dpi = 300)
  ggsave(here("output", "figures", "fig4_meta_forest.png"), fig4,
         width = 11, height = 7, dpi = 300)
  
  message("  Figure 4 saved")
} else {
  message("  WARNING: No meta-analysis data for Figure 4")
}

# ==============================================================================
# FIGURE 5: Cumulative Incidence of Score Positivity
# ==============================================================================

message("\n== Creating Figure 5: Cumulative Incidence ==")

if (exists("cuminc_final") && nrow(cuminc_final) > 0) {
  
  fig5_data = copy(cuminc_final)
  
  # Check which column name is used for score
  if ("score_name" %in% names(fig5_data)) {
    fig5_data[, score_lab := factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)]
  } else if ("score" %in% names(fig5_data)) {
    fig5_data[, score_lab := factor(score, levels = names(SCORE_LABS), labels = SCORE_LABS)]
  } else if ("score_lab" %in% names(fig5_data)) {
    # Already has score_lab
  } else {
    message("  WARNING: Cannot find score column in cuminc_final")
    fig5_data = data.table()
  }
  
  # Check for ca_lab or create it
  if (!"ca_lab" %in% names(fig5_data) && "ca_01" %in% names(fig5_data)) {
    fig5_data[, ca_lab := fifelse(ca_01 == 1, "Cancer", "Non-cancer")]
  }
  
  # Check for outcome_lab or create it
  if (!"outcome_lab" %in% names(fig5_data) && "o_primary_01" %in% names(fig5_data)) {
    fig5_data[, outcome_lab := fifelse(o_primary_01 == 1, "With deterioration", "Without deterioration")]
  }
}

if (exists("fig5_data") && nrow(fig5_data) > 0 && "score_lab" %in% names(fig5_data)) {
  
  fig5 = ggplot(fig5_data, 
                aes(x = time_bin_start / 24, y = cum_inc, 
                    color = ca_lab, linetype = outcome_lab)) +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
    scale_x_continuous(breaks = seq(0, 14, 2)) +
    scale_color_manual(values = pal_cancer, name = NULL) +
    scale_linetype_manual(
      values = c("With deterioration" = "solid", "Without deterioration" = "dashed"),
      name   = NULL
    ) +
    facet_wrap(~ score_lab, ncol = 3) +
    labs(
      x        = "Days from Ward Admission",
      y        = "Cumulative Proportion Ever Positive",
      title    = "Time to Score Positivity",
      subtitle = "Stratified by cancer status and outcome"
    ) +
    theme_paper +
    theme(legend.position = "bottom")
  
  ggsave(here("output", "figures", "fig5_cumulative_incidence.pdf"), fig5,
         width = 10, height = 7, dpi = 300)
  ggsave(here("output", "figures", "fig5_cumulative_incidence.png"), fig5,
         width = 10, height = 7, dpi = 300)
  
  message("  Figure 5 saved")
} else {
  message("  WARNING: No cumulative incidence data for Figure 5")
}

# ==============================================================================
# FIGURE 6: Site-Level Heterogeneity
# ==============================================================================

message("\n== Creating Figure 6: Site Heterogeneity ==")

if (exists("auroc_site_enc") && nrow(auroc_site_enc) > 0) {
  
  fig6_data = auroc_site_enc[analysis == "main" & !is.na(auroc)]
  fig6_data[, ca_lab := fifelse(ca_01 == 1, "Cancer", "Non-cancer")]
  fig6_data[, score_lab := factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)]
  fig6_data[, site_upper := toupper(site_clean)]
  
  # Order sites by mean AUROC
  site_order = fig6_data[, .(mean_auc = mean(auroc, na.rm = TRUE)), by = site_upper]
  setorder(site_order, -mean_auc)
  fig6_data[, site_upper := factor(site_upper, levels = site_order$site_upper)]
  
  fig6 = ggplot(fig6_data, aes(x = site_upper, y = auroc, color = ca_lab)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60") +
    geom_pointrange(
      aes(ymin = ci_lower, ymax = ci_upper),
      position = position_dodge(width = 0.6),
      size     = 0.4,
      fatten   = 2
    ) +
    scale_color_manual(values = pal_cancer, name = NULL) +
    scale_y_continuous(limits = c(0.45, 0.9), breaks = seq(0.5, 0.9, 0.1)) +
    facet_wrap(~ score_lab, ncol = 3) +
    labs(
      x        = "Site",
      y        = "AUROC (95% CI)",
      title    = "Site-Level Discrimination",
      subtitle = "Encounter-level AUROCs with DeLong confidence intervals"
    ) +
    theme_paper +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "top"
    )
  
  ggsave(here("output", "figures", "fig6_site_heterogeneity.pdf"), fig6,
         width = 12, height = 8, dpi = 300)
  ggsave(here("output", "figures", "fig6_site_heterogeneity.png"), fig6,
         width = 12, height = 8, dpi = 300)
  
  message("  Figure 6 saved")
} else {
  message("  WARNING: No site-level AUROC data for Figure 6")
}

# ==============================================================================
# FIGURE 7: UpSet Plots (Score Combinations)
# ==============================================================================

message("\n== Creating Figure 7: UpSet Plots ==")

if (exists("upset_raw") && nrow(upset_raw) > 0) {
  
  # Prepare data for UpSetR - need binary matrix format
  upset_data = copy(upset_raw)
  
  # Cancer patients
  upset_ca = upset_data[ca_01 == 1, .(
    n = sum(n, na.rm = TRUE)
  ), by = .(sirs, qsofa, mews, news, mews_sf)]
  
  # Expand for UpSetR (it needs individual rows or use fromExpression)
  # Create expression list for UpSetR
  make_upset_expression = function(dt) {
    expr_list = list()
    for (i in seq_len(nrow(dt))) {
      row = dt[i]
      scores_positive = c()
      if (row$sirs == 1)    scores_positive = c(scores_positive, "SIRS")
      if (row$qsofa == 1)   scores_positive = c(scores_positive, "qSOFA")
      if (row$mews == 1)    scores_positive = c(scores_positive, "MEWS")
      if (row$news == 1)    scores_positive = c(scores_positive, "NEWS")
      if (row$mews_sf == 1) scores_positive = c(scores_positive, "MEWS-SF")
      
      if (length(scores_positive) > 0) {
        combo_name = paste(scores_positive, collapse = "&")
        if (combo_name %in% names(expr_list)) {
          expr_list[[combo_name]] = expr_list[[combo_name]] + row$n
        } else {
          expr_list[[combo_name]] = row$n
        }
      }
    }
    return(expr_list)
  }
  
  expr_ca = make_upset_expression(upset_ca)
  
  # Save UpSet plot for cancer patients
  pdf(here("output", "figures", "fig7a_upset_cancer.pdf"), width = 10, height = 6)
  print(upset(fromExpression(expr_ca),
        sets = c("SIRS", "qSOFA", "MEWS", "NEWS", "MEWS-SF"),
        order.by = "freq",
        decreasing = TRUE,
        mb.ratio = c(0.6, 0.4),
        number.angles = 0,
        text.scale = 1.2,
        point.size = 3,
        line.size = 1,
        mainbar.y.label = "Number of Encounters",
        sets.x.label = "Ever Positive",
        main.bar.color = "#D55E00",
        sets.bar.color = "#0072B2"))
  dev.off()
  
  # Non-cancer patients
  upset_nonca = upset_data[ca_01 == 0, .(
    n = sum(n, na.rm = TRUE)
  ), by = .(sirs, qsofa, mews, news, mews_sf)]
  
  expr_nonca = make_upset_expression(upset_nonca)
  
  pdf(here("output", "figures", "fig7b_upset_noncancer.pdf"), width = 10, height = 6)
  print(upset(fromExpression(expr_nonca),
        sets = c("SIRS", "qSOFA", "MEWS", "NEWS", "MEWS-SF"),
        order.by = "freq",
        decreasing = TRUE,
        mb.ratio = c(0.6, 0.4),
        number.angles = 0,
        text.scale = 1.2,
        point.size = 3,
        line.size = 1,
        mainbar.y.label = "Number of Encounters",
        sets.x.label = "Ever Positive",
        main.bar.color = "#0072B2",
        sets.bar.color = "#0072B2"))
  dev.off()
  
  # PNG versions
  png(here("output", "figures", "fig7a_upset_cancer.png"), width = 10, height = 6, 
      units = "in", res = 300)
  print(upset(fromExpression(expr_ca),
        sets = c("SIRS", "qSOFA", "MEWS", "NEWS", "MEWS-SF"),
        order.by = "freq",
        decreasing = TRUE,
        mb.ratio = c(0.6, 0.4),
        number.angles = 0,
        text.scale = 1.2,
        point.size = 3,
        line.size = 1,
        mainbar.y.label = "Number of Encounters",
        sets.x.label = "Ever Positive",
        main.bar.color = "#D55E00",
        sets.bar.color = "#0072B2"))
  dev.off()
  
  png(here("output", "figures", "fig7b_upset_noncancer.png"), width = 10, height = 6,
      units = "in", res = 300)
  print(upset(fromExpression(expr_nonca),
        sets = c("SIRS", "qSOFA", "MEWS", "NEWS", "MEWS-SF"),
        order.by = "freq",
        decreasing = TRUE,
        mb.ratio = c(0.6, 0.4),
        number.angles = 0,
        text.scale = 1.2,
        point.size = 3,
        line.size = 1,
        mainbar.y.label = "Number of Encounters",
        sets.x.label = "Ever Positive",
        main.bar.color = "#0072B2",
        sets.bar.color = "#0072B2"))
  dev.off()
  
  message("  Figure 7 (UpSet plots) saved")
} else {
  message("  WARNING: No upset data for Figure 7")
}

# ==============================================================================
# FIGURE S1: Hematologic vs Solid Comparison
# ==============================================================================

message("\n== Creating Figure S1: Cancer Type Comparison ==")

if (exists("liquid_aurocs_final") && nrow(liquid_aurocs_final) > 0) {
  
  figS1_data = copy(liquid_aurocs_final)
  figS1_data[, score_lab := factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)]
  
  figS1 = ggplot(figS1_data, aes(x = auroc, y = score_lab, color = cancer_type)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray60") +
    geom_pointrange(
      aes(xmin = ci_lower, xmax = ci_upper),
      position = position_dodge(width = 0.5),
      size     = 0.7,
      fatten   = 3
    ) +
    scale_color_manual(values = pal_type, name = "Cancer Type") +
    scale_x_continuous(limits = c(0.5, 0.85), breaks = seq(0.5, 0.85, 0.05)) +
    labs(
      x        = "AUROC (95% CI)",
      y        = NULL,
      title    = "Hematologic vs Solid Malignancies",
      subtitle = "Encounter-level discrimination"
    ) +
    theme_paper +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position    = "top"
    )
  
  ggsave(here("output", "figures", "figS1_cancer_type.pdf"), figS1,
         width = 8, height = 5, dpi = 300)
  ggsave(here("output", "figures", "figS1_cancer_type.png"), figS1,
         width = 8, height = 5, dpi = 300)
  
  message("  Figure S1 saved")
} else {
  message("  WARNING: No liquid/solid data for Figure S1")
}

# ==============================================================================
# FIGURE S2: Slope Comparison (Score-Outcome Relationship)
# ==============================================================================

message("\n== Creating Figure S2: Slope Comparison ==")

if (exists("slope_preds_final") && nrow(slope_preds_final) > 0) {
  
  figS2 = ggplot(slope_preds_final, 
                 aes(x = score_value, y = predicted_prob,
                     color = cancer_status, fill = cancer_status)) +
    geom_ribbon(aes(ymin = predicted_prob - 0.02, ymax = predicted_prob + 0.02),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 1) +
    scale_y_continuous(labels = percent_format(), limits = c(0, NA)) +
    scale_color_manual(values = pal_cancer, name = NULL) +
    scale_fill_manual(values = pal_cancer, name = NULL) +
    facet_wrap(~ score_label, scales = "free_x", ncol = 3) +
    labs(
      x        = "Score Value",
      y        = "Predicted Probability of Deterioration",
      title    = "Score-Outcome Relationship",
      subtitle = "Predicted probability from pooled logistic regression"
    ) +
    theme_paper +
    theme(legend.position = "bottom")
  
  ggsave(here("output", "figures", "figS2_slopes.pdf"), figS2,
         width = 10, height = 7, dpi = 300)
  ggsave(here("output", "figures", "figS2_slopes.png"), figS2,
         width = 10, height = 7, dpi = 300)
  
  message("  Figure S2 saved")
} else {
  message("  WARNING: No slope data for Figure S2")
}

# ==============================================================================
# FIGURE S3: Comparison of AUROC Methods (Meta vs Weighted)
# ==============================================================================

message("\n== Creating Figure S3: Method Comparison ==")

if (exists("auroc_comparison_final") && nrow(auroc_comparison_final) > 0) {
  
  figS3_data = auroc_comparison_final[analysis == "main" & !is.na(auroc_meta) & !is.na(auroc_weighted)]
  figS3_data[, score_lab := factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)]
  figS3_data[, ca_lab := fifelse(ca_01 == 1, "Cancer", "Non-cancer")]
  
  figS3 = ggplot(figS3_data, aes(x = auroc_weighted, y = auroc_meta, color = ca_lab)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_point(size = 3, alpha = 0.8) +
    geom_errorbar(aes(ymin = ci_lo_meta, ymax = ci_hi_meta), width = 0, alpha = 0.5) +
    geom_errorbarh(aes(xmin = ci_lo_weighted, xmax = ci_hi_weighted), height = 0, alpha = 0.5) +
    scale_color_manual(values = pal_cancer, name = NULL) +
    coord_equal(xlim = c(0.55, 0.8), ylim = c(0.55, 0.8)) +
    facet_wrap(~ metric, ncol = 2) +
    labs(
      x        = "Pooled Weighted AUC (Hanley-McNeil CI)",
      y        = "Meta-analyzed Site AUROCs (Random Effects)",
      title    = "Comparison of AUROC Estimation Methods",
      subtitle = "Main analysis; diagonal = perfect agreement"
    ) +
    theme_paper +
    theme(legend.position = "top")
  
  ggsave(here("output", "figures", "figS3_method_comparison.pdf"), figS3,
         width = 9, height = 5, dpi = 300)
  ggsave(here("output", "figures", "figS3_method_comparison.png"), figS3,
         width = 9, height = 5, dpi = 300)
  
  message("  Figure S3 saved")
} else {
  message("  WARNING: No comparison data for Figure S3")
}

# ==============================================================================
# COMPLETE
# ==============================================================================

message("\n== Figure generation complete ==")

figures_manifest = data.table(
  figure = c("Fig 1", "Fig 2", "Fig 3", "Fig 4", "Fig 5", "Fig 6", "Fig 7a/b",
             "Fig S1", "Fig S2", "Fig S3"),
  content = c(
    "Score distribution histograms",
    "AUROC forest plot (main analysis)",
    "Sensitivity analyses forest plot (faceted by score)",
    "Meta-analysis forest plot (interaction terms)",
    "Cumulative incidence of positivity",
    "Site-level heterogeneity",
    "UpSet plots (cancer/non-cancer)",
    "Hematologic vs solid comparison",
    "Score-outcome slopes",
    "AUROC method comparison (meta vs weighted)"
  ),
  file_base = c(
    "fig1_score_distributions",
    "fig2_auroc_forest",
    "fig3_sensitivity_forest",
    "fig4_meta_forest",
    "fig5_cumulative_incidence",
    "fig6_site_heterogeneity",
    "fig7a/b_upset",
    "figS1_cancer_type",
    "figS2_slopes",
    "figS3_method_comparison"
  )
)

message("\n  Figures saved to: ", here("output", "figures"))
print(figures_manifest)
