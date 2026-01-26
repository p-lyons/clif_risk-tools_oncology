# ==============================================================================
# 06_figures.R
# Publication-ready figure generation for pooled analyses
# Assumes data objects from 00_load through 05_subgroups are already in environment
# ==============================================================================

# setup ------------------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(scales)
library(data.table)
library(tidytable)

# output directory -------------------------------------------------------------

fig_dir = here::here("output", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

message("\n== Figure output directory: ", fig_dir, " ==")

# theme setup ------------------------------------------------------------------

theme_paper = theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.border       = element_rect(color = "gray50", fill = NA, linewidth = 0.5),
    strip.text         = element_text(face = "bold", size = 11),
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle      = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title         = element_text(face = "bold", size = 11)
  )

# color palettes
colors_cancer = c("Non-cancer" = "#3498DB", "Cancer" = "#E74C3C")
colors_type   = c("Solid" = "#7570b3", "Hematologic" = "#e7298a")

# ==============================================================================
# FIGURE 1: AUROC Forest Plot (Main Result)
# ==============================================================================

message("\n== Creating Figure 1: AUROC Forest Plot ==")

fig1_data = auroc_results_final[analysis == "main"]

# separate by metric for cleaner dodging
fig1_enc = fig1_data[metric == "Encounter max"]
fig1_24h = fig1_data[metric == "24-hour horizon"]

fig1 = ggplot(fig1_data, aes(x = auroc, y = score_lab, color = ca_lab)) +

  geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.4, color = "gray50") +
  geom_pointrange(
    aes(xmin = ci_lower, xmax = ci_upper),
    position = position_dodge(width = 0.6),
    size     = 0.7,
    fatten   = 3
  ) +
  scale_color_manual(values = colors_cancer, name = NULL) +
  scale_x_continuous(
    limits = c(0.55, 0.80), 
    breaks = seq(0.55, 0.80, 0.05),
    expand = c(0, 0)
  ) +
  facet_wrap(~ metric_lab, ncol = 2) +
  labs(
    x     = "AUROC (95% CI)",
    y     = NULL,
    title = "Discrimination Performance by Cancer Status"
  ) +
  theme_paper +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position    = "top"
  )

ggsave(file.path(fig_dir, "fig1_auroc_forest.pdf"), fig1, width = 9, height = 5)
ggsave(file.path(fig_dir, "fig1_auroc_forest.png"), fig1, width = 9, height = 5, dpi = 300)

message("  Figure 1 saved")

# ==============================================================================
# FIGURE 2: Sensitivity/Specificity at Standard Thresholds
# ==============================================================================

message("\n== Creating Figure 2: Sensitivity/Specificity ==")

fig2_data = copy(sesp_final)
fig2_data[, `:=`(
  ca_lab    = fifelse(ca_01 == 1, "Cancer", "Non-cancer"),
  score_lab = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
)]

# reshape for faceting
fig2_long = melt(
  fig2_data,
  id.vars       = c("score_lab", "ca_lab"),
  measure.vars  = c("sensitivity", "specificity"),
  variable.name = "metric",
  value.name    = "value"
)
fig2_long[, metric := fifelse(metric == "sensitivity", "Sensitivity", "Specificity")]

fig2 = ggplot(fig2_long, aes(x = score_lab, y = value, fill = ca_lab)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(
    aes(label = sprintf("%.0f%%", value * 100)),
    position = position_dodge(width = 0.8),
    vjust    = -0.5,
    size     = 3,
    fontface = "bold"
  ) +
  scale_fill_manual(values = colors_cancer, name = NULL) +
  scale_y_continuous(
    labels = percent_format(), 
    limits = c(0, 1.05),
    expand = c(0, 0)
  ) +
  facet_wrap(~ metric) +
  labs(
    x        = NULL, 
    y        = NULL,
    title    = "Score Performance at Standard Thresholds",
    subtitle = "SIRS ≥2, qSOFA ≥2, MEWS ≥5, NEWS ≥5, MEWS-SF ≥7"
  ) +
  theme_paper +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(fig_dir, "fig2_sens_spec.pdf"), fig2, width = 9, height = 5)
ggsave(file.path(fig_dir, "fig2_sens_spec.png"), fig2, width = 9, height = 5, dpi = 300)

message("  Figure 2 saved")

# ==============================================================================
# FIGURE 3: Meta-analysis Forest Plot (Interaction Terms)
# ==============================================================================

message("\n== Creating Figure 3: Meta-analysis Forest Plot ==")

fig3_data = forest_data_final[!is.na(or_int)]

# order sites by whether pooled
fig3_data[, site_label := fifelse(is_pooled, "Pooled", toupper(site))]

fig3 = ggplot(fig3_data, aes(x = or_int, y = reorder(site_label, is_pooled), 
                              color = is_pooled, size = is_pooled)) +
  geom_vline(xintercept = 1, linetype = "solid", alpha = 0.3) +
  geom_vline(xintercept = c(0.98, 1.02), linetype = "dashed", alpha = 0.4, color = "#E74C3C") +
  geom_pointrange(aes(xmin = or_int_lower, xmax = or_int_upper)) +
  scale_color_manual(values = c("TRUE" = "#2C3E50", "FALSE" = "#7F8C8D"), guide = "none") +
  scale_size_manual(values = c("TRUE" = 1.0, "FALSE" = 0.5), guide = "none") +
  scale_x_continuous(
    trans  = "log",
    breaks = c(0.90, 0.95, 1.0, 1.05, 1.10),
    limits = c(0.88, 1.12)
  ) +
  facet_wrap(~ score_lab, ncol = 3) +
  labs(
    x        = "Interaction OR per Score Point (95% CI)",
    y        = NULL,
    title    = "Score × Cancer Interaction by Site",
    subtitle = "Red dashed lines indicate equivalence bounds (0.98–1.02)"
  ) +
  theme_paper +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text         = element_text(size = 10)
  )

ggsave(file.path(fig_dir, "fig3_meta_forest.pdf"), fig3, width = 10, height = 8)
ggsave(file.path(fig_dir, "fig3_meta_forest.png"), fig3, width = 10, height = 8, dpi = 300)

message("  Figure 3 saved")

# ==============================================================================
# FIGURE 4: Cumulative Incidence of Score Positivity
# ==============================================================================

message("\n== Creating Figure 4: Cumulative Incidence ==")

fig4_data = cuminc_final[o_primary_01 == 0]  # patients without deterioration

fig4 = ggplot(fig4_data, aes(x = time_bin_start / 24, y = cum_inc, color = ca_lab)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = colors_cancer, name = NULL) +
  scale_y_continuous(
    labels = percent_format(),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_continuous(
    limits = c(0, 7),
    breaks = 0:7
  ) +
  facet_wrap(~ score_lab, ncol = 3, scales = "free_y") +
  labs(
    x        = "Days Since Ward Admission",
    y        = "Cumulative Incidence",
    title    = "Cumulative Incidence of Score Positivity",
    subtitle = "Among patients who did not deteriorate"
  ) +
  theme_paper +
  theme(
    legend.position = "top",
    panel.spacing   = unit(1, "lines")
  )

ggsave(file.path(fig_dir, "fig4_cumulative_incidence.pdf"), fig4, width = 10, height = 7)
ggsave(file.path(fig_dir, "fig4_cumulative_incidence.png"), fig4, width = 10, height = 7, dpi = 300)

message("  Figure 4 saved")

# ==============================================================================
# FIGURE 5: Slope Comparison (Score-Outcome Relationship)
# ==============================================================================

message("\n== Creating Figure 5: Slope Comparison ==")

if (exists("slope_preds_final") && nrow(slope_preds_final) > 0) {
  
  fig5 = ggplot(slope_preds_final, aes(x = score_value, y = predicted_prob,
                                        color = cancer_status, linetype = cancer_status)) +
    geom_line(linewidth = 1.2) +
    scale_y_continuous(labels = percent_format(), limits = c(0, NA)) +
    scale_color_manual(
      values = c("Non-cancer" = "#3498DB", "Cancer" = "#E74C3C"),
      name   = NULL
    ) +
    scale_linetype_manual(
      values = c("Non-cancer" = "solid", "Cancer" = "dashed"),
      name   = NULL
    ) +
    facet_wrap(~ score_label, scales = "free_x", ncol = 3) +
    labs(
      x        = "Score Value",
      y        = "Predicted Probability of Deterioration",
      title    = "Score-Outcome Relationships by Cancer Status",
      subtitle = "Pooled across all sites"
    ) +
    theme_paper +
    theme(
      legend.position = "top",
      strip.text      = element_text(size = 9)
    )
  
  ggsave(file.path(fig_dir, "fig5_slopes.pdf"), fig5, width = 10, height = 7)
  ggsave(file.path(fig_dir, "fig5_slopes.png"), fig5, width = 10, height = 7, dpi = 300)
  
  message("  Figure 5 saved")
  
} else {
  message("  WARNING: No slope data available for Figure 5")
}

# ==============================================================================
# FIGURE S1: Inclusion Flow Diagram
# ==============================================================================

message("\n== Creating Figure S1: Flow Diagram ==")

if (exists("flow_diagram_final")) {
  
  # save using DiagrammeRsvg
  library(DiagrammeRsvg)
  
  export_svg(flow_diagram_final) |>
    charToRaw() |>
    rsvg::rsvg_pdf(file.path(fig_dir, "figS1_flow.pdf"))
  
  export_svg(flow_diagram_final) |>
    charToRaw() |>
    rsvg::rsvg_png(file.path(fig_dir, "figS1_flow.png"), width = 1200)
  
  message("  Figure S1 saved")
  
} else {
  message("  WARNING: No flow diagram available for Figure S1")
}

# ==============================================================================
# FIGURE S2: Site-Level Heterogeneity
# ==============================================================================

message("\n== Creating Figure S2: Site Heterogeneity ==")

figS2 = ggplot(site_aurocs_final, aes(x = site_ordered, y = auroc, color = ca_lab)) +
  geom_pointrange(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = position_dodge(width = 0.6),
    size     = 0.4,
    fatten   = 2
  ) +
  geom_hline(
    data     = site_heterogeneity_final,
    aes(yintercept = mean_auroc, color = ca_lab),
    linetype = "dashed",
    alpha    = 0.6
  ) +
  scale_color_manual(values = colors_cancer, name = NULL) +
  scale_y_continuous(limits = c(0.45, 0.85), breaks = seq(0.45, 0.85, 0.1)) +
  facet_wrap(~ score_lab, ncol = 3) +
  labs(
    x        = "Site (ordered by mean rank across scores)",
    y        = "AUROC (95% CI)",
    title    = "Site-Level Score Performance",
    subtitle = "Dashed lines indicate pooled mean"
  ) +
  theme_paper +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "top"
  )

ggsave(file.path(fig_dir, "figS2_site_heterogeneity.pdf"), figS2, width = 11, height = 8)
ggsave(file.path(fig_dir, "figS2_site_heterogeneity.png"), figS2, width = 11, height = 8, dpi = 300)

message("  Figure S2 saved")

# ==============================================================================
# FIGURE S3: Hematologic vs Solid Tumor Comparison
# ==============================================================================

message("\n== Creating Figure S3: Cancer Type Comparison ==")

if (exists("liquid_aurocs_final") && nrow(liquid_aurocs_final) > 0) {
  
  figS3 = ggplot(liquid_aurocs_final, aes(x = auroc, y = score_lab, color = cancer_type)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.4) +
    geom_pointrange(
      aes(xmin = ci_lower, xmax = ci_upper),
      position = position_dodge(width = 0.5),
      size     = 0.7,
      fatten   = 3
    ) +
    scale_color_manual(values = colors_type, name = "Cancer Type") +
    scale_x_continuous(limits = c(0.50, 0.80), breaks = seq(0.50, 0.80, 0.05)) +
    labs(
      x        = "AUROC (95% CI)",
      y        = NULL,
      title    = "Score Performance: Hematologic vs Solid Tumors"
    ) +
    theme_paper +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position    = "top"
    )
  
  ggsave(file.path(fig_dir, "figS3_cancer_type.pdf"), figS3, width = 8, height = 5)
  ggsave(file.path(fig_dir, "figS3_cancer_type.png"), figS3, width = 8, height = 5, dpi = 300)
  
  message("  Figure S3 saved")
  
} else {
  message("  WARNING: No liquid/solid data available for Figure S3")
}

# ==============================================================================
# FIGURE S4: AUROC Difference Heatmap (Sensitivity Analyses)
# ==============================================================================

message("\n== Creating Figure S4: AUROC Difference Heatmap ==")

figS4_data = unique(auroc_results_final[!is.na(diff_auc), 
                                         .(metric_lab, analysis_lab, score_lab, diff_auc, p_delong)])

figS4_data[, sig := fcase(
  is.na(p_delong), "",
  p_delong < 0.001, "***",
  p_delong < 0.01,  "**",
  p_delong < 0.05,  "*",
  default = ""
)]

figS4 = ggplot(figS4_data, aes(x = analysis_lab, y = score_lab, fill = diff_auc)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sig), vjust = 0.5, size = 5, fontface = "bold") +
  scale_fill_gradient2(
    name     = "ΔAUROC\n(Cancer − Non-cancer)",
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = c(-0.08, 0.08)
  ) +
  facet_wrap(~ metric_lab, nrow = 2) +
  labs(
    x        = NULL, 
    y        = NULL,
    title    = "AUROC Differences Across Sensitivity Analyses",
    subtitle = "* p<0.05, ** p<0.01, *** p<0.001"
  ) +
  theme_paper +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1),
    legend.position  = "right",
    legend.key.height = unit(1.5, "cm"),
    panel.grid       = element_blank()
  )

ggsave(file.path(fig_dir, "figS4_auroc_diff_heatmap.pdf"), figS4, width = 9, height = 8)
ggsave(file.path(fig_dir, "figS4_auroc_diff_heatmap.png"), figS4, width = 9, height = 8, dpi = 300)

message("  Figure S4 saved")

# ==============================================================================
# FIGURE S5: Time to First Positivity
# ==============================================================================

message("\n== Creating Figure S5: Time to First Positivity ==")

figS5_data = first_pos_final[outcome == 1]

figS5 = ggplot(figS5_data, aes(x = score_lab, y = mean_hours / 24, fill = ca_lab)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "white", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = (mean_hours - sd_hours) / 24, 
        ymax = (mean_hours + sd_hours) / 24),
    position = position_dodge(width = 0.8),
    width    = 0.25,
    linewidth = 0.6
  ) +
  scale_fill_manual(values = colors_cancer, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x        = NULL,
    y        = "Days to First Positive (Mean ± SD)",
    title    = "Time to First Score Positivity",
    subtitle = "Among patients who deteriorated"
  ) +
  theme_paper +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(fig_dir, "figS5_time_to_positive.pdf"), figS5, width = 8, height = 5)
ggsave(file.path(fig_dir, "figS5_time_to_positive.png"), figS5, width = 8, height = 5, dpi = 300)

message("  Figure S5 saved")

# ==============================================================================
# FIGURE S6: Sensitivity Analysis Summary (Cancer Patients Only)
# ==============================================================================

message("\n== Creating Figure S6: Sensitivity Analysis Summary ==")

figS6_data = auroc_results_final[metric == "Encounter max" & ca_01 == 1]

figS6 = ggplot(figS6_data, aes(x = auroc, y = score_lab, color = analysis_lab)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.4) +
  geom_pointrange(
    aes(xmin = ci_lower, xmax = ci_upper),
    position = position_dodge(width = 0.7),
    size     = 0.5,
    fatten   = 2.5
  ) +
  scale_color_brewer(palette = "Set1", name = "Analysis") +
  scale_x_continuous(limits = c(0.55, 0.80), breaks = seq(0.55, 0.80, 0.05)) +
  labs(
    x        = "AUROC (95% CI)",
    y        = NULL,
    title    = "Sensitivity Analyses: Cancer Patients",
    subtitle = "Encounter-level maximum scores"
  ) +
  theme_paper +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position    = "right"
  )

ggsave(file.path(fig_dir, "figS6_sensitivity.pdf"), figS6, width = 10, height = 5)
ggsave(file.path(fig_dir, "figS6_sensitivity.png"), figS6, width = 10, height = 5, dpi = 300)

message("  Figure S6 saved")

# ==============================================================================
# COMBINED FIGURE: Main Results Panel (for manuscript)
# ==============================================================================

message("\n== Creating Combined Main Figure ==")

combined_main = (fig1 + labs(title = NULL)) / 
                (fig2 + labs(title = NULL, subtitle = NULL)) + 
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 0.9))

ggsave(file.path(fig_dir, "fig_main_combined.pdf"), combined_main, width = 10, height = 9)
ggsave(file.path(fig_dir, "fig_main_combined.png"), combined_main, width = 10, height = 9, dpi = 300)

message("  Combined main figure saved")

# ==============================================================================
# FIGURE S7: UpSet Plots (Score Combination Patterns)
# ==============================================================================

message("\n== Creating Figure S7: UpSet Plots ==")

library(UpSetR)

# prepare data for UpSetR (needs binary matrix format)
prepare_upset_matrix = function(upset_data, cancer_val) {
  
  agg = upset_data[ca_01 == cancer_val, 
                   .(n = sum(n)), 
                   by = .(sirs, qsofa, mews, news, mews_sf)]
  
  # expand to individual rows for UpSetR
  expanded_list = vector("list", nrow(agg))
  
  for (i in seq_len(nrow(agg))) {
    row_df = data.frame(
      SIRS    = agg$sirs[i],
      qSOFA   = agg$qsofa[i],
      MEWS    = agg$mews[i],
      NEWS    = agg$news[i],
      MEWS_SF = agg$mews_sf[i]
    )
    expanded_list[[i]] = row_df[rep(1, agg$n[i]), ]
  }
  
  do.call(rbind, expanded_list)
}

if (exists("upset_final") && nrow(upset_final) > 0) {
  
  upset_cancer    = prepare_upset_matrix(upset_final, 1)
  upset_noncancer = prepare_upset_matrix(upset_final, 0)
  
  n_cancer    = nrow(upset_cancer)
  n_noncancer = nrow(upset_noncancer)
  
  message("  Cancer encounters: ", format(n_cancer, big.mark = ","))
  message("  Non-cancer encounters: ", format(n_noncancer, big.mark = ","))
  
  # Figure S7A: Cancer patients
  pdf(file.path(fig_dir, "figS7a_upset_cancer.pdf"), width = 10, height = 6)
  print(
    upset(
      upset_cancer,
      sets           = c("SIRS", "qSOFA", "MEWS", "NEWS", "MEWS_SF"),
      order.by       = "freq",
      number.angles  = 0,
      main.bar.color = "#E74C3C",
      sets.bar.color = "#E74C3C",
      matrix.color   = "#2C3E50",
      point.size     = 3.5,
      line.size      = 1,
      text.scale     = c(1.3, 1.2, 1.2, 1.1, 1.3, 1.1),
      mainbar.y.label = "Number of Encounters",
      sets.x.label    = "Set Size"
    )
  )
  dev.off()
  
  png(file.path(fig_dir, "figS7a_upset_cancer.png"), width = 10, height = 6, units = "in", res = 300)
  print(
    upset(
      upset_cancer,
      sets           = c("SIRS", "qSOFA", "MEWS", "NEWS", "MEWS_SF"),
      order.by       = "freq",
      number.angles  = 0,
      main.bar.color = "#E74C3C",
      sets.bar.color = "#E74C3C",
      matrix.color   = "#2C3E50",
      point.size     = 3.5,
      line.size      = 1,
      text.scale     = c(1.3, 1.2, 1.2, 1.1, 1.3, 1.1),
      mainbar.y.label = "Number of Encounters",
      sets.x.label    = "Set Size"
    )
  )
  dev.off()
  
  # Figure S7B: Non-cancer patients
  pdf(file.path(fig_dir, "figS7b_upset_noncancer.pdf"), width = 10, height = 6)
  print(
    upset(
      upset_noncancer,
      sets           = c("SIRS", "qSOFA", "MEWS", "NEWS", "MEWS_SF"),
      order.by       = "freq",
      number.angles  = 0,
      main.bar.color = "#3498DB",
      sets.bar.color = "#3498DB",
      matrix.color   = "#2C3E50",
      point.size     = 3.5,
      line.size      = 1,
      text.scale     = c(1.3, 1.2, 1.2, 1.1, 1.3, 1.1),
      mainbar.y.label = "Number of Encounters",
      sets.x.label    = "Set Size"
    )
  )
  dev.off()
  
  png(file.path(fig_dir, "figS7b_upset_noncancer.png"), width = 10, height = 6, units = "in", res = 300)
  print(
    upset(
      upset_noncancer,
      sets           = c("SIRS", "qSOFA", "MEWS", "NEWS", "MEWS_SF"),
      order.by       = "freq",
      number.angles  = 0,
      main.bar.color = "#3498DB",
      sets.bar.color = "#3498DB",
      matrix.color   = "#2C3E50",
      point.size     = 3.5,
      line.size      = 1,
      text.scale     = c(1.3, 1.2, 1.2, 1.1, 1.3, 1.1),
      mainbar.y.label = "Number of Encounters",
      sets.x.label    = "Set Size"
    )
  )
  dev.off()
  
  message("  Figure S7 saved (A: cancer, B: non-cancer)")
  
} else {
  message("  WARNING: No upset data available for Figure S7")
}

# ==============================================================================
# SUMMARY
# ==============================================================================

message("\n== Figure generation complete ==")

figures_manifest = data.table(
  figure    = c("Fig 1", "Fig 2", "Fig 3", "Fig 4", "Fig 5",
                "Fig S1", "Fig S2", "Fig S3", "Fig S4", "Fig S5", "Fig S6", 
                "Fig S7A", "Fig S7B", "Combined"),
  content   = c("AUROC forest plot",
                "Sensitivity/specificity at thresholds",
                "Meta-analysis forest plot",
                "Cumulative incidence of positivity",
                "Score-outcome slope comparison",
                "Inclusion flow diagram",
                "Site-level heterogeneity",
                "Hematologic vs solid comparison",
                "AUROC difference heatmap",
                "Time to first positivity",
                "Sensitivity analyses",
                "UpSet plot (cancer)",
                "UpSet plot (non-cancer)",
                "Panels A-B for manuscript"),
  filename  = c("fig1_auroc_forest", "fig2_sens_spec", "fig3_meta_forest",
                "fig4_cumulative_incidence", "fig5_slopes",
                "figS1_flow", "figS2_site_heterogeneity", "figS3_cancer_type",
                "figS4_auroc_diff_heatmap", "figS5_time_to_positive", 
                "figS6_sensitivity", "figS7a_upset_cancer", "figS7b_upset_noncancer",
                "fig_main_combined")
)

message("\nFigures saved to: ", fig_dir)
print(figures_manifest)
