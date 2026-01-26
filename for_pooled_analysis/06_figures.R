# ==============================================================================
# 06_figures.R
# All figure generation for pooled analyses
# Sources data objects from 01_tables.R through 05_subgroups.R
# ==============================================================================

# setup ------------------------------------------------------------------------

library(ggplot2)
library(patchwork)
library(scales)
library(data.table)

# create output directory ------------------------------------------------------

if (!dir.exists(here("output", "figures"))) {
  dir.create(here("output", "figures"), recursive = TRUE)
}

# theme setup ------------------------------------------------------------------

theme_paper = theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray95"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom"
  )

# Color palettes
colors_cancer = c("Non-cancer" = "#1b9e77", "Cancer" = "#d95f02")
colors_type   = c("Solid" = "#7570b3", "Hematologic" = "#e7298a")

# ==============================================================================
# FIGURE 1: AUROC Forest Plot (Main Result)
# ==============================================================================

message("\n== Creating Figure 1: AUROC Forest Plot ==")

fig1_data = auroc_results_final[analysis == "main"]

fig1 = ggplot(fig1_data, aes(x = auroc, y = score_lab, color = ca_lab)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.3) +
  geom_pointrange(
    aes(xmin = ci_lower, xmax = ci_upper),
    position = position_dodge(width = 0.5),
    size     = 0.8
  ) +
  scale_color_manual(values = colors_cancer, name = NULL) +
  scale_x_continuous(limits = c(0.5, 0.85), breaks = seq(0.5, 0.85, 0.05)) +
  facet_wrap(~ metric_lab, ncol = 2) +
  labs(
    x = "AUROC (95% CI)",
    y = NULL
  ) +
  theme_paper +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position    = "top"
  )

ggsave(here("output", "figures", "fig1_auroc_forest.png"), fig1, 
       width = 10, height = 5, dpi = 300)
ggsave(here("output", "figures", "fig1_auroc_forest.pdf"), fig1, 
       width = 10, height = 5)

message("  Figure 1 saved")

# ==============================================================================
# FIGURE 2: AUROC Difference Heatmap
# ==============================================================================

message("\n== Creating Figure 2: AUROC Difference Heatmap ==")

fig2_data = unique(auroc_results_final[, .(metric_lab, analysis_lab, score_lab, diff_auc, p_delong)])

fig2_data[, sig := fcase(
  p_delong < 0.001, "***",
  p_delong < 0.01,  "**",
  p_delong < 0.05,  "*",
  default = ""
)]

fig2 = ggplot(fig2_data, aes(x = analysis_lab, y = score_lab, fill = diff_auc)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig), vjust = 0.5, size = 4) +
  scale_fill_gradient2(
    name     = "ΔAUROC (Cancer − Non-cancer)",
    low      = "#b2182b",
    mid      = "white",
    high     = "#2166ac",
    midpoint = 0,
    limits   = c(-0.1, 0.1)
  ) +
  facet_wrap(~ metric_lab, nrow = 2) +
  labs(x = NULL, y = NULL) +
  theme_paper +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1),
    legend.position  = "bottom",
    legend.key.width = unit(1.5, "cm")
  )

ggsave(here("output", "figures", "fig2_auroc_diff_heatmap.png"), fig2, 
       width = 8, height = 7, dpi = 300)
ggsave(here("output", "figures", "fig2_auroc_diff_heatmap.pdf"), fig2, 
       width = 8, height = 7)

message("  Figure 2 saved")

# ==============================================================================
# FIGURE 3: Sensitivity/Specificity at Standard Thresholds
# ==============================================================================

message("\n== Creating Figure 3: Sensitivity/Specificity ==")

fig3_data = sesp_final[, .(score_name, ca_01, sensitivity, specificity, positivity)]
fig3_data[, `:=`(
  ca_lab    = fifelse(ca_01 == 1, "Cancer", "Non-cancer"),
  score_lab = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
)]

# Reshape for plotting
fig3_long = melt(
  fig3_data,
  id.vars       = c("score_lab", "ca_lab"),
  measure.vars  = c("sensitivity", "specificity"),
  variable.name = "metric",
  value.name    = "value"
)
fig3_long[, metric := fifelse(metric == "sensitivity", "Sensitivity", "Specificity")]

fig3 = ggplot(fig3_long, aes(x = score_lab, y = value, fill = ca_lab)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = colors_cancer, name = NULL) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  facet_wrap(~ metric) +
  labs(x = NULL, y = NULL) +
  theme_paper +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave(here("output", "figures", "fig3_sens_spec.png"), fig3, 
       width = 8, height = 5, dpi = 300)
ggsave(here("output", "figures", "fig3_sens_spec.pdf"), fig3, 
       width = 8, height = 5)

message("  Figure 3 saved")

# ==============================================================================
# FIGURE 4: Meta-analysis Forest Plot (Interaction Terms)
# ==============================================================================

message("\n== Creating Figure 4: Meta-analysis Forest Plot ==")

fig4_data = forest_data_final[!is.na(or_int)]

fig4 = ggplot(fig4_data, aes(x = or_int, y = reorder(site, -is_pooled), 
                              color = is_pooled, size = is_pooled)) +
  geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = c(0.98, 1.02), linetype = "dotted", alpha = 0.3, color = "red") +
  geom_pointrange(aes(xmin = or_int_lower, xmax = or_int_upper)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray50"), guide = "none") +
  scale_size_manual(values = c("TRUE" = 1.2, "FALSE" = 0.6), guide = "none") +
  scale_x_continuous(trans = "log", breaks = c(0.9, 0.95, 1, 1.05, 1.1)) +
  facet_wrap(~ score_lab, ncol = 3) +
  labs(
    x = "Interaction OR per score point (95% CI)",
    y = NULL
  ) +
  theme_paper

ggsave(here("output", "figures", "fig4_meta_forest.png"), fig4, 
       width = 10, height = 8, dpi = 300)
ggsave(here("output", "figures", "fig4_meta_forest.pdf"), fig4, 
       width = 10, height = 8)

message("  Figure 4 saved")

# ==============================================================================
# FIGURE 5: Slope Comparison (Score-Outcome Relationship)
# ==============================================================================

message("\n== Creating Figure 5: Slope Comparison ==")

if (nrow(slope_preds_final) > 0) {
  
  fig5 = ggplot(slope_preds_final, aes(x = score_value, y = predicted_prob,
                                        color = cancer_status, linetype = cancer_status)) +
    geom_line(linewidth = 1) +
    scale_y_continuous(labels = percent_format(), limits = c(0, NA)) +
    scale_color_manual(values = c("Non-cancer" = "#2166ac", "Cancer" = "#d95f02")) +
    scale_linetype_manual(values = c("Non-cancer" = "solid", "Cancer" = "dashed")) +
    facet_wrap(~ score_label, scales = "free_x", ncol = 3) +
    labs(
      x        = "Score Value",
      y        = "Predicted Probability of Deterioration",
      color    = NULL,
      linetype = NULL
    ) +
    theme_paper +
    theme(
      legend.position = "bottom",
      strip.text      = element_text(size = 9)
    )
  
  ggsave(here("output", "figures", "fig5_slopes.png"), fig5, 
         width = 10, height = 7, dpi = 300)
  ggsave(here("output", "figures", "fig5_slopes.pdf"), fig5, 
         width = 10, height = 7)
  
  message("  Figure 5 saved")
  
} else {
  message("  WARNING: No slope prediction data available for Figure 5")
}

# ==============================================================================
# FIGURE S1: Inclusion Flow Diagram
# ==============================================================================

message("\n== Creating Figure S1: Flow Diagram ==")

library(DiagrammeR)

# Build flow diagram from flow_data_final (created in 01_tables.R)
ca_start = flow_data_final[1]$n_remaining_ca
no_start = flow_data_final[1]$n_remaining_no
ca_final = flow_data_final[nrow(flow_data_final)]$n_remaining_ca
no_final = flow_data_final[nrow(flow_data_final)]$n_remaining_no

diagram_code = paste0(
  "digraph flowchart {
    graph [rankdir = TB, splines = ortho]
    node [shape = box, fontname = Helvetica, fontsize = 10, width = 3]
    
    A [label = 'Adult inpatient admissions during study period\\nPatients with Cancer: ", format_n(ca_start),
  "\\nPatients without Cancer: ", format_n(no_start), "']"
)

step_labels = flow_data_final$step

for (i in 2:nrow(flow_data_final)) {
  current_node = LETTERS[i]
  prev_node    = LETTERS[i-1]
  excl_node    = paste0("E", i-1)
  
  excl_label = gsub("After excluding ", "", step_labels[i])
  excl_label = gsub("patients ", "", excl_label)
  
  diagram_code = paste0(
    diagram_code, 
    "\n  ", excl_node,
    " [label = 'Excluded: ", excl_label,
    "\\nWith Cancer: ",    format_n(flow_data_final$n_excluded_ca[i]),
    "\\nWithout Cancer: ", format_n(flow_data_final$n_excluded_no[i]), 
    "', style = filled, fillcolor = '#ffcccc']",
    "\n  ", current_node,
    " [label = 'Remaining\\nWith Cancer: ",    format_n(flow_data_final$n_remaining_ca[i]),
    "\\nWithout Cancer: ", format_n(flow_data_final$n_remaining_no[i]), "']",
    "\n  ", prev_node, " -> ", excl_node, " [style = dashed]",
    "\n  ", prev_node, " -> ", current_node
  )
}

# Add final analysis cohort box
final_node = LETTERS[nrow(flow_data_final) + 1]
diagram_code = paste0(
  diagram_code,
  "\n  ", final_node, " [label = 'Final Analysis Cohort\\nWith Cancer: ", format_n(ca_final),
  "\\nWithout Cancer: ", format_n(no_final), "', style = filled, fillcolor = '#ccffcc']",
  "\n  ", LETTERS[nrow(flow_data_final)], " -> ", final_node,
  "\n}"
)

flow_diagram_s1 = grViz(diagram_code)

# Export to PNG using DiagrammeRsvg and rsvg if available
tryCatch({
  library(DiagrammeRsvg)
  library(rsvg)
  
  svg_text = export_svg(flow_diagram_s1)
  rsvg_png(charToRaw(svg_text), here("output", "figures", "figS1_flow.png"), width = 1200)
  message("  Figure S1 saved")
  
}, error = function(e) {
  message("  WARNING: Could not export flow diagram as PNG. Install DiagrammeRsvg and rsvg packages.")
  message("  Flow diagram object saved but not rendered to file.")
})

# ==============================================================================
# FIGURE S2: Site-Level AUROC Heterogeneity
# ==============================================================================

message("\n== Creating Figure S2: Site-Level Heterogeneity ==")

figS2 = ggplot(site_aurocs_final, 
               aes(x = site_ordered, y = auroc, color = ca_lab)) +
  geom_pointrange(
    aes(ymin = ci_lower, ymax = ci_upper),
    position = position_dodge(width = 0.6),
    size     = 0.5
  ) +
  geom_hline(
    data = site_heterogeneity_final,
    aes(yintercept = mean_auroc, color = ca_lab),
    linetype = "dashed",
    alpha    = 0.5
  ) +
  scale_color_manual(values = colors_cancer, name = NULL) +
  scale_y_continuous(limits = c(0.4, 0.9), breaks = seq(0.4, 0.9, 0.1)) +
  facet_wrap(~ score_lab, ncol = 3) +
  labs(
    x = "Site (ordered by mean rank across scores)",
    y = "AUROC (95% CI)"
  ) +
  theme_paper +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "top"
  )

ggsave(here("output", "figures", "figS2_site_heterogeneity.png"), figS2, 
       width = 12, height = 8, dpi = 300)
ggsave(here("output", "figures", "figS2_site_heterogeneity.pdf"), figS2, 
       width = 12, height = 8)

message("  Figure S2 saved")

# ==============================================================================
# FIGURE S3: Hematologic vs Solid Tumor Comparison
# ==============================================================================

message("\n== Creating Figure S3: Cancer Type Comparison ==")

if (nrow(liquid_aurocs_final) > 0) {
  
  figS3 = ggplot(liquid_aurocs_final, aes(x = auroc, y = score_lab, color = cancer_type)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.3) +
    geom_pointrange(
      aes(xmin = ci_lower, xmax = ci_upper),
      position = position_dodge(width = 0.5),
      size     = 0.8
    ) +
    scale_color_manual(values = colors_type, name = "Cancer Type") +
    scale_x_continuous(limits = c(0.5, 0.85), breaks = seq(0.5, 0.85, 0.05)) +
    labs(
      x = "AUROC (95% CI)",
      y = NULL,
      title = "Score Performance: Hematologic vs Solid Tumors"
    ) +
    theme_paper +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position    = "top"
    )
  
  ggsave(here("output", "figures", "figS3_cancer_type.png"), figS3, 
         width = 8, height = 5, dpi = 300)
  ggsave(here("output", "figures", "figS3_cancer_type.pdf"), figS3, 
         width = 8, height = 5)
  
  message("  Figure S3 saved")
  
} else {
  message("  WARNING: No liquid/solid data available for Figure S3")
}

# ==============================================================================
# FIGURE S4: Cumulative Incidence of Score Positivity
# ==============================================================================

message("\n== Creating Figure S4: Cumulative Incidence ==")

figS4 = ggplot(cuminc_final, 
               aes(x = time_bin_start / 24, y = cum_inc, 
                   color = ca_lab, linetype = outcome_lab)) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  scale_color_manual(values = colors_cancer, name = NULL) +
  scale_linetype_manual(values = c("With deterioration" = "solid", 
                                    "Without deterioration" = "dashed"),
                        name = NULL) +
  facet_wrap(~ score_lab, ncol = 3) +
  labs(
    x = "Days from Ward Admission",
    y = "Cumulative Proportion Ever Positive"
  ) +
  theme_paper +
  theme(legend.position = "bottom")

ggsave(here("output", "figures", "figS4_cumulative_incidence.png"), figS4, 
       width = 10, height = 7, dpi = 300)
ggsave(here("output", "figures", "figS4_cumulative_incidence.pdf"), figS4, 
       width = 10, height = 7)

message("  Figure S4 saved")

# ==============================================================================
# FIGURE S5: Time to First Positivity Distribution
# ==============================================================================

message("\n== Creating Figure S5: Time to First Positivity ==")

figS5_data = first_pos_final[outcome == 1]

figS5 = ggplot(figS5_data, aes(x = score_lab, y = mean_hours / 24, fill = ca_lab)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = (mean_hours - sd_hours) / 24, 
        ymax = (mean_hours + sd_hours) / 24),
    position = position_dodge(width = 0.8),
    width    = 0.2
  ) +
  scale_fill_manual(values = colors_cancer, name = NULL) +
  labs(
    x = NULL,
    y = "Days to First Positive Score (Mean ± SD)",
    title = "Time to First Score Positivity Among Patients Who Deteriorated"
  ) +
  theme_paper +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave(here("output", "figures", "figS5_time_to_positive.png"), figS5, 
       width = 8, height = 5, dpi = 300)
ggsave(here("output", "figures", "figS5_time_to_positive.pdf"), figS5, 
       width = 8, height = 5)

message("  Figure S5 saved")

# ==============================================================================
# FIGURE S6: Sensitivity Analysis Summary
# ==============================================================================

message("\n== Creating Figure S6: Sensitivity Analysis Summary ==")

# Filter to encounter-level max and cancer patients only
figS6_data = auroc_results_final[metric == "Encounter max" & ca_01 == 1]

figS6 = ggplot(figS6_data, aes(x = auroc, y = score_lab, color = analysis_lab)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.3) +
  geom_pointrange(
    aes(xmin = ci_lower, xmax = ci_upper),
    position = position_dodge(width = 0.7),
    size     = 0.6
  ) +
  scale_color_brewer(palette = "Set1", name = "Analysis") +
  scale_x_continuous(limits = c(0.5, 0.85), breaks = seq(0.5, 0.85, 0.05)) +
  labs(
    x     = "AUROC (95% CI)",
    y     = NULL,
    title = "Sensitivity Analyses: Cancer Patients"
  ) +
  theme_paper +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position    = "right"
  )

ggsave(here("output", "figures", "figS6_sensitivity.png"), figS6, 
       width = 10, height = 5, dpi = 300)
ggsave(here("output", "figures", "figS6_sensitivity.pdf"), figS6, 
       width = 10, height = 5)

message("  Figure S6 saved")

# ==============================================================================
# COMBINED FIGURE: Main Results Panel (for manuscript)
# ==============================================================================

message("\n== Creating Combined Main Figure ==")

# Combine Figures 1 and 3 for a comprehensive main figure
combined_main = fig1 / fig3 + 
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 0.8))

ggsave(here("output", "figures", "fig_main_combined.png"), combined_main, 
       width = 10, height = 9, dpi = 300)
ggsave(here("output", "figures", "fig_main_combined.pdf"), combined_main, 
       width = 10, height = 9)

message("  Combined main figure saved")

# ==============================================================================
# EXPORTS LIST
# ==============================================================================

message("\n== Figure generation complete ==")

figures_manifest = data.table(
  figure    = c("Fig 1", "Fig 2", "Fig 3", "Fig 4", "Fig 5",
                "Fig S1", "Fig S2", "Fig S3", "Fig S4", "Fig S5", "Fig S6",
                "Combined Main"),
  content   = c("AUROC forest plot (main analysis)",
                "AUROC difference heatmap (all analyses)",
                "Sensitivity/specificity at thresholds",
                "Meta-analysis forest plot (interaction terms)",
                "Score-outcome slope comparison",
                "Inclusion flow diagram",
                "Site-level AUROC heterogeneity",
                "Hematologic vs solid tumor comparison",
                "Cumulative incidence of positivity",
                "Time to first positivity",
                "Sensitivity analyses summary",
                "Panels A-B for manuscript"),
  file_base = c("fig1_auroc_forest", "fig2_auroc_diff_heatmap", "fig3_sens_spec",
                "fig4_meta_forest", "fig5_slopes",
                "figS1_flow", "figS2_site_heterogeneity", "figS3_cancer_type",
                "figS4_cumulative_incidence", "figS5_time_to_positive", 
                "figS6_sensitivity", "fig_main_combined")
)

message("  Figures saved to: ", here("output", "figures"))
print(figures_manifest)
