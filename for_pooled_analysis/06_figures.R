# ==============================================================================
# figures_clean.R
# Publication-quality figures for pooled EWS analyses
# ==============================================================================


# ==============================================================================
# SETUP
# ==============================================================================

library(patchwork)
library(ggplot2)
library(scales)
library(purrr)

## output directory ------------------------------------------------------------

if (!dir.exists(here("output", "figures"))) {
  dir.create(here("output", "figures"), recursive = TRUE)
}

## shared theme and palette ----------------------------------------------------

pal_cancer = c("Non-Cancer" = "#4575b4", "Cancer" = "#d73027")

theme_ews = function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      strip.text       = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "gray85", color = "black"),
      legend.position  = "none"
    )
}

## score labels ----------------------------------------------------------------

score_levels = c("sirs", "qsofa", "mews", "mews_sf", "news")
score_labels = c("SIRS", "QSOFA", "MEWS", "MEWS-SF", "NEWS")
names(score_labels) = score_levels

## helper functions ------------------------------------------------------------

format_score = function(x) {
  x = tolower(x)
  x = fifelse(x == "mews-sf", "mews_sf", x)
  factor(x, levels = score_levels, labels = score_labels)
}

format_cohort = function(ca_01, n_by_cohort = NULL) {
  # If n_by_cohort provided, create labels with sample sizes
  if (!is.null(n_by_cohort)) {
    lab_noca = sprintf("Non-Cancer (n=%s)", format(n_by_cohort["0"], big.mark = ","))
    lab_ca   = sprintf("Cancer (n=%s)", format(n_by_cohort["1"], big.mark = ","))
    factor(
      fifelse(ca_01 == 1, lab_ca, lab_noca),
      levels = c(lab_noca, lab_ca)
    )
  } else {
    factor(
      fifelse(ca_01 == 1, "Cancer", "Non-Cancer"),
      levels = c("Non-Cancer", "Cancer")
    )
  }
}

# Helper to build palette matching cohort labels (with or without n)
build_cohort_palette = function(cohort_labels) {
  setNames(c("#4575b4", "#d73027"), cohort_labels)
}

# Verify constants from 00_load.R
if (!exists("COHORT_N")) stop("COHORT_N not found. Source 00_load.R first.")
if (!exists("VARIANT_N")) message("WARNING: VARIANT_N not found")
if (!exists("SITE_N")) message("WARNING: SITE_N not found")
message("Using COHORT_N: ", format(COHORT_N["0"], big.mark=","), " non-cancer, ",
        format(COHORT_N["1"], big.mark=","), " cancer")

calc_auroc_from_counts = function(data) {
  
  safe_fsum = function(x) {
    result = fsum(as.numeric(x))
    if (length(result) == 0) return(0)
    result
  }
  
  events     = fsubset(data, outcome == 1)
  non_events = fsubset(data, outcome == 0)
  
  n_events     = safe_fsum(events$n)
  n_non_events = safe_fsum(non_events$n)
  
  concordant = 0
  tied       = 0
  
  for (i in seq_len(nrow(events))) {
    ev_score   = events$value[i]
    ev_n       = as.numeric(events$n[i])
    lower_idx  = which(non_events$value < ev_score)
    concordant = concordant + ev_n * safe_fsum(non_events$n[lower_idx])
    tied_idx   = which(non_events$value == ev_score)
    tied       = tied + ev_n * safe_fsum(non_events$n[tied_idx])
  }
  
  auroc = (concordant + 0.5 * tied) / (n_events * n_non_events)
  q1    = auroc / (2 - auroc)
  q2    = 2 * auroc^2 / (1 + auroc)
  se    = sqrt((auroc * (1 - auroc) + (n_events - 1) * (q1 - auroc^2) + 
                  (n_non_events - 1) * (q2 - auroc^2)) / (n_events * n_non_events))
  
  tidytable(
    auroc    = auroc,
    se       = se,
    ci_lower = auroc - 1.96 * se,
    ci_upper = auroc + 1.96 * se,
    n_events = as.integer(n_events),
    n_total  = as.integer(n_events + n_non_events)
  )
}


# ==============================================================================
# FIGURE 1: Risk by Score Value
# ==============================================================================

## data prep -------------------------------------------------------------------

fig1_data = maxscores_ca_raw |>
  fsubset(tolower(analysis) == "main") |>
  fmutate(n_events = n * outcome) |>
  fgroup_by(score_name, ca_01, max_value) |>
  fsummarise(n_at_score = fsum(n), n_outcomes = fsum(n_events)) |>
  fungroup() |>
  fmutate(
    prob         = n_outcomes / n_at_score,
    score_label  = format_score(score_name),
    cohort_label = format_cohort(ca_01, COHORT_N)
  ) |>
  fsubset(n_outcomes >= 100)

fig1_pal = build_cohort_palette(levels(fig1_data$cohort_label))

## panel function --------------------------------------------------------------

fig1_panel = function(data, score, show_y = TRUE, x_breaks = NULL, pal = fig1_pal) {
  
  panel_data = fsubset(data, score_label == score)
  
  p = ggplot(panel_data, aes(x = max_value, y = prob, color = cohort_label)) +
    geom_line(linewidth = 0.6, alpha = 0.3) +
    geom_point(aes(size = n_at_score), alpha = 0.75) +
    facet_wrap(~score_label) +
    scale_color_manual(values = pal) +
    scale_size_area(labels = label_comma(), breaks = c(100000, 300000, 500000)) +
    scale_y_continuous(labels = label_percent(), limits = c(0, 0.9), breaks = seq(0, 0.9, 0.15)) +
    theme_ews() +
    theme(panel.grid.major.x = element_blank()) +
    labs(x = NULL)
  
  if (!is.null(x_breaks)) p = p + scale_x_continuous(breaks = x_breaks)
  if (show_y) p = p + labs(y = "Observed Deterioration Rate (%)") else p = p + labs(y = NULL)
  
  p
}

## build panels ----------------------------------------------------------------

fig1_sirs   = fig1_panel(fig1_data, "SIRS",    show_y = TRUE,  x_breaks = 0:4)
fig1_qsofa  = fig1_panel(fig1_data, "QSOFA",   show_y = FALSE, x_breaks = 0:3)
fig1_mews   = fig1_panel(fig1_data, "MEWS",    show_y = FALSE, x_breaks = seq(0, 10, 2))
fig1_mewssf = fig1_panel(fig1_data, "MEWS-SF", show_y = TRUE,  x_breaks = seq(0, 10, 2))
fig1_news   = fig1_panel(fig1_data, "NEWS",    show_y = FALSE, x_breaks = seq(0, 10, 2))

## legend ----------------------------------------------------------------------

fig1_legend_data = tidytable(
  cohort_label = factor(levels(fig1_data$cohort_label), levels = levels(fig1_data$cohort_label)),
  n_at_score   = c(100000, 500000),
  x = 1, y = 1
)

fig1_legend = ggplot(fig1_legend_data, aes(x = x, y = y, color = cohort_label, size = n_at_score)) +
  geom_point() +
  scale_color_manual(values = fig1_pal, name = "Cohort") +
  scale_size_area(labels = label_comma(), breaks = c(100000, 300000, 500000)) +
  xlim(10, 20) + ylim(10, 20) +
  guides(
    color = guide_legend(title = "Cohort", order = 1),
    size  = guide_legend(title = "Total Patients", order = 2)
  ) +
  theme_void() +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.5, 0.5),
    legend.box             = "vertical",
    legend.spacing.y       = unit(0.8, "cm")
  )

## assemble and save -----------------------------------------------------------

fig1 = (fig1_sirs | fig1_qsofa | fig1_mews) /
  (fig1_mewssf | fig1_news | fig1_legend) +
  plot_annotation(
    caption = "Maximum Score Value Within Encounter",
    theme   = theme(plot.caption = element_text(hjust = 0.5, size = 11))
  )

fig1

ggsave(here("output", "figures", "figure_01_risk_by_score.pdf"), fig1, width = 10, height = 6)

# ==============================================================================
# FIGURE 2: AUROC Comparison (Main Analysis)
# ==============================================================================

## data prep -------------------------------------------------------------------

fig2_data = auroc_results_final |>
  fsubset(analysis == "main" & metric == "Encounter max") |>
  fmutate(
    score_label  = format_score(score_name),
    cohort_label = format_cohort(ca_01, COHORT_N)
  )

fig2_pal = build_cohort_palette(levels(fig2_data$cohort_label))

## plot and save ---------------------------------------------------------------

fig2 = ggplot(fig2_data, aes(x = score_label, y = auroc, fill = cohort_label)) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.2, linewidth = 0.5, color = "black",
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    shape = 21, color = "black", size = 4, stroke = 0.7,
    position = position_dodge(width = 0.6)
  ) +
  scale_fill_manual(values = fig2_pal, name = "Cohort") +
  scale_y_continuous(limits = c(0.58, 0.80), breaks = seq(0.60, 0.80, 0.05)) +
  labs(x = NULL, y = "AUROC (95% CI)") +
  theme_ews() +
  theme(legend.position = "top", axis.text.x = element_text(face = "bold", size = 11))

fig2

ggsave(here("output", "figures", "figure_02_auroc_main.pdf"), fig2, width = 7, height = 5)

# ==============================================================================
# FIGURE 3: Cumulative Incidence of Positivity
# ==============================================================================

## data prep -------------------------------------------------------------------

fig3_data = cuminc_raw |>
  fsubset(analysis == "main" & o_primary_01 == 1) |>
  fgroup_by(score, ca_01, time_bin_start) |>
  fsummarise(n_at_risk = fsum(n_at_risk), n_became_pos = fsum(n_became_pos)) |>
  fungroup() |>
  fmutate(
    cum_inc      = n_became_pos / n_at_risk,
    score_label  = format_score(score),
    cohort_label = format_cohort(ca_01, COHORT_N),
    time_days    = time_bin_start / 24
  )

fig3_pal = build_cohort_palette(levels(fig3_data$cohort_label))

## plot and save ---------------------------------------------------------------

fig3 = ggplot(fig3_data, aes(x = time_days, y = cum_inc, color = cohort_label)) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~score_label, nrow = 1) +
  scale_color_manual(values = fig3_pal, name = "Cohort") +
  scale_x_continuous(breaks = seq(0, 7, 2)) +
  scale_y_continuous(labels = label_percent()) +
  labs(x = "Days from Ward Admission", y = "Cumulative Incidence of Positivity") +
  theme_ews() +
  theme(legend.position = "top")

fig3

ggsave(here("output", "figures", "figure_03_cuminc.pdf"), fig3, width = 10, height = 4)

# ==============================================================================
# FIGURE 4: Threshold-Performance Plots (SIRS and NEWS)
# ==============================================================================

## data prep -------------------------------------------------------------------

fig4_prep = maxscores_ca_raw |>
  fsubset(tolower(analysis) == "main") |>
  fgroup_by(score_name, ca_01, max_value, outcome) |>
  fsummarise(n = fsum(n)) |>
  fungroup()

calc_threshold_metrics = function(data, score, cohort) {
  
  d = fsubset(data, score_name == score & ca_01 == cohort)
  
  totals    = d |> fgroup_by(outcome) |> fsummarise(n = fsum(n)) |> fungroup()
  total_pos = totals$n[totals$outcome == 1]
  total_neg = totals$n[totals$outcome == 0]
  if (length(total_pos) == 0) total_pos = 0
  if (length(total_neg) == 0) total_neg = 0
  
  thresholds = sort(unique(d$max_value))
  
  safe_sum = function(x) {
    if (length(x) == 0) return(0)
    s = fsum(x)
    if (length(s) == 0 || is.na(s)) return(0)
    s
  }
  
  results = lapply(thresholds, function(thresh) {
    above = fsubset(d, max_value >= thresh)
    below = fsubset(d, max_value < thresh)
    tp = safe_sum(above$n[above$outcome == 1])
    fp = safe_sum(above$n[above$outcome == 0])
    fn = safe_sum(below$n[below$outcome == 1])
    tn = safe_sum(below$n[below$outcome == 0])
    tidytable(
      score_name = score, ca_01 = cohort, threshold = thresh,
      sens = ifelse(tp + fn > 0, tp / (tp + fn), NA),
      spec = ifelse(tn + fp > 0, tn / (tn + fp), NA),
      ppv  = ifelse(tp + fp > 0, tp / (tp + fp), NA),
      npv  = ifelse(tn + fn > 0, tn / (tn + fn), NA)
    )
  })
  
  bind_rows(results)
}

fig4_combos = expand.grid(
  score  = unique(fig4_prep$score_name),
  cohort = unique(fig4_prep$ca_01),
  stringsAsFactors = FALSE
)

fig4_metrics = map2(fig4_combos$score, fig4_combos$cohort,
                    ~calc_threshold_metrics(fig4_prep, .x, .y)) |> 
  bind_rows()

fig4_data = fig4_metrics |>
  pivot_longer(cols = c(sens, spec, ppv, npv), names_to = "metric", values_to = "value") |>
  fmutate(
    value        = fifelse(metric == "npv" & is.na(value) & threshold == 0, 1, value),
    score_label  = format_score(score_name),
    cohort_label = format_cohort(ca_01, COHORT_N),
    metric_label = factor(metric,
                          levels = c("sens", "spec", "ppv", "npv"),
                          labels = c("Sensitivity", "Specificity", "PPV", "NPV"))
  )

fig4_pal = build_cohort_palette(levels(fig4_data$cohort_label))

fig4_data_main = fsubset(fig4_data, score_label %in% c("SIRS", "NEWS"))

## panel function --------------------------------------------------------------

fig4_panel = function(data, score, show_y = TRUE, pal = fig4_pal) {
  
  panel_data = fsubset(data, score_label == score)
  
  p = ggplot(panel_data, aes(x = threshold, y = value, color = cohort_label)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.5) +
    facet_wrap(~metric_label, ncol = 1) +
    scale_color_manual(values = pal, name = "Cohort") +
    scale_y_continuous(labels = label_percent(), limits = c(0, 1)) +
    labs(x = NULL, title = score) +
    theme_ews() +
    theme(plot.title = element_text(hjust = 0, face = "bold", size = 12))
  
  if (!show_y) p = p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  p
}

## build panels ----------------------------------------------------------------

fig4_sirs = fig4_panel(fig4_data_main, "SIRS", show_y = TRUE)
fig4_news = fig4_panel(fig4_data_main, "NEWS", show_y = FALSE)

## legend ----------------------------------------------------------------------

fig4_legend = ggplot(
  tidytable(cohort_label = factor(levels(fig4_data$cohort_label), levels = levels(fig4_data$cohort_label)),
            x = c(1, 2), y = c(1, 1)), 
  aes(x = x, y = y, color = cohort_label)
) +
  geom_point(size = 3) +
  scale_color_manual(values = fig4_pal, name = "Cohort") +
  xlim(10, 20) + ylim(10, 20) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_void() +
  theme(legend.position = "inside", legend.position.inside = c(0.5, 0.5), legend.direction = "horizontal")

## assemble and save -----------------------------------------------------------

fig4 = fig4_legend / (fig4_sirs | fig4_news) +
  plot_layout(heights = c(0.05, 1)) +
  plot_annotation(
    caption = "Threshold (Positive if Score ≥ Value)",
    theme   = theme(plot.caption = element_text(hjust = 0.5, size = 11))
  )

fig4

ggsave(here("output", "figures", "figure_04_threshold_performance.pdf"), fig4, width = 8, height = 9)


# ==============================================================================
# SUPP FIGURE 2: Score Distribution Histograms (Mirrored)
# ==============================================================================

## data prep -------------------------------------------------------------------

sf2_data = maxscores_ca_raw |>
  fsubset(tolower(analysis) == "main") |>
  fgroup_by(score_name, ca_01, max_value) |>
  fsummarise(n = fsum(n)) |>
  fmutate(
    score_label  = format_score(score_name),
    cohort_label = format_cohort(ca_01, COHORT_N),
    n_mirror     = fifelse(ca_01 == 1, -n, n)
  ) |>
  fsubset(max_value <= 11)

## plot and save ---------------------------------------------------------------

sf2 = ggplot(sf2_data, aes(x = n_mirror, y = factor(max_value), fill = cohort_label)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "black") +
  facet_wrap(~score_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = build_cohort_palette(levels(sf2_data$cohort_label)), name = "Cohort") +
  scale_x_continuous(labels = function(x) label_comma()(abs(x))) +
  labs(y = "Maximum Score Value", x = "Number of Encounters") +
  theme_ews() +
  theme(legend.position = "top")

sf2

ggsave(here("output", "figures", "figure_s02_score_histograms.pdf"), sf2, width = 6, height = 11)


# ==============================================================================
# SUPP FIGURE 3: AUROC Comparison (24h Analysis)
# ==============================================================================

## data prep -------------------------------------------------------------------

sf3_pooled = counts_h24_raw |>
  fsubset(analysis == "main") |>
  fgroup_by(score_name, ca_01, value, outcome) |>
  fsummarise(n = fsum(n)) |>
  fungroup()

sf3_combos = sf3_pooled |> fselect(score_name, ca_01) |> funique()

sf3_aurocs = map2(sf3_combos$score_name, sf3_combos$ca_01, function(sc, ca) {
  d = fsubset(sf3_pooled, score_name == sc & ca_01 == ca)
  result = calc_auroc_from_counts(d)
  result$score_name = sc
  result$ca_01 = ca
  result
}) |>
  bind_rows() |>
  fmutate(score_label = format_score(score_name), cohort_label = format_cohort(ca_01, COHORT_N))

sf3_pal = build_cohort_palette(levels(sf3_aurocs$cohort_label))

## plot and save ---------------------------------------------------------------

sf3 = ggplot(sf3_aurocs, aes(x = score_label, y = auroc, fill = cohort_label)) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.2, linewidth = 0.5, color = "black",
    position = position_dodge(width = 0.4)
  ) +
  geom_point(
    shape = 21, color = "black", size = 4, stroke = 0.7,
    position = position_dodge(width = 0.4)
  ) +
  scale_fill_manual(values = sf3_pal, name = "Cohort") +
  scale_y_continuous(limits = c(0.62, 0.72), breaks = seq(0.62, 0.72, 0.02)) +
  labs(x = NULL, y = "AUROC (95% CI)") +
  theme_ews() +
  theme(legend.position = "top", axis.text.x = element_text(face = "bold", size = 11))

sf3

ggsave(here("output", "figures", "figure_s03_auroc_24h.pdf"), sf3, width = 7, height = 5)


# ==============================================================================
# SUPP FIGURE 4: Site-level AUROC Caterpillar
# ==============================================================================

## data prep -------------------------------------------------------------------

sf4_prep = maxscores_ca_raw |>
  fsubset(analysis == "main") |>
  fmutate(value = max_value) |>
  fgroup_by(site, score_name, ca_01, value, outcome) |>
  fsummarise(n = fsum(n)) |>
  fungroup()

sf4_combos = sf4_prep |> fselect(site, score_name, ca_01) |> funique()

sf4_aurocs = pmap(list(sf4_combos$site, sf4_combos$score_name, sf4_combos$ca_01),
                  function(st, sc, ca) {
                    d = fsubset(sf4_prep, site == st & score_name == sc & ca_01 == ca)
                    result = calc_auroc_from_counts(d)
                    result$site = st
                    result$score_name = sc
                    result$ca_01 = ca
                    result
                  }) |> 
  bind_rows()

site_order  = sf4_aurocs |> fsubset(score_name == "news" & ca_01 == 0) |> arrange(auroc) |> pull(site)

# Get site Ns from SITE_N constant
if (exists("SITE_N")) {
  site_labels = setNames(
    sprintf("%s (n=%s)", LETTERS[seq_along(site_order)], format(SITE_N[site_order], big.mark = ",")),
    site_order
  )
  pooled_label = sprintf("Pooled (n=%s)", format(sum(SITE_N), big.mark = ","))
} else {
  site_labels = setNames(LETTERS[seq_along(site_order)], site_order)
  pooled_label = "Pooled"
}

sf4_pooled_prep   = sf4_prep |> fgroup_by(score_name, ca_01, value, outcome) |> fsummarise(n = fsum(n)) |> fungroup()
sf4_pooled_combos = sf4_pooled_prep |> fselect(score_name, ca_01) |> funique()

sf4_pooled = map2(sf4_pooled_combos$score_name, sf4_pooled_combos$ca_01, function(sc, ca) {
  d = fsubset(sf4_pooled_prep, score_name == sc & ca_01 == ca)
  result = calc_auroc_from_counts(d)
  result$site = "pooled"
  result$score_name = sc
  result$ca_01 = ca
  result
}) |> bind_rows()

sf4_all = bind_rows(sf4_aurocs, sf4_pooled) |>
  fmutate(
    score_label  = format_score(score_name),
    cohort_label = format_cohort(ca_01, COHORT_N),
    site_label   = fifelse(site == "pooled", pooled_label, site_labels[site]),
    site_label   = factor(site_label, levels = c(pooled_label, site_labels[site_order])),
    is_pooled    = site == "pooled"
  )

sf4_pal = build_cohort_palette(levels(sf4_all$cohort_label))

## plot and save ---------------------------------------------------------------

sf4 = ggplot(sf4_all, aes(x = site_label, y = auroc, color = cohort_label)) +
  geom_hline(
    data = sf4_all |> fsubset(is_pooled),
    aes(yintercept = auroc, color = cohort_label),
    linetype = "dashed", linewidth = 0.4, alpha = 0.5
  ) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.3, linewidth = 0.5,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(aes(shape = is_pooled), size = 2.5, position = position_dodge(width = 0.6)) +
  facet_wrap(~score_label, nrow = 1) +
  scale_color_manual(values = sf4_pal, name = "Cohort") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
  scale_y_continuous(limits = c(0.60, 0.82), breaks = seq(0.60, 0.80, 0.05)) +
  coord_flip() +
  labs(x = NULL, y = "AUROC (95% CI)") +
  theme_ews() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

sf4

ggsave(here("output", "figures", "figure_s04_auroc_by_site.pdf"), sf4, width = 14, height = 6)


# ==============================================================================
# SUPP FIGURE 5: Threshold-Performance Plots (QSOFA, MEWS, MEWS-SF)
# ==============================================================================

## data prep -------------------------------------------------------------------

sf5_data = fsubset(fig4_data, score_label %in% c("QSOFA", "MEWS", "MEWS-SF"))

## build panels ----------------------------------------------------------------

sf5_qsofa  = fig4_panel(sf5_data, "QSOFA",   show_y = TRUE)
sf5_mews   = fig4_panel(sf5_data, "MEWS",    show_y = FALSE)
sf5_mewssf = fig4_panel(sf5_data, "MEWS-SF", show_y = FALSE)

## legend ----------------------------------------------------------------------

sf5_legend = ggplot(
  tidytable(cohort_label = factor(levels(fig4_data$cohort_label), levels = levels(fig4_data$cohort_label)),
            x = c(1, 2), y = c(1, 1)), 
  aes(x = x, y = y, color = cohort_label)
) +
  geom_point(size = 3) +
  scale_color_manual(values = fig4_pal, name = "Cohort") +
  xlim(10, 20) + ylim(10, 20) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_void() +
  theme(legend.position = "inside", legend.position.inside = c(0.5, 0.5), legend.direction = "horizontal")

## assemble and save -----------------------------------------------------------

sf5 = sf5_legend / (sf5_qsofa | sf5_mews | sf5_mewssf) +
  plot_layout(heights = c(0.05, 1)) +
  plot_annotation(
    caption = "Threshold (Positive if Score ≥ Value)",
    theme   = theme(plot.caption = element_text(hjust = 0.5, size = 11))
  )

sf5

ggsave(here("output", "figures", "figure_s05_threshold_performance_other.pdf"), sf5, width = 11, height = 9)


# ==============================================================================
# SUPP FIGURE 6: Interaction Forest Plot
# ==============================================================================

## data prep -------------------------------------------------------------------

site_order_sf6 = sort(unique(forest_data_final$site[forest_data_final$site != "Pooled"]))

# Add site Ns to labels
if (exists("SITE_N")) {
  site_labels_sf6 = setNames(
    sprintf("%s (n=%s)", LETTERS[seq_along(site_order_sf6)], format(SITE_N[site_order_sf6], big.mark = ",")),
    site_order_sf6
  )
  pooled_label_sf6 = sprintf("Pooled (n=%s)", format(sum(SITE_N), big.mark = ","))
} else {
  site_labels_sf6 = setNames(LETTERS[seq_along(site_order_sf6)], site_order_sf6)
  pooled_label_sf6 = "Pooled"
}

sf6_data = forest_data_final |>
  fsubset(analysis == "main" | is.na(analysis)) |>
  fmutate(
    site_label  = fifelse(is_pooled, pooled_label_sf6, site_labels_sf6[site]),
    site_label  = factor(site_label, levels = c(site_labels_sf6[site_order_sf6], pooled_label_sf6)),
    score_label = format_score(score),
    abs_log_or  = abs(log(or_int))
  )

## plot and save ---------------------------------------------------------------

sf6 = ggplot(sf6_data, aes(x = site_label, y = or_int, ymin = or_int_lower, ymax = or_int_upper)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbar(width = 0.2, linewidth = 0.5, color = "gray30") +
  geom_point(aes(fill = abs_log_or, shape = is_pooled, size = is_pooled)) +
  facet_wrap(~score_label, nrow = 1) +
  scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 23), guide = "none") +
  scale_size_manual(values = c("FALSE" = 3, "TRUE" = 5), guide = "none") +
  scale_fill_viridis_c(option = "magma", begin = 0.2, end = 0.9, guide = "none") +
  scale_y_continuous(
    trans  = "log",
    breaks = c(0.8, 0.9, 1.0, 1.1, 1.2),
    labels = function(x) sprintf("%.1f", x),
    limits = c(0.75, 1.25)
  ) +
  coord_flip(clip = "off") +
  labs(
    x = NULL, y = "Interaction OR (95% CI)",
    caption = "← Weaker in cancer
OR = 1 (no difference)
Stronger in cancer →"
  ) +
  theme_ews() +
  theme(
    panel.grid.major.y = element_blank(),
    plot.caption = element_text(hjust = 0.5, size = 10, face = "italic", color = "gray40")
  )

sf6

ggsave(here("output", "figures", "figure_s06_interaction_forest.pdf"), sf6, width = 14, height = 6)


# ==============================================================================
# SUPP FIGURE 7: Sensitivity Analysis AUROC Differences
# ==============================================================================

## data prep -------------------------------------------------------------------

# Use VARIANT_N for proper encounter counts per analysis
if (exists("VARIANT_N")) {
  sf7_labels = c(
    main           = sprintf("Main (n=%s)", format(VARIANT_N["main"], big.mark = ",")),
    one_enc_per_pt = sprintf("One enc/patient (n=%s)", format(VARIANT_N["one_enc_per_pt"], big.mark = ",")),
    win0_96h       = sprintf("0-96h window (n=%s)", format(VARIANT_N["win0_96h"], big.mark = ",")),
    fullcode_only  = sprintf("Full code only (n=%s)", format(VARIANT_N["fullcode_only"], big.mark = ",")),
    no_ed_req      = sprintf("No ED req (n=%s)", format(VARIANT_N["no_ed_req"], big.mark = ","))
  )
} else {
  sf7_labels = c(
    main           = "Main",
    one_enc_per_pt = "One encounter/patient",
    win0_96h       = "0-96h window",
    fullcode_only  = "Full code only",
    no_ed_req      = "No ED requirement"
  )
}

sf7_wide = auroc_results_final |>
  fsubset(metric == "Encounter max") |>
  fselect(score_name, ca_01, analysis, auroc, auroc_se) |>
  pivot_wider(names_from = ca_01, values_from = c(auroc, auroc_se), names_sep = "_") |>
  fmutate(
    diff       = auroc_1 - auroc_0,
    diff_se    = sqrt(auroc_se_1^2 + auroc_se_0^2),
    diff_lower = diff - 1.96 * diff_se,
    diff_upper = diff + 1.96 * diff_se,
    score_label = format_score(score_name),
    analysis_label = factor(
      analysis,
      levels = names(sf7_labels),
      labels = sf7_labels
    )
  )

## plot and save ---------------------------------------------------------------

sf7 = ggplot(sf7_wide, aes(x = analysis_label, y = diff, ymin = diff_lower, ymax = diff_upper)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(width = 0.2, linewidth = 0.5, color = "gray30") +
  geom_point(aes(fill = abs(diff)), shape = 21, size = 3) +
  facet_wrap(~score_label, nrow = 1) +
  scale_fill_viridis_c(option = "magma", begin = 0.2, end = 0.9, guide = "none") +
  scale_y_continuous(
    limits = c(-0.10, 0.06),
    breaks = seq(-0.10, 0.06, 0.04),
    labels = function(x) sprintf("%.2f", x)
  ) +
  coord_flip() +
  labs(
    x = NULL, y = "AUROC Difference (Cancer − Non-Cancer)",
    caption = "Positive values indicate higher discrimination in cancer patients"
  ) +
  theme_ews() +
  theme(
    panel.grid.major.y = element_blank(),
    plot.caption = element_text(hjust = 0.5, size = 10, face = "italic", color = "gray40")
  )

sf7

ggsave(here("output", "figures", "figure_s07_sensitivity_auroc_diff.pdf"), sf7, width = 12, height = 5)


# ==============================================================================
# SUPP FIGURE 8: Hematologic vs Solid Tumor AUROC Comparison
# ==============================================================================

## data prep -------------------------------------------------------------------

sf8_noncancer_prep = maxscores_ca_raw |>
  fsubset(analysis == "main" & ca_01 == 0) |>
  fmutate(value = max_value) |>
  fgroup_by(score_name, value, outcome) |>
  fsummarise(n = fsum(n)) |>
  fungroup()

sf8_noncancer = lapply(unique(sf8_noncancer_prep$score_name), function(sc) {
  d = fsubset(sf8_noncancer_prep, score_name == sc)
  result = calc_auroc_from_counts(d)
  result$score_name = sc
  result$cancer_type = "Non-Cancer"
  result
}) |>
  bind_rows() |>
  fselect(score_name, auroc, ci_lower, ci_upper, cancer_type)

sf8_data = liquid_aurocs_final |>
  fselect(score_name, auroc, ci_lower, ci_upper, cancer_type) |>
  bind_rows(sf8_noncancer) |>
  fmutate(
    score_label = format_score(score_name),
    cancer_type = factor(cancer_type, levels = c("Non-Cancer", "Solid", "Hematologic"))
  )

pal_cancer_type = c("Non-Cancer" = "#4575b4", "Solid" = "#f4a582", "Hematologic" = "#d73027")

## plot and save ---------------------------------------------------------------

sf8 = ggplot(sf8_data, aes(x = score_label, y = auroc, fill = cancer_type)) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.2, linewidth = 0.5, color = "black",
    position = position_dodge(width = 0.7)
  ) +
  geom_point(
    shape = 21, color = "black", size = 4, stroke = 0.7,
    position = position_dodge(width = 0.7)
  ) +
  scale_fill_manual(values = pal_cancer_type, name = "Cohort") +
  scale_y_continuous(limits = c(0.62, 0.81), breaks = seq(0.65, 0.80, 0.05)) +
  labs(x = NULL, y = "AUROC (95% CI)") +
  theme_ews() +
  theme(legend.position = "top", axis.text.x = element_text(face = "bold", size = 11))

sf8

ggsave(here("output", "figures", "figure_s08_heme_solid_auroc.pdf"), sf8, width = 8, height = 5)


# ==============================================================================
# SUPP FIGURE 9: Score Co-Positivity (Mirrored UpSet-style)
# ==============================================================================

## data prep -------------------------------------------------------------------

score_cols = c("sirs", "qsofa", "mews", "mews_sf", "news")

upset_agg = upset |>
  fgroup_by(ca_01, sirs, qsofa, mews, mews_sf, news) |>
  fsummarise(n = fsum(n)) |>
  fungroup() |>
  fmutate(
    combo = paste(sirs, qsofa, mews, mews_sf, news, sep = "-"),
    combo_label = pmap_chr(
      list(sirs, qsofa, mews, mews_sf, news),
      function(s, q, m, ms, n) {
        scores = c("SIRS", "QSOFA", "MEWS", "MEWS-SF", "NEWS")[c(s, q, m, ms, n) == 1]
        if (length(scores) == 0) return("None")
        paste(scores, collapse = " + ")
      }
    ),
    n_positive   = sirs + qsofa + mews + mews_sf + news,
    cohort_label = format_cohort(ca_01, COHORT_N)
  )

top_combos = upset_agg |>
  fsubset(n_positive > 0) |>
  fgroup_by(combo, combo_label, sirs, qsofa, mews, mews_sf, news, n_positive) |>
  fsummarise(total_n = fsum(n)) |>
  fungroup() |>
  arrange(-total_n) |>
  head(15) |>
  pull(combo)

upset_plot_data = upset_agg |>
  fsubset(combo %in% top_combos) |>
  fgroup_by(ca_01) |>
  fmutate(total_cohort = fsum(n), pct = n / total_cohort * 100) |>
  fungroup() |>
  fmutate(
    n_mirror   = fifelse(ca_01 == 1, -n, n),
    pct_mirror = fifelse(ca_01 == 1, -pct, pct)
  )

combo_order = upset_plot_data |>
  fgroup_by(combo_label) |>
  fsummarise(total = fsum(n)) |>
  arrange(-total) |>
  pull(combo_label)

upset_plot_data = upset_plot_data |>
  fmutate(combo_label = factor(combo_label, levels = rev(combo_order)))

sf9_pal = build_cohort_palette(levels(upset_plot_data$cohort_label))

## bar plot --------------------------------------------------------------------

p_bars = ggplot(upset_plot_data, aes(x = combo_label, y = n_mirror, fill = cohort_label)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = sf9_pal, name = "Cohort") +
  scale_y_continuous(labels = function(x) scales::comma(abs(x)), breaks = scales::pretty_breaks(n = 6)) +
  coord_flip() +
  labs(x = NULL, y = "Number of Encounters", title = "Score Co-Positivity Patterns") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position    = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_blank(),
    plot.title         = element_text(face = "bold", hjust = 0.5)
  )

## matrix plot -----------------------------------------------------------------

matrix_data = upset_plot_data |>
  fselect(combo_label, sirs, qsofa, mews, mews_sf, news) |>
  funique() |>
  pivot_longer(cols = c(sirs, qsofa, mews, mews_sf, news), names_to = "score", values_to = "positive") |>
  fmutate(
    score_label = factor(score,
                         levels = c("sirs", "qsofa", "mews", "mews_sf", "news"),
                         labels = c("SIRS", "QSOFA", "MEWS", "MEWS-SF", "NEWS"))
  )

p_matrix = ggplot(matrix_data, aes(x = score_label, y = combo_label)) +
  geom_point(aes(color = factor(positive)), size = 3) +
  scale_color_manual(values = c("0" = "gray80", "1" = "gray20"), guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid  = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 9)
  )

matrix_lines = matrix_data |>
  fsubset(positive == 1) |>
  fgroup_by(combo_label) |>
  fmutate(x_num = as.numeric(score_label), x_min = fmin(x_num), x_max = fmax(x_num)) |>
  fungroup() |>
  fsubset(x_min != x_max) |>
  fselect(combo_label, x_min, x_max) |>
  funique()

if (nrow(matrix_lines) > 0) {
  p_matrix = p_matrix +
    geom_segment(
      data = matrix_lines,
      aes(x = x_min, xend = x_max, y = combo_label, yend = combo_label),
      linewidth = 0.8, color = "gray20"
    )
}

## assemble and save -----------------------------------------------------------

sf9 = p_matrix + p_bars + plot_layout(widths = c(1, 3))

sf9

ggsave(here("output", "figures", "figure_s09_upset_mirrored.pdf"), sf9, width = 10, height = 8)

# ==============================================================================
# DONE
# ==============================================================================

message("All figures saved to: ", here("output", "figures"))
