# ==============================================================================
# figures_clean.R
# Publication-quality figures for pooled EWS analyses
# ==============================================================================

# setup ------------------------------------------------------------------------

## additional libraries --------------------------------------------------------

library(patchwork)
library(ggplot2)
library(scales)

##  output directory -----------------------------------------------------------

if (!dir.exists(here("output", "figures"))) {
  dir.create(here("output", "figures"), recursive = TRUE)
}

## shared theme and palette ----------------------------------------------------

pal_cancer = c("Non-Cancer" = "#4575b4", "Cancer" = "#d73027")

theme_ews = function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor   = element_blank(),
      strip.text         = element_text(face = "bold", size = 10),
      strip.background   = element_rect(fill = "gray85", color = "black"),
      legend.position    = "none"
    )
}

## score labels (reused across figures) ----------------------------------------

score_levels = c("sirs", "qsofa", "mews", "mews_sf", "news")
score_labels = c("SIRS", "QSOFA", "MEWS", "MEWS-SF", "NEWS")
names(score_labels) = score_levels

## helpers: format factor vars -------------------------------------------------

### score 
format_score = function(x) {
  x = tolower(x)
  x = fifelse(x == "mews-sf", "mews_sf", x)
  factor(x, levels = score_levels, labels = score_labels)
}

### cohort 
format_cohort = function(ca_01) {
  factor(
    fifelse(ca_01 == 1, "Cancer", "Non-Cancer"),
    levels = c("Non-Cancer", "Cancer")
  )
}

# ==============================================================================
# FIGURE 1: Risk by Score Value ------------------------------------------------
# ==============================================================================

## data prep -------------------------------------------------------------------

# NOTE: assumes maxscores_ca_raw exists with columns:
#   analysis, score_name, ca_01, max_value, n, outcome

fig1_data = maxscores_ca_raw |>
  fsubset(tolower(analysis) == "main") |>
  fmutate(n_events = n * outcome) |>
  fgroup_by(score_name, ca_01, max_value) |>
  fsummarise(
    n_at_score = fsum(n),
    n_outcomes = fsum(n_events)
  ) |>
  fungroup() |>
  fmutate(
    prob         = n_outcomes / n_at_score,
    score_label  = format_score(score_name),
    cohort_label = format_cohort(ca_01)
  ) |>
  fsubset(n_outcomes >= 100)

## panel function --------------------------------------------------------------

fig1_panel = function(data, score, show_y = TRUE, x_breaks = NULL) {
  
  panel_data = fsubset(data, score_label == score)
  
  p = 
    ggplot(panel_data, aes(x = max_value, y = prob, color = cohort_label)) +
    geom_line(linewidth = 0.6, alpha = 0.3) +
    geom_point(aes(size = n_at_score), alpha = 0.75) +
    facet_wrap(~score_label) +
    scale_color_manual(values = pal_cancer) +
    scale_size_area(
      labels = label_comma(),
      breaks = c(100000, 300000, 500000)
    ) +
    scale_y_continuous(
      labels = label_percent(),
      limits = c(0, 0.9),
      breaks = seq(0, 0.9, 0.15)
    ) +
    theme_ews() +
    theme(panel.grid.major.x = element_blank()) +
    labs(x = NULL)
  
  # x-axis breaks
  if (!is.null(x_breaks)) {
    p = p + scale_x_continuous(breaks = x_breaks)
  }
  
  # y-axis label
  if (show_y) {
    p = p + labs(y = "Deterioration Risk (%)")
  } else {
    p = p + labs(y = NULL)
  }
  p
}

## build panels ----------------------------------------------------------------

fig1_sirs   = fig1_panel(fig1_data, "SIRS",    show_y = TRUE,  x_breaks = 0:4)
fig1_qsofa  = fig1_panel(fig1_data, "QSOFA",   show_y = FALSE, x_breaks = 0:3)
fig1_mews   = fig1_panel(fig1_data, "MEWS",    show_y = FALSE, x_breaks = seq(0, 10, 2))
fig1_mewssf = fig1_panel(fig1_data, "MEWS-SF", show_y = TRUE,  x_breaks = seq(0, 10, 2))
fig1_news   = fig1_panel(fig1_data, "NEWS",    show_y = FALSE, x_breaks = seq(0, 10, 2))

## legend panel ----------------------------------------------------------------

fig1_legend_data = tidytable(
  cohort_label = factor(c("Non-Cancer", "Cancer"), levels = c("Non-Cancer", "Cancer")),
  n_at_score   = c(100000, 500000),
  x            = 1,
  y            = 1
)

fig1_legend = 
  ggplot(fig1_legend_data, aes(x = x, y = y, color = cohort_label, size = n_at_score)) +
  geom_point() +
  scale_color_manual(values = pal_cancer) +
  scale_size_area(labels = label_comma(), breaks = c(100000, 300000, 500000)) +
  xlim(10, 20) +
  ylim(10, 20) +
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

## assemble --------------------------------------------------------------------

fig1 = (fig1_sirs | fig1_qsofa | fig1_mews) /
  (fig1_mewssf | fig1_news | fig1_legend) +
  plot_annotation(
    caption = "Maximum Score Value Within Encounter",
    theme   = theme(plot.caption = element_text(hjust = 0.5, size = 11))
  )

fig1

ggsave(
  here("output", "figures", "fig1_risk_by_score.pdf"),
  fig1,
  width  = 10,
  height = 6,
  device = quartz
)

# ==============================================================================
# FIGURE 2: AUROC Comparison (Main Analysis) -----------------------------------
# ==============================================================================

## data prep -------------------------------------------------------------------

# NOTE: assumes auroc_results_final exists with columns:
#   analysis, metric, score_name, ca_01, auroc, ci_lower, ci_upper

fig2_data = auroc_results_final |>
  fsubset(analysis == "main" & metric == "Encounter max") |>
  fmutate(
    score_label  = format_score(score_name),
    cohort_label = format_cohort(ca_01)
  )

## plot ------------------------------------------------------------------------

fig2 = 
  ggplot(fig2_data, aes(x = score_label, y = auroc, fill = cohort_label)) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width     = 0.2,
    linewidth = 0.5,
    color     = "black",
    position  = position_dodge(width = 0.6)
  ) +
  geom_point(
    shape    = 21,
    color    = "black",
    size     = 4,
    stroke   = 0.7,
    position = position_dodge(width = 0.6)
  ) +
  scale_fill_manual(values = pal_cancer, name = "Cohort") +
  scale_y_continuous(
    limits = c(0.58, 0.80),
    breaks = seq(0.60, 0.80, 0.05)
  ) +
  labs(
    x = NULL,
    y = "AUROC (95% CI)"
  ) +
  theme_ews() +
  theme(
    legend.position = "top",
    axis.text.x     = element_text(face = "bold", size = 11)
  )

fig2

ggsave(
  here("output", "figures", "fig2_auroc_main.pdf"),
  fig2,
  width  = 7,
  height = 5,
  device = quartz
)

# ==============================================================================
# FIGURE 3: Cumulative Incidence of Positivity ---------------------------------
# ==============================================================================

## data prep -------------------------------------------------------------------

# NOTE: assumes cuminc_raw exists with columns:
#   analysis, score, ca_01, time_bin_start, n_at_risk, n_became_pos, o_primary_01

fig3_data = 
  fsubset(cuminc_raw, analysis == "main" & o_primary_01 == 1) |>
  fgroup_by(score, ca_01, time_bin_start) |>
  fsummarise(
    n_at_risk    = fsum(n_at_risk),
    n_became_pos = fsum(n_became_pos)
  ) |>
  fungroup() |>
  fmutate(
    cum_inc      = n_became_pos / n_at_risk,
    score_label  = format_score(score),
    cohort_label = format_cohort(ca_01),
    time_days    = time_bin_start / 24
  )

## plot ------------------------------------------------------------------------

fig3 = 
  ggplot(fig3_data, aes(x = time_days, y = cum_inc, color = cohort_label)) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~score_label, nrow = 1) +
  scale_color_manual(values = pal_cancer, name = "Cohort") +
  scale_x_continuous(breaks = seq(0, 7, 2)) +
  scale_y_continuous(labels = label_percent()) +
  labs(
    x = "Days from Ward Admission",
    y = "Cumulative Incidence of Positivity"
  ) +
  theme_ews() +
  theme(legend.position = "top")

fig3

ggsave(
  here("output", "figures", "fig3_cuminc.pdf"),
  fig3,
  width  = 10,
  height = 4,
  device = quartz
)

# ==============================================================================
# SUPP FIGURE 2: Score Distribution Histograms (Mirrored) ---------------------
# ==============================================================================

sf2_data = 
  maxscores_ca_raw |>
  fsubset(tolower(analysis) == "main") |>
  fgroup_by(score_name, ca_01, max_value) |>
  fsummarise(n = fsum(n)) |>
  fmutate(
    score_label  = format_score(score_name),
    cohort_label = format_cohort(ca_01),
    n_mirror     = fifelse(ca_01 == 1, -n, n)
  ) |>
  fsubset(max_value <= 11)

## plot ------------------------------------------------------------------------

sf2 = 
  ggplot(sf2_data, aes(x = n_mirror, y = factor(max_value), fill = cohort_label)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "black") +
  facet_wrap(~score_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = pal_cancer, name = "Cohort") +
  scale_x_continuous(labels = function(x) label_comma()(abs(x))) +
  labs(
    y = "Maximum Score Value",
    x = "Number of Encounters"
  ) +
  theme_ews() +
  theme(legend.position = "top")

sf2

ggsave(
  here("output", "figures", "sf2_score_histograms.pdf"),
  sf2,
  width  = 4,
  height = 11,
  device = quartz
)


# ==============================================================================
# FIGURE 4: Threshold-Performance Plots (SIRS and NEWS) ------------------------
# ==============================================================================

## data prep -------------------------------------------------------------------

# Calculate Se/Sp/PPV/NPV at each threshold for each score and cohort
# Threshold = "positive if score >= threshold"

fig4_prep = maxscores_ca_raw |>
  fsubset(tolower(analysis) == "main") |>
  fgroup_by(score_name, ca_01, max_value, outcome) |>
  fsummarise(n = fsum(n)) |>
  fungroup()

# Get all possible thresholds per score
fig4_thresholds = fig4_prep |>
  fgroup_by(score_name, ca_01) |>
  fsummarise(
    min_val = fmin(max_value),
    max_val = fmax(max_value)
  ) |>
  fungroup()

# Function to calculate metrics at each threshold
calc_threshold_metrics = function(data, score, cohort) {
  
  d = fsubset(data, score_name == score & ca_01 == cohort)
  
  # Total outcomes and non-outcomes
  totals = d |>
    fgroup_by(outcome) |>
    fsummarise(n = fsum(n)) |>
    fungroup()
  
  total_pos = totals$n[totals$outcome == 1]
  total_neg = totals$n[totals$outcome == 0]
  
  if (length(total_pos) == 0) total_pos = 0
  if (length(total_neg) == 0) total_neg = 0
  
  # Get range of thresholds
  thresholds = sort(unique(d$max_value))
  
  # Helper to safely sum
  safe_sum = function(x) {
    if (length(x) == 0) return(0)
    s = fsum(x)
    if (length(s) == 0 || is.na(s)) return(0)
    s
  }
  
  # Calculate metrics at each threshold
  results = lapply(thresholds, function(thresh) {
    
    # Positive = score >= threshold
    above = fsubset(d, max_value >= thresh)
    below = fsubset(d, max_value < thresh)
    
    tp = safe_sum(above$n[above$outcome == 1])
    fp = safe_sum(above$n[above$outcome == 0])
    fn = safe_sum(below$n[below$outcome == 1])
    tn = safe_sum(below$n[below$outcome == 0])
    
    sens = ifelse(tp + fn > 0, tp / (tp + fn), NA)
    spec = ifelse(tn + fp > 0, tn / (tn + fp), NA)
    ppv  = ifelse(tp + fp > 0, tp / (tp + fp), NA)
    npv  = ifelse(tn + fn > 0, tn / (tn + fn), NA)
    
    tidytable(
      score_name = score,
      ca_01      = cohort,
      threshold  = thresh,
      sens       = sens,
      spec       = spec,
      ppv        = ppv,
      npv        = npv
    )
  })
  
  bind_rows(results)
}

# Calculate for all score/cohort combinations
fig4_combos = expand.grid(
  score  = unique(fig4_prep$score_name),
  cohort = unique(fig4_prep$ca_01),
  stringsAsFactors = FALSE
)

fig4_metrics = map2(
  fig4_combos$score,
  fig4_combos$cohort,
  ~calc_threshold_metrics(fig4_prep, .x, .y)
) |>
  bind_rows()

# Reshape to long format for plotting
fig4_data = fig4_metrics |>
  pivot_longer(
    cols      = c(sens, spec, ppv, npv),
    names_to  = "metric",
    values_to = "value"
  ) |>
  fmutate(
    # Impute NPV to 1 at threshold 0 (undefined but functionally 100%)
    value        = fifelse(metric == "npv" & is.na(value) & threshold == 0, 1, value),
    score_label  = format_score(score_name),
    cohort_label = format_cohort(ca_01),
    metric_label = factor(
      metric,
      levels = c("sens", "spec", "ppv", "npv"),
      labels = c("Sensitivity", "Specificity", "PPV", "NPV")
    )
  )

# Filter to SIRS and NEWS for main figure
fig4_data_main = fsubset(fig4_data, score_label %in% c("SIRS", "NEWS"))

## panel function --------------------------------------------------------------

fig4_panel = function(data, score, show_y = TRUE) {
  
  panel_data = fsubset(data, score_label == score)
  
  p = ggplot(panel_data, aes(x = threshold, y = value, color = cohort_label)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.5) +
    facet_wrap(~metric_label, ncol = 1) +
    scale_color_manual(values = pal_cancer, name = "Cohort") +
    scale_y_continuous(labels = label_percent(), limits = c(0, 1)) +
    labs(x = NULL) +
    theme_ews() +
    theme(plot.title = element_text(hjust = 0, face = "bold", size = 12))
  
  if (show_y) {
    p = p + labs(y = NULL, title = score)
  } else {
    p = p + labs(y = NULL, title = score) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  
  p
}

## build panels ----------------------------------------------------------------

fig4_sirs = fig4_panel(fig4_data_main, "SIRS", show_y = TRUE)
fig4_news = fig4_panel(fig4_data_main, "NEWS", show_y = FALSE)

## shared legend ---------------------------------------------------------------

fig4_legend = ggplot(
  tidytable(
    cohort_label = factor(c("Non-Cancer", "Cancer"), levels = c("Non-Cancer", "Cancer")),
    x = c(1, 2),
    y = c(1, 1)
  ), 
  aes(x = x, y = y, color = cohort_label)
) +
  geom_point(size = 3) +
  scale_color_manual(values = pal_cancer, name = "Cohort") +
  xlim(10, 20) +
  ylim(10, 20) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_void() +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.5, 0.5),
    legend.direction       = "horizontal"
  )

## assemble --------------------------------------------------------------------

fig4 = fig4_legend / (fig4_sirs | fig4_news) +
  plot_layout(heights = c(0.05, 1)) +
  plot_annotation(
    caption = "Threshold (Positive if Score ≥ Value)",
    theme   = theme(plot.caption = element_text(hjust = 0.5, size = 11))
  )

fig4

ggsave(
  here("output", "figures", "fig4_threshold_performance.pdf"),
  fig4,
  width  = 8,
  height = 9,
  device = quartz
)


# ==============================================================================
# SUPP FIGURE 5: Threshold-Performance Plots (QSOFA, MEWS, MEWS-SF) ------------
# ==============================================================================

## data prep -------------------------------------------------------------------

sf5_data = fsubset(fig4_data, score_label %in% c("QSOFA", "MEWS", "MEWS-SF"))

## build panels ----------------------------------------------------------------

sf5_qsofa  = fig4_panel(sf5_data, "QSOFA",   show_y = TRUE)
sf5_mews   = fig4_panel(sf5_data, "MEWS",    show_y = FALSE)
sf5_mewssf = fig4_panel(sf5_data, "MEWS-SF", show_y = FALSE)

## shared legend ---------------------------------------------------------------

sf5_legend = ggplot(
  tidytable(
    cohort_label = factor(c("Non-Cancer", "Cancer"), levels = c("Non-Cancer", "Cancer")),
    x = c(1, 2),
    y = c(1, 1)
  ), 
  aes(x = x, y = y, color = cohort_label)
) +
  geom_point(size = 3) +
  scale_color_manual(values = pal_cancer, name = "Cohort") +
  xlim(10, 20) +
  ylim(10, 20) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_void() +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.5, 0.5),
    legend.direction       = "horizontal"
  )

## assemble --------------------------------------------------------------------

sf5 = sf5_legend / (sf5_qsofa | sf5_mews | sf5_mewssf) +
  plot_layout(heights = c(0.05, 1)) +
  plot_annotation(
    caption = "Threshold (Positive if Score ≥ Value)",
    theme   = theme(plot.caption = element_text(hjust = 0.5, size = 11))
  )

sf5

ggsave(
  here("output", "figures", "sf5_threshold_performance_other.pdf"),
  sf5,
  width  = 11,
  height = 9,
  device = quartz
)


# ==============================================================================
# SUPP FIGURE 3: AUROC Comparison (24h Analysis) -------------------------------
# ==============================================================================

## data prep -------------------------------------------------------------------

# Pool counts across sites, then calculate AUROC using Wilcoxon-Mann-Whitney
# AUROC = P(score_event > score_nonevent) + 0.5 * P(score_event = score_nonevent)

sf3_pooled = counts_h24_raw |>
  fsubset(analysis == "main") |>
  fgroup_by(score_name, ca_01, value, outcome) |>
  fsummarise(n = fsum(n)) |>
  fungroup()

# Function to calculate AUROC from count data
calc_auroc_from_counts = function(data) {
  
  # Helper to safely sum (fsum returns integer(0) on empty input)
  # Also convert to numeric to avoid integer overflow
  safe_fsum = function(x) {
    result = fsum(as.numeric(x))
    if (length(result) == 0) return(0)
    result
  }
  
  # Separate events and non-events
  events     = fsubset(data, outcome == 1)
  non_events = fsubset(data, outcome == 0)
  
  # Total counts
  n_events     = safe_fsum(events$n)
  n_non_events = safe_fsum(non_events$n)
  
  # Calculate concordance
  # For each event score value, count non-events with lower scores and ties
  concordant = 0
  tied       = 0
  
  for (i in seq_len(nrow(events))) {
    ev_score = events$value[i]
    ev_n     = as.numeric(events$n[i])
    
    # Non-events with lower scores (concordant)
    lower_idx  = which(non_events$value < ev_score)
    concordant = concordant + ev_n * safe_fsum(non_events$n[lower_idx])
    
    # Non-events with same score (tied)
    tied_idx = which(non_events$value == ev_score)
    tied     = tied + ev_n * safe_fsum(non_events$n[tied_idx])
  }
  
  # AUROC = (concordant + 0.5 * tied) / (n_events * n_non_events)
  auroc = (concordant + 0.5 * tied) / (n_events * n_non_events)
  
  # SE using Hanley-McNeil approximation
  q1 = auroc / (2 - auroc)
  q2 = 2 * auroc^2 / (1 + auroc)
  se = sqrt((auroc * (1 - auroc) + (n_events - 1) * (q1 - auroc^2) + 
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

# Calculate for each score/cohort combination
sf3_combos = sf3_pooled |>
  fselect(score_name, ca_01) |>
  funique()

sf3_aurocs = map2(
  sf3_combos$score_name,
  sf3_combos$ca_01,
  function(sc, ca) {
    d = fsubset(sf3_pooled, score_name == sc & ca_01 == ca)
    result = calc_auroc_from_counts(d)
    result$score_name = sc
    result$ca_01 = ca
    result
  }
) |>
  bind_rows() |>
  fmutate(
    score_label  = format_score(score_name),
    cohort_label = format_cohort(ca_01)
  )

## plot ------------------------------------------------------------------------

sf3 = 
  ggplot(sf3_aurocs, aes(x = score_label, y = auroc, fill = cohort_label)) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width     = 0.2,
    linewidth = 0.5,
    color     = "black",
    position  = position_dodge(width = 0.6)
  ) +
  geom_point(
    shape    = 21,
    color    = "black",
    size     = 4,
    stroke   = 0.7,
    position = position_dodge(width = 0.6)
  ) +
  scale_fill_manual(values = pal_cancer, name = "Cohort") +
  scale_y_continuous(
    limits = c(0.62, 0.72),
    breaks = seq(0.62, 0.72, 0.02)
  ) +
  labs(
    x = NULL,
    y = "AUROC (95% CI)"
  ) +
  theme_ews() +
  theme(
    legend.position = "top",
    axis.text.x     = element_text(face = "bold", size = 11)
  )

sf3

ggsave(
  here("output", "figures", "sf3_auroc_24h.pdf"),
  sf3,
  width  = 7,
  height = 5,
  device = quartz
)

# ==============================================================================
# SUPP FIGURE 4: Site-level AUROC Caterpillar ----------------------------------
# ==============================================================================

## data prep -------------------------------------------------------------------

# Calculate site-level AUROCs from maxscores_ca_raw

sf4_prep = maxscores_ca_raw |>
  fsubset(analysis == "main") |>
  fmutate(value = max_value) |>
  fgroup_by(site, score_name, ca_01, value, outcome) |>
  fsummarise(n = fsum(n)) |>
  fungroup()

# Get all site/score/cohort combinations
sf4_combos = sf4_prep |>
  fselect(site, score_name, ca_01) |>
  funique()

# Calculate AUROC for each combination
sf4_aurocs = pmap(
  list(sf4_combos$site, sf4_combos$score_name, sf4_combos$ca_01),
  function(st, sc, ca) {
    d = fsubset(sf4_prep, site == st & score_name == sc & ca_01 == ca)
    result = calc_auroc_from_counts(d)
    result$site = st
    result$score_name = sc
    result$ca_01 = ca
    result
  }
) |>
  bind_rows()

# Create site ordering based on NEWS non-cancer AUROC
site_order = sf4_aurocs |>
  fsubset(score_name == "news" & ca_01 == 0) |>
  arrange(auroc) |>
  pull(site)

# Create anonymous site labels (A, B, C, ...)
site_labels = setNames(LETTERS[seq_along(site_order)], site_order)

# Calculate pooled estimates (across all sites)
sf4_pooled_prep = sf4_prep |>
  fgroup_by(score_name, ca_01, value, outcome) |>
  fsummarise(n = fsum(n)) |>
  fungroup()

sf4_pooled_combos = sf4_pooled_prep |>
  fselect(score_name, ca_01) |>
  funique()

sf4_pooled = map2(
  sf4_pooled_combos$score_name,
  sf4_pooled_combos$ca_01,
  function(sc, ca) {
    d = fsubset(sf4_pooled_prep, score_name == sc & ca_01 == ca)
    result = calc_auroc_from_counts(d)
    result$site = "pooled"
    result$score_name = sc
    result$ca_01 = ca
    result
  }
) |>
  bind_rows()

# Combine site-level and pooled
sf4_all = bind_rows(sf4_aurocs, sf4_pooled) |>
  fmutate(
    score_label  = format_score(score_name),
    cohort_label = format_cohort(ca_01),
    site_label   = fifelse(site == "pooled", "Pooled", site_labels[site]),
    site_label   = factor(site_label, levels = c("Pooled", LETTERS[seq_along(site_order)])),
    is_pooled    = site == "pooled"
  )

## plot ------------------------------------------------------------------------

sf4 = 
  ggplot(sf4_all, aes(x = site_label, y = auroc, color = cohort_label)) +
  geom_hline(
    data = sf4_all |> fsubset(is_pooled),
    aes(yintercept = auroc, color = cohort_label),
    linetype = "dashed",
    linewidth = 0.4,
    alpha = 0.5
  ) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width     = 0.3,
    linewidth = 0.5,
    position  = position_dodge(width = 0.6)
  ) +
  geom_point(
    aes(shape = is_pooled),
    size     = 2.5,
    position = position_dodge(width = 0.6)
  ) +
  facet_wrap(~score_label, nrow = 1) +
  scale_color_manual(values = pal_cancer, name = "Cohort") +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18), guide = "none") +
  scale_y_continuous(
    limits = c(0.60, 0.82),
    breaks = seq(0.60, 0.80, 0.05)
  ) +
  coord_flip() +
  labs(
    x = NULL,
    y = "AUROC (95% CI)"
  ) +
  theme_ews() +
  theme(
    legend.position    = "top",
    panel.grid.major.y = element_blank()
  )

sf4

ggsave(
  here("output", "figures", "sf4_auroc_by_site.pdf"),
  sf4,
  width  = 12,
  height = 6,
  device = quartz
)



# ==============================================================================
# SUPP FIGURE 6: Interaction Forest Plot ---------------------------------------
# ==============================================================================

## data prep -------------------------------------------------------------------

# NOTE: assumes forest_data_final exists from 04_meta.R with columns:
#   site, score, score_lab, or_int, or_int_lower, or_int_upper, is_pooled

# Create anonymous site labels
site_order_sf6 = sort(unique(forest_data_final$site[forest_data_final$site != "Pooled"]))
site_labels_sf6 = setNames(LETTERS[seq_along(site_order_sf6)], site_order_sf6)

sf6_data = forest_data_final |>
  fsubset(analysis == "main" | is.na(analysis)) |>
  fmutate(
    site_label  = fifelse(is_pooled, "Pooled", site_labels_sf6[site]),
    site_label  = factor(site_label, levels = c(LETTERS[seq_along(site_order_sf6)], "Pooled")),
    score_label = format_score(score),
    # Absolute distance from null (log scale)
    abs_log_or  = abs(log(or_int))
  )

# Get magma colors at the endpoints for annotation
magma_low  = viridisLite::magma(10, begin = 0.2, end = 0.9)[1]
magma_high = viridisLite::magma(10, begin = 0.2, end = 0.9)[10]

## plot ------------------------------------------------------------------------

sf6 = 
  ggplot(sf6_data, aes(x = site_label, y = or_int, ymin = or_int_lower, ymax = or_int_upper)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  geom_errorbar(
    width     = 0.2,
    linewidth = 0.5,
    color     = "gray30"
  ) +
  geom_point(
    aes(fill = abs_log_or, shape = is_pooled, size = is_pooled)
  ) +
  facet_wrap(~score_label, nrow = 1) +
  scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 23), guide = "none") +
  scale_size_manual(values = c("FALSE" = 3, "TRUE" = 5), guide = "none") +
  scale_fill_viridis_c(
    option = "magma",
    begin  = 0.2,
    end    = 0.9,
    name   = "|log(OR)|",
    guide  = "none"
  ) +
  scale_y_continuous(
    trans  = "log",
    breaks = c(0.8, 0.9, 1.0, 1.1, 1.2),
    labels = function(x) sprintf("%.1f", x),
    limits = c(0.75, 1.25)
  ) +
  coord_flip(clip = "off") +
  labs(
    x       = NULL,
    y       = "Interaction OR (95% CI)"
  ) +
  theme_ews() +
  theme(
    panel.grid.major.y = element_blank(),
    plot.margin        = margin(t = 10, r = 10, b = 25, l = 10)
  )

# Plain caption (color = intensity, not direction)
sf6 = sf6 +
  labs(
    caption = "← Weaker in cancer
OR = 1 (no difference)
Stronger in cancer →"
  ) +
  theme(
    plot.caption = element_text(hjust = 0.5, size = 10, face = "italic", color = "gray40")
  )

sf6

ggsave(
  here("output", "figures", "sf6_interaction_forest.pdf"),
  sf6,
  width  = 12,
  height = 6,
  device = quartz
)

# ==============================================================================
# SUPP FIGURE 7: Sensitivity Analysis AUROC Differences ------------------------
# ==============================================================================

## data prep -------------------------------------------------------------------

# Calculate AUROC difference (cancer - non-cancer) for each score × analysis

sf7_wide = auroc_results_final |>
  fsubset(metric == "Encounter max") |>
  fselect(score_name, ca_01, analysis, auroc, auroc_se) |>
  pivot_wider(
    names_from  = ca_01,
    values_from = c(auroc, auroc_se),
    names_sep   = "_"
  ) |>
  fmutate(
    diff        = auroc_1 - auroc_0,
    diff_se     = sqrt(auroc_se_1^2 + auroc_se_0^2),
    diff_lower  = diff - 1.96 * diff_se,
    diff_upper  = diff + 1.96 * diff_se,
    score_label = format_score(score_name),
    analysis_label = factor(
      analysis,
      levels = c("main", "one_enc_per_pt", "win0_96h", "fullcode_only", "no_ed_req"),
      labels = c("Main", "One encounter/patient", "0-96h window", "Full code only", "No ED requirement")
    )
  )

## plot ------------------------------------------------------------------------

sf7 = 
  ggplot(sf7_wide, aes(x = analysis_label, y = diff, ymin = diff_lower, ymax = diff_upper)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(
    width     = 0.2,
    linewidth = 0.5,
    color     = "gray30"
  ) +
  geom_point(
    aes(fill = abs(diff)),
    shape = 21,
    size  = 3
  ) +
  facet_wrap(~score_label, nrow = 1) +
  scale_fill_viridis_c(
    option = "magma",
    begin  = 0.2,
    end    = 0.9,
    guide  = "none"
  ) +
  scale_y_continuous(
    limits = c(-0.10, 0.06),
    breaks = seq(-0.10, 0.06, 0.04),
    labels = function(x) sprintf("%.2f", x)
  ) +
  coord_flip() +
  labs(
    x       = NULL,
    y       = "AUROC Difference (Cancer − Non-Cancer)",
    caption = "Positive values indicate higher discrimination in cancer patients"
  ) +
  theme_ews() +
  theme(
    panel.grid.major.y = element_blank(),
    plot.caption       = element_text(hjust = 0.5, size = 10, face = "italic", color = "gray40")
  )

sf7

ggsave(
  here("output", "figures", "sf7_sensitivity_auroc_diff.pdf"),
  sf7,
  width  = 12,
  height = 5,
  device = quartz
)


# ==============================================================================
# SUPP FIGURE 8: Hematologic vs Solid Tumor AUROC Comparison -------------------
# ==============================================================================

## data prep -------------------------------------------------------------------

# NOTE: assumes liquid_aurocs_final exists from 05_subgroups.R with columns:
#   score_name, liquid_01, auroc, ci_lower, ci_upper, cancer_type
# Need to calculate non-cancer AUROCs from pooled maxscores_ca_raw to match

# Calculate pooled non-cancer AUROCs using same method as liquid_aurocs
sf8_noncancer_prep = maxscores_ca_raw |>
  fsubset(analysis == "main" & ca_01 == 0) |>
  fmutate(value = max_value) |>
  fgroup_by(score_name, value, outcome) |>
  fsummarise(n = fsum(n)) |>
  fungroup()

sf8_noncancer_combos = unique(sf8_noncancer_prep$score_name)

sf8_noncancer = lapply(sf8_noncancer_combos, function(sc) {
  d = fsubset(sf8_noncancer_prep, score_name == sc)
  result = calc_auroc_from_counts(d)
  result$score_name = sc
  result$cancer_type = "Non-Cancer"
  result
}) |>
  bind_rows() |>
  fselect(score_name, auroc, ci_lower, ci_upper, cancer_type)

# Combine with heme/solid
sf8_data = liquid_aurocs_final |>
  fselect(score_name, auroc, ci_lower, ci_upper, cancer_type) |>
  bind_rows(sf8_noncancer) |>
  fmutate(
    score_label = format_score(score_name),
    cancer_type = factor(cancer_type, levels = c("Non-Cancer", "Solid", "Hematologic"))
  )

# Color palette - non-cancer blue, solid light red, heme dark red
pal_cancer_type = c("Non-Cancer" = "#4575b4", "Solid" = "#8A1F19", "Hematologic" = "#EDA09C")

## plot ------------------------------------------------------------------------

sf8 = 
  ggplot(sf8_data, aes(x = score_label, y = auroc, fill = cancer_type)) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width     = 0.2,
    linewidth = 0.5,
    color     = "black",
    position  = position_dodge(width = 0.2)
  ) +
  geom_point(
    shape    = 21,
    color    = "black",
    size     = 4,
    stroke   = 0.7,
    position = position_dodge(width = 0.2)
  ) +
  scale_fill_manual(values = pal_cancer_type, name = "Cohort") +
  scale_y_continuous(
    limits = c(0.62, 0.81),
    breaks = seq(0.60, 0.80, 0.05)
  ) +
  labs(
    x = NULL,
    y = "AUROC (95% CI)"
  ) +
  theme_ews() +
  theme(
    legend.position = "top",
    axis.text.x     = element_text(face = "bold", size = 11)
  )

sf8

ggsave(
  here("output", "figures", "sf8_heme_solid_auroc.pdf"),
  sf8,
  width  = 8,
  height = 5,
  device = quartz
)


