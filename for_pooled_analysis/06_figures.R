# ==============================================================================
# 06_figures.R
# Publication-quality figures for pooled EWS analyses.
#
# Run order assumes 00_load -> 01_tables -> 02_discrimination -> 03_threshold
# -> 04_threshold_equivalence -> 05_subgroups -> 06_figures.
#
# Main-text figures:
#   Figure 1  Calibration curves (observed deterioration rate by score value)
#   Figure 2  Cumulative incidence of positivity
#   Figure 3  Efficiency curves (sensitivity vs threshold positivity, all 5 scores)
#
# Supplement figures (numbered in the order cited in Results):
#   Figure S1  Flow diagram (built in 01_tables.R)
#   Figure S2  Score distributions + risk differences (2-panel combined)
#   Figure S3  Score co-positivity (mirrored UpSet-style)
#   Figure S4  Component contributions to score positivity (dumbbell)
#   Figure S5  Threshold performance characteristics (sens, spec, PPV, NPV
#              across integer thresholds; 5 scores x 4 measures grid)
#   Figure S6  AUROC heterogeneity: site-level + cancer subtype (2-panel combined)
#   Figure S7  Robustness of cancer-noncancer AUROC gap (sensitivity analyses,
#              24h naive, and 24h clustered bootstrap, unified rows)
# ==============================================================================

# SETUP ------------------------------------------------------------------------

library(patchwork)
library(ggplot2)
library(scales)
library(purrr)

if (!dir.exists(here("output", "figures"))) {
  dir.create(here("output", "figures"), recursive = TRUE)
}

## shared theme and palette ---------------------------------------------------

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

## score labels ---------------------------------------------------------------

score_levels = c("sirs", "qsofa", "mews", "mews_sf", "news")
score_labels = c("SIRS", "QSOFA", "MEWS", "MEWS-SF", "NEWS")
names(score_labels) = score_levels

## helper functions -----------------------------------------------------------

format_score = function(x) {
  x = tolower(x)
  x = fifelse(x == "mews-sf", "mews_sf", x)
  factor(x, levels = score_levels, labels = score_labels)
}

format_cohort = function(ca_01, n_by_cohort = NULL) {
  if (!is.null(n_by_cohort)) {
    lab_noca = sprintf("Non-Cancer (n=%s)", format(n_by_cohort["0"], big.mark = ","))
    lab_ca   = sprintf("Cancer (n=%s)",     format(n_by_cohort["1"], big.mark = ","))
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

build_cohort_palette = function(cohort_labels) {
  setNames(c("#4575b4", "#d73027"), cohort_labels)
}

## environment checks ---------------------------------------------------------

if (!exists("COHORT_N"))  stop("COHORT_N not found. Source 00_load.R first.")
if (!exists("VARIANT_N")) message("WARNING: VARIANT_N not found")
if (!exists("SITE_N"))    message("WARNING: SITE_N not found")
message(
  "Using COHORT_N: ", format(COHORT_N["0"], big.mark = ","), " non-cancer, ",
  format(COHORT_N["1"], big.mark = ","), " cancer"
)

## weighted AUROC from pooled counts (used by some supp figures) --------------
# Implementation note: we coerce to plain data.frame and pull columns as base
# numeric vectors before indexing. Earlier versions using $n[idx] directly on
# data.table/tidytable inputs threw "subscript out of bounds" when the
# subsetted object retained non-standard column indexing semantics.

calc_auroc_from_counts = function(data) {
  
  data = as.data.frame(data)
  
  if (nrow(data) == 0) {
    return(tidytable(
      auroc    = NA_real_,
      se       = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      n_events = 0L,
      n_total  = 0L
    ))
  }
  
  events     = data[data$outcome == 1L, , drop = FALSE]
  non_events = data[data$outcome == 0L, , drop = FALSE]
  
  ev_vals  = as.numeric(events$value)
  ev_ns    = as.numeric(events$n)
  nev_vals = as.numeric(non_events$value)
  nev_ns   = as.numeric(non_events$n)
  
  n_events     = sum(ev_ns,  na.rm = TRUE)
  n_non_events = sum(nev_ns, na.rm = TRUE)
  
  if (n_events == 0 || n_non_events == 0) {
    return(tidytable(
      auroc    = NA_real_,
      se       = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      n_events = as.integer(n_events),
      n_total  = as.integer(n_events + n_non_events)
    ))
  }
  
  concordant = 0
  tied       = 0
  
  for (i in seq_along(ev_vals)) {
    sv = ev_vals[i]
    en = ev_ns[i]
    if (!is.finite(en)) next
    
    lower_mask = !is.na(nev_vals) & nev_vals < sv
    tied_mask  = !is.na(nev_vals) & nev_vals == sv
    
    concordant = concordant + en * sum(nev_ns[lower_mask], na.rm = TRUE)
    tied       = tied       + en * sum(nev_ns[tied_mask],  na.rm = TRUE)
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
# FIGURE 1: Calibration curves (observed deterioration rate by score value)
# ==============================================================================
#
# Observed deterioration rate at each integer score value, by cohort. The
# complementary risk-difference panel (with parallel-shift trend line) lives
# in Supp Figure 2.

message("\n== Building Figure 1 (calibration curves) ==")

fig1_data = rd_counts_final |>
  fmutate(
    score_lab  = format_score(score_name),
    cohort_lab = format_cohort(ca_01, COHORT_N)
  ) |>
  fsubset(n >= 100) |>
  # suppress visually distracting tiny non-cancer points at score=0 for MEWS
  # and MEWS-SF; retained in underlying data but hidden on the plot
  fsubset(!(ca_01 == 0L & max_value == 0L & score_name %in% c("mews", "mews_sf")))

fig1_pal = build_cohort_palette(levels(fig1_data$cohort_lab))

fig1_panel = function(score, show_y = TRUE, x_breaks = NULL) {
  
  d = fsubset(fig1_data, score_lab == score)
  
  p = ggplot(d, aes(x = max_value, y = p, color = cohort_lab)) +
    geom_line(linewidth = 0.6, alpha = 0.35) +
    geom_point(aes(size = n), alpha = 0.75) +
    facet_wrap(~score_lab) +
    scale_color_manual(values = fig1_pal) +
    scale_size_area(labels = label_comma(), breaks = c(1e5, 3e5, 5e5), max_size = 5) +
    scale_y_continuous(labels = label_percent(), limits = c(0, 0.9),
                       breaks = seq(0, 0.9, 0.15)) +
    theme_ews() +
    theme(panel.grid.major.x = element_blank()) +
    labs(x = NULL, y = if (show_y) "Observed Deterioration Rate" else NULL)
  
  if (!is.null(x_breaks)) {
    p = p + scale_x_continuous(
      breaks = x_breaks,
      limits = c(min(x_breaks), max(x_breaks))
    )
  }
  p
}

fig1_sirs   = fig1_panel("SIRS",    TRUE,  0:4)
fig1_qsofa  = fig1_panel("QSOFA",   FALSE, 0:3)
fig1_mews   = fig1_panel("MEWS",    FALSE, seq(0, 12, 2))
fig1_mewssf = fig1_panel("MEWS-SF", TRUE,  seq(0, 12, 2))
fig1_news   = fig1_panel("NEWS",    FALSE, seq(0, 16, 2))

## legend panel ---------------------------------------------------------------

fig1_legend_data = tidytable(
  cohort_lab = factor(levels(fig1_data$cohort_lab), levels = levels(fig1_data$cohort_lab)),
  n = c(1e5, 5e5), x = 1, y = 1
)

fig1_legend = ggplot(fig1_legend_data,
                     aes(x = x, y = y, color = cohort_lab, size = n)) +
  geom_point() +
  scale_color_manual(values = fig1_pal, name = "Cohort") +
  scale_size_area(labels = label_comma(), breaks = c(1e5, 3e5, 5e5), max_size = 5) +
  xlim(10, 20) + ylim(10, 20) +
  guides(
    color = guide_legend(title = "Cohort",        order = 1),
    size  = guide_legend(title = "Total Patients", order = 2)
  ) +
  theme_void() +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.5, 0.5),
    legend.box             = "vertical",
    legend.spacing.y       = unit(0.8, "cm")
  )

fig1 = (fig1_sirs | fig1_qsofa | fig1_mews) /
  (fig1_mewssf | fig1_news | fig1_legend) +
  plot_annotation(
    caption = "Score Value (encounter maximum)",
    theme   = theme(plot.caption = element_text(hjust = 0.5, size = 11))
  )

ggsave(here("output", "figures", "figure_01_calibration.pdf"),
       fig1, width = 10, height = 6)

# ==============================================================================
# FIGURE 2: Cumulative incidence of positivity
# ==============================================================================

message("\n== Building Figure 2 (cumulative incidence) ==")

fig2_data = cuminc_raw |>
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

fig2_pal = build_cohort_palette(levels(fig2_data$cohort_label))

fig2 = ggplot(fig2_data, aes(x = time_days, y = cum_inc, color = cohort_label)) +
  geom_step(linewidth = 0.7) +
  facet_wrap(~score_label, nrow = 1) +
  scale_color_manual(values = fig2_pal, name = "Cohort") +
  scale_x_continuous(breaks = seq(0, 7, 2)) +
  scale_y_continuous(labels = label_percent()) +
  labs(x = "Days from Ward Admission", y = "Cumulative Incidence of Positivity") +
  theme_ews() +
  theme(legend.position = "top")

ggsave(here("output", "figures", "figure_02_cuminc.pdf"),
       fig2, width = 10, height = 4)

# ==============================================================================
# FIGURE 3: Efficiency curves (alert rate vs sensitivity)
# ==============================================================================
#
# Each panel shows alert rate (y) against sensitivity (x) across integer
# thresholds, by cohort. Sensitivity-as-target framing: read rightward along x
# to "achieve X% sensitivity" and read y to see "% of encounters alerted".
# Dotted diagonal = no-skill reference (sens = alert rate).
#
# Standard-threshold markers:
#   - Circles (shape 21) mark each score's standard integer threshold
#     (SIRS >=2, qSOFA >=2, MEWS >=5, MEWS-SF >=7).
#   - Diamonds (shape 23) mark the field-conventional NEWS rule used in the
#     primary analysis: "aggregate >=5 OR any single parameter >=3". This rule
#     cannot be represented as a single integer threshold, so it sits off the
#     aggregate-only curve that the integer sweep traces.

message("\n== Building Figure 3 (efficiency curves) ==")

fig3_data = eff_data_final |>
  fmutate(
    score_lab  = format_score(score_name),
    cohort_lab = format_cohort(ca_01, COHORT_N)
  )

fig3_pal = build_cohort_palette(levels(fig3_data$cohort_lab))

# Integer-threshold standard operating points for SIRS, qSOFA, MEWS, MEWS-SF.
# NEWS is excluded here because its field-conventional rule incorporates the
# single-parameter escape clause and is plotted separately below.
fig3_std = fig3_data[is_std == TRUE & score_name != "news"]

# ---- Field-conventional NEWS operating point --------------------------------
#
# The "aggregate score >=5 OR any single parameter >=3" rule drives the NEWS
# numbers reported in Table 2. It is not a single integer threshold and
# therefore cannot be swept alongside the other scores; instead, we compute
# its sensitivity and alert rate once per cohort and plot it as a separate
# marker shape.
#
# These values must be pulled from the upstream object that produced the
# Table 2 NEWS row (generated in 01_tables.R from the field-conventional
# rule applied to the same cohort and outcome definitions used elsewhere).
# Do not hardcode from the manuscript; source from the computed object to
# keep figure and table in lockstep.
#
# Expected upstream object: `news_fc_operating` with columns
#   ca_01 (0 = non-cancer, 1 = cancer), sens, alert_rate
# If your upstream object is named differently, adjust the line below.

if (!exists("news_fc_operating")) {
  stop(
    "Figure 3 requires an upstream object `news_fc_operating` with NEWS ",
    "field-conventional (aggregate>=5 OR any parameter>=3) sensitivity and ",
    "alert_rate by cohort (ca_01 in 0/1). Produce this in 01_tables.R or ",
    "04_threshold_equivalence.R alongside the Table 2 NEWS computation."
  )
}

fig3_news_fc = news_fc_operating |>
  fmutate(
    score_name = "news",
    score_lab  = factor("NEWS", levels = levels(fig3_data$score_lab)),
    cohort_lab = format_cohort(ca_01, COHORT_N)
  )

fig3_panel = function(score, show_y = TRUE) {
  
  d       = fsubset(fig3_data,     score_lab == score)
  d_std   = fsubset(fig3_std,      score_lab == score)
  d_fc    = fsubset(fig3_news_fc,  score_lab == score)
  
  p = ggplot(d, aes(x = sens, y = alert_rate, color = cohort_lab)) +
    geom_abline(slope = 1, intercept = 0, color = "gray80", linetype = "dotted") +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    geom_point(size = 1.3, alpha = 0.7) +
    # Standard integer-threshold operating points (SIRS, qSOFA, MEWS, MEWS-SF):
    geom_point(data = d_std, size = 3, shape = 21, stroke = 0.6,
               aes(fill = cohort_lab), color = "black") +
    # Field-conventional NEWS operating point (NEWS panel only):
    geom_point(data = d_fc, size = 3.2, shape = 23, stroke = 0.6,
               aes(fill = cohort_lab), color = "black") +
    facet_wrap(~score_lab) +
    scale_color_manual(values = fig3_pal) +
    scale_fill_manual(values  = fig3_pal) +
    scale_x_continuous(labels = label_percent(), limits = c(0, 1),
                       breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(labels = label_percent(), limits = c(0, 1),
                       breaks = seq(0, 1, 0.25)) +
    theme_ews() +
    labs(x = "Sensitivity", y = if (show_y) "% Threshold-positive encounters" else NULL)
  
  p
}

fig3_sirs   = fig3_panel("SIRS",    TRUE)
fig3_qsofa  = fig3_panel("QSOFA",   FALSE)
fig3_mews   = fig3_panel("MEWS",    FALSE)
fig3_mewssf = fig3_panel("MEWS-SF", TRUE)
fig3_news   = fig3_panel("NEWS",    FALSE)

fig3_legend = ggplot(
  tidytable(
    cohort_lab = factor(levels(fig3_data$cohort_lab), levels = levels(fig3_data$cohort_lab)),
    x = 1, y = 1
  ),
  aes(x = x, y = y, color = cohort_lab)
) +
  geom_point(size = 3) +
  scale_color_manual(values = fig3_pal, name = "Cohort") +
  xlim(10, 20) + ylim(10, 20) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme_void() +
  theme(
    legend.position        = "inside",
    legend.position.inside = c(0.5, 0.5),
    legend.direction       = "vertical"
  )

fig3 = (fig3_sirs | fig3_qsofa | fig3_mews) /
  (fig3_mewssf | fig3_news | fig3_legend)

ggsave(here("output", "figures", "figure_03_efficiency_curves.pdf"),
       fig3, width = 12, height = 7)

# ==============================================================================
# SUPP FIGURE 2: Score distributions + risk differences (2-panel combined)
# ==============================================================================
#
# Companion to Figure 1. Two-panel vertical stack with shared x-axis facets
# (free_x scales since score ranges differ).
#
#   Panel A: distribution of encounter-level maximum score values by cohort.
#            Vertical bars (count on y, score on x), faceted by score,
#            dodged across cohorts.
#
#   Panel B: risk difference (cancer - noncancer) at each integer score value,
#            with Newcombe 95% CIs and a weighted-regression trend line
#            anchored at the data's weighted center of mass. A near-flat
#            trend indicates the gap is approximately constant across the
#            score range (parallel-shift argument).
#
# Layout rationale: stacking with shared score facets lets a reader scan
# vertically within one score to see "where do encounters sit" (top) and
# "what's the risk gap at each value" (bottom), without reorienting axes.

message("\n== Building Supp Figure 2 (score distributions + risk differences) ==")

## Panel A: histograms rotated (score value on x-axis) ------------------------

sf2a_data = maxscores_ca_raw |>
  fsubset(tolower(analysis) == "main") |>
  fgroup_by(score_name, ca_01, max_value) |>
  fsummarise(n = fsum(n)) |>
  fungroup() |>
  fmutate(
    score_lab    = format_score(score_name),
    cohort_label = format_cohort(ca_01, COHORT_N)
  ) |>
  fsubset(max_value <= 14)

sf2a_pal = build_cohort_palette(levels(sf2a_data$cohort_label))

sf2a = ggplot(sf2a_data, aes(x = max_value, y = n, fill = cohort_label)) +
  geom_col(position = position_dodge(width = 0.45), width = 0.4) +
  facet_wrap(~score_lab, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = sf2a_pal, name = "Cohort") +
  scale_y_continuous(labels = label_comma()) +
  scale_x_continuous(breaks = function(x) seq(floor(x[1]), ceiling(x[2]), by = 1L)) +
  labs(x = NULL, y = "Number of Encounters") +
  theme_ews() +
  theme(
    legend.position    = "top",
    panel.grid.major.x = element_blank(),
    strip.text         = element_text(face = "bold", size = 10)
  )

## Panel B: risk differences with trend lines --------------------------------

sf2b_data = rd_wide_final |>
  fsubset(has_both == TRUE) |>
  fmutate(score_lab = format_score(score_name))

sf2b_trend_anchors = sf2b_data[, .(
  rd_mean = weighted.mean(rd, w = n_nc + n_ca, na.rm = TRUE),
  x_mean  = weighted.mean(max_value, w = n_nc + n_ca, na.rm = TRUE),
  x_min   = min(max_value, na.rm = TRUE),
  x_max   = max(max_value, na.rm = TRUE)
), by = score_name]

sf2b_trend = merge(
  rd_trend_final[, .(score_name, slope_rd_per_pt)],
  sf2b_trend_anchors,
  by = "score_name"
)
sf2b_trend[, `:=`(
  y_at_xmin = rd_mean + slope_rd_per_pt * (x_min - x_mean),
  y_at_xmax = rd_mean + slope_rd_per_pt * (x_max - x_mean),
  score_lab = format_score(score_name)
)]

sf2b = ggplot(sf2b_data, aes(x = max_value, y = rd)) +
  geom_hline(yintercept = 0, color = "gray60", linewidth = 0.3) +
  geom_errorbar(aes(ymin = rd_lo, ymax = rd_hi), width = 0, color = "gray30") +
  geom_point(color = "black", size = 1.8) +
  geom_segment(
    data = sf2b_trend,
    aes(x = x_min, xend = x_max, y = y_at_xmin, yend = y_at_xmax),
    linetype = "dashed", color = "gray40", linewidth = 0.5
  ) +
  facet_wrap(~score_lab, nrow = 1, scales = "free_x") +
  scale_y_continuous(labels = label_percent(), limits = c(-0.1, 0.5)) +
  scale_x_continuous(breaks = function(x) seq(floor(x[1]), ceiling(x[2]), by = 1L)) +
  labs(x = "Score value (encounter maximum)",
       y = "Absolute risk difference\n(cancer − non-cancer)") +
  theme_ews() +
  theme(
    panel.grid.major.x = element_blank(),
    strip.text         = element_text(face = "bold", size = 10)
  )

## assemble: panel A on top (taller), panel B below (shorter) ----------------
# strip labels on panel A only (panel B shares same facet columns; redundant
# strip on B is dropped to save vertical space and reinforce shared x-axis).

sf2b_nostrip = sf2b + theme(strip.text = element_blank(),
                            strip.background = element_blank())

sf2 = (sf2a / sf2b_nostrip) +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 12)))

ggsave(here("output", "figures", "figure_s02_score_distributions.pdf"),
       sf2, width = 16, height = 7)

# ==============================================================================
# SUPP FIGURE 6 (Panel A): Forest plot of per-site AUROC differences
# ==============================================================================
# Prepared here; combined with Panel B (heme/solid) below into sf6.
#
# Design rationale: an earlier version of this panel plotted absolute
# site-level AUROCs (one row per hospital, two cohort points per row,
# faceted by score). That layout invited the reader to compare sites to
# one another, which is not the paper's claim. Here we plot, for each
# score, the site-level AUROC difference (cancer minus non-cancer) with
# 95% CIs, plus the pooled random-effects estimate as a diamond. The
# figure directly shows what the paper actually claims: that the
# cancer-vs-noncancer gap is directionally consistent across sites and
# is not driven by outlier hospitals.
#
# Pooled estimates come from auroc_diff_final (built in 02_discrimination.R
# via logit-scale random-effects meta-analysis), keeping the pooled diamond
# aligned with Table 2's reported cancer-vs-noncancer difference rather
# than recomputing from sf6a_aurocs' site-level counts.

message("\n== Building Supp Figure 6 Panel A (AUROC-difference forest) ==")

sf6a_prep = maxscores_ca_raw |>
  fsubset(analysis == "main") |>
  fmutate(value = max_value) |>
  fgroup_by(site, score_name, ca_01, value, outcome) |>
  fsummarise(n = fsum(n)) |>
  fungroup() |>
  as.data.frame()

sf6a_combos = unique(sf6a_prep[, c("site", "score_name", "ca_01")])

sf6a_list = vector("list", nrow(sf6a_combos))
for (i in seq_len(nrow(sf6a_combos))) {
  st = sf6a_combos$site[i]
  sc = sf6a_combos$score_name[i]
  ca = sf6a_combos$ca_01[i]
  d  = sf6a_prep[sf6a_prep$site == st & sf6a_prep$score_name == sc & sf6a_prep$ca_01 == ca, ,
                 drop = FALSE]
  auc_res = as.data.frame(calc_auroc_from_counts(d))
  sf6a_list[[i]] = data.frame(
    site       = st,
    score_name = sc,
    ca_01      = ca,
    auroc      = auc_res$auroc,
    se         = auc_res$se,
    n_events   = auc_res$n_events,
    n_total    = auc_res$n_total,
    stringsAsFactors = FALSE
  )
}
sf6a_aurocs = do.call(rbind, sf6a_list) |> as_tidytable()

# Compute per-site cancer-minus-noncancer AUROC differences with propagated
# SE. Each site contributes one row per score; sites with missing estimates
# in either cohort are dropped from the plot.

sf6a_ca = sf6a_aurocs |> fsubset(ca_01 == 1L) |>
  fmutate(auroc_ca = auroc, se_ca = se) |>
  fselect(site, score_name, auroc_ca, se_ca)

sf6a_nc = sf6a_aurocs |> fsubset(ca_01 == 0L) |>
  fmutate(auroc_nc = auroc, se_nc = se) |>
  fselect(site, score_name, auroc_nc, se_nc)

sf6a_site_diff = merge(sf6a_ca, sf6a_nc, by = c("site", "score_name"),
                       all.x = TRUE) |>
  fmutate(
    diff      = auroc_ca - auroc_nc,
    se_diff   = sqrt(se_ca^2 + se_nc^2),
    ci_lower  = diff - 1.96 * se_diff,
    ci_upper  = diff + 1.96 * se_diff,
    is_pooled = FALSE
  ) |>
  fsubset(!is.na(diff) & !is.na(se_diff) & is.finite(se_diff)) |>
  fselect(site, score_name, diff, se_diff, ci_lower, ci_upper, is_pooled)

# Pooled estimates: pull from auroc_diff_final, the random-effects
# meta-analysis output. Use encounter-level main analysis.
sf6a_pooled_diff = auroc_diff_final |>
  as_tidytable() |>
  fsubset(analysis == "main" & metric == "Encounter max") |>
  fmutate(
    site      = "pooled",
    diff      = diff_auc,
    ci_lower  = diff_auc - 1.96 * se_diff,
    ci_upper  = diff_auc + 1.96 * se_diff,
    is_pooled = TRUE
  ) |>
  fselect(site, score_name, diff, se_diff, ci_lower, ci_upper, is_pooled)

sf6a_all = rbind(sf6a_site_diff, sf6a_pooled_diff) |>
  fmutate(score_label = format_score(score_name))

# Order scores on the y-axis by descending pooled gap magnitude (most-negative
# diff at top so the largest cancer-vs-noncancer gap reads first).
sf6a_score_order = sf6a_pooled_diff |>
  arrange(diff) |>
  pull(score_name)

sf6a_all[, score_label := factor(
  format_score(score_name),
  levels = format_score(sf6a_score_order)
)]

# One row per score. Site-level points are jittered vertically within each
# score row for readability; the pooled diamond sits at the score's y-anchor.
sf6a_all[, y_pos := as.numeric(score_label)]
set.seed(20260421)
sf6a_all[is_pooled == FALSE, y_pos := y_pos + runif(.N, -0.18, 0.18)]

sf6a = ggplot(sf6a_all, aes(x = diff, y = y_pos)) +
  geom_vline(xintercept = 0, color = "gray70", linetype = "dashed", linewidth = 0.4) +
  geom_errorbarh(
    data = sf6a_all[is_pooled == FALSE],
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0, color = "gray65", linewidth = 0.3, alpha = 0.7
  ) +
  geom_point(
    data = sf6a_all[is_pooled == FALSE],
    color = "gray45", fill = "gray81", shape = 21, size = 1.8, stroke = 0.3
  ) +
  geom_errorbarh(
    data = sf6a_all[is_pooled == TRUE],
    aes(xmin = ci_lower, xmax = ci_upper),
    height = 0, color = "black", linewidth = 0.7
  ) +
  geom_point(
    data = sf6a_all[is_pooled == TRUE],
    color = "black", fill = "black", shape = 23, size = 3.2, stroke = 0.5
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(sf6a_all$score_label)),
    labels = levels(sf6a_all$score_label),
    expand = expansion(add = 0.5)
  ) +
  scale_x_continuous(
    limits = c(-0.10, 0.10),
    breaks = seq(-0.10, 0.10, 0.05),
    labels = function(x) sprintf("%+.2f", x)
  ) +
  labs(x = "AUROC difference (cancer - non-cancer)", y = NULL) +
  theme_ews() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(face = "bold")
  )

# ==============================================================================
# SUPP FIGURE 6 (Panel B): Hematologic vs solid tumor AUROC comparison
# ==============================================================================
# Prepared here, then combined with Panel A into sf6 below.
#
# Non-cancer baseline: meta-analyze site-level AUROCs from auroc_enc_raw to
# match the CI construction used for solid/hematologic subgroups (both use
# meta_analyze_aurocs via site-level DeLong SEs), so CIs are visually
# comparable across cancer-type categories. Using calc_auroc_from_counts on
# pooled counts would produce CIs one order of magnitude narrower because
# between-site heterogeneity would not be reflected.

message("\n== Building Supp Figure 6 Panel B (heme vs solid AUROC) ==")

sf6b_nc_combos = unique(auroc_enc_raw[analysis == "main" & ca_01 == 0L,
                                      .(score_name)])

sf6b_nc_list = vector("list", nrow(sf6b_nc_combos))
for (i in seq_len(nrow(sf6b_nc_combos))) {
  sc  = sf6b_nc_combos$score_name[i]
  sub = auroc_enc_raw[analysis == "main" & ca_01 == 0L & score_name == sc]
  ma  = meta_analyze_aurocs(sub)
  sf6b_nc_list[[i]] = data.frame(
    score_name  = sc,
    auroc       = ma$auroc,
    ci_lower    = ma$ci_lower,
    ci_upper    = ma$ci_upper,
    cancer_type = "Non-Cancer",
    stringsAsFactors = FALSE
  )
}

sf6b_noncancer = do.call(rbind, sf6b_nc_list)

sf6b_data = liquid_aurocs_final |>
  fselect(score_name, auroc, ci_lower, ci_upper, cancer_type) |>
  as.data.frame() |>
  rbind(sf6b_noncancer) |>
  as_tidytable() |>
  fmutate(
    score_label = format_score(score_name),
    cancer_type = factor(cancer_type, levels = c("Hematologic", "Solid", "Non-Cancer"))
  )

# Palette: Non-Cancer blue matches Panel A; Solid/Hematologic pair with warm
# tones to signal the cancer subgroup split.
sf6b_pal = c("Non-Cancer" = "#4575b4", "Solid" = "#f4a582", "Hematologic" = "#d73027")

# Faceted by score with coord_flip, to mirror Panel A's layout. Each facet
# shows three horizontal rows (Non-Cancer, Solid, Hematologic) with AUROC on
# the shared x-axis. Y-axis range standardized to 0.60-0.90 to match Panel A.
sf6b = ggplot(sf6b_data, aes(x = cancer_type, y = auroc, color = cancer_type)) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.2, linewidth = 0.5
  ) +
  geom_point(size = 2.5) +
  facet_wrap(~score_label, nrow = 1) +
  scale_color_manual(values = sf6b_pal, name = "Cancer type") +
  scale_y_continuous(limits = c(0.65, 0.85), breaks = seq(0.65, 0.85, 0.05)) +
  coord_flip() +
  labs(x = NULL, y = "AUROC (95% CI)") +
  theme_ews() +
  theme(legend.position = "top", panel.grid.major.y = element_blank())

## Assemble Supp Figure 6: Panel A (site differences) / Panel B (heme-solid)
# Panel A is now a single-panel forest plot with one row per score; Panel B
# is a 5-facet strip. Use equal vertical heights since both panels convey a
# compact horizontal strip of estimates rather than a dense grid. Keep
# Panel B's facet strips visible since the new Panel A doesn't share facets.

sf6 = (sf6a / sf6b) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "A",
                  theme = theme(plot.tag = element_text(face = "bold", size = 12)))

ggsave(here("output", "figures", "figure_s06_auroc_heterogeneity.pdf"),
       sf6, width = 11, height = 7)

# ==============================================================================
# SUPP FIGURE 7: Robustness of the cancer-noncancer AUROC gap
# ==============================================================================
#
# Unified view of how the cancer-noncancer AUROC difference behaves under
# varying conditions. Each row on the y-axis is a condition; within each
# score facet, the point is the point estimate of (cancer - non-cancer)
# AUROC and the error bar is the 95% CI. Dashed vertical line at 0.
#
# Rows (bottom to top after coord_flip):
#   - Main analysis (primary sample)
#   - One encounter per patient
#   - 0-96 hour observation window
#   - Full code only
#   - No ED admission requirement
#   - 24-hour horizon, naive (treats repeated obs as independent)
#   - 24-hour horizon, clustered bootstrap (accounts for non-independence)
#
# The two 24-hour rows illustrate how accounting for within-encounter
# correlation widens the CI enough that it crosses zero, per the Results.

message("\n== Building Supp Figure 7 (robustness of AUROC gap) ==")

## Row labels for variants --------------------------------------------------

if (exists("VARIANT_N")) {
  sf7_variant_labels = c(
    main           = sprintf("Main (n=%s)",            format(VARIANT_N["main"],           big.mark = ",")),
    one_enc_per_pt = sprintf("One enc/patient (n=%s)", format(VARIANT_N["one_enc_per_pt"], big.mark = ",")),
    win0_96h       = sprintf("0-96h window (n=%s)",    format(VARIANT_N["win0_96h"],       big.mark = ",")),
    fullcode_only  = sprintf("Full code only (n=%s)",  format(VARIANT_N["fullcode_only"],  big.mark = ",")),
    no_ed_req      = sprintf("No ED req (n=%s)",       format(VARIANT_N["no_ed_req"],      big.mark = ","))
  )
} else {
  sf7_variant_labels = c(
    main           = "Main",
    one_enc_per_pt = "One encounter/patient",
    win0_96h       = "0-96h window",
    fullcode_only  = "Full code only",
    no_ed_req      = "No ED requirement"
  )
}

sf7_h24_naive_label = "24h horizon (naive)"
sf7_h24_boot_label  = "24h horizon (bootstrap)"

## Part 1: sensitivity-analysis variant diffs (encounter-level) --------------
# Use auroc_results_final (site-level meta-analyzed AUROCs); compute the
# (cancer - non-cancer) diff per variant per score with analytic SEs.

sf7_sens = auroc_results_final |>
  fsubset(metric == "Encounter max") |>
  fselect(score_name, ca_01, analysis, auroc, auroc_se) |>
  pivot_wider(names_from = ca_01, values_from = c(auroc, auroc_se), names_sep = "_") |>
  fmutate(
    diff       = auroc_1 - auroc_0,
    diff_se    = sqrt(auroc_se_1^2 + auroc_se_0^2),
    diff_lower = diff - 1.96 * diff_se,
    diff_upper = diff + 1.96 * diff_se,
    row_key    = analysis
  ) |>
  fselect(score_name, row_key, diff, diff_lower, diff_upper) |>
  as.data.table()

## Part 2: 24-hour horizon, naive (pooled counts, no clustering adjustment) --
# Compute per-cohort AUROC from counts_h24_raw using calc_auroc_from_counts,
# then take the (cancer - non-cancer) difference and combine analytic SEs.

sf7_h24n_pooled = counts_h24_raw |>
  fsubset(analysis == "main") |>
  fgroup_by(score_name, ca_01, value, outcome) |>
  fsummarise(n = fsum(n)) |>
  fungroup() |>
  as.data.frame()

sf7_h24n_combos = unique(sf7_h24n_pooled[, c("score_name", "ca_01")])

sf7_h24n_list = vector("list", nrow(sf7_h24n_combos))
for (i in seq_len(nrow(sf7_h24n_combos))) {
  sc = sf7_h24n_combos$score_name[i]
  ca = sf7_h24n_combos$ca_01[i]
  d  = sf7_h24n_pooled[sf7_h24n_pooled$score_name == sc & sf7_h24n_pooled$ca_01 == ca, ,
                       drop = FALSE]
  auc_res = as.data.frame(calc_auroc_from_counts(d))
  sf7_h24n_list[[i]] = data.frame(
    score_name = sc,
    ca_01      = ca,
    auroc      = auc_res$auroc,
    se         = auc_res$se,
    stringsAsFactors = FALSE
  )
}
sf7_h24n_aurocs = as_tidytable(do.call(rbind, sf7_h24n_list))

sf7_h24n = sf7_h24n_aurocs |>
  pivot_wider(names_from = ca_01, values_from = c(auroc, se), names_sep = "_") |>
  fmutate(
    diff       = auroc_1 - auroc_0,
    diff_se    = sqrt(se_1^2 + se_0^2),
    diff_lower = diff - 1.96 * diff_se,
    diff_upper = diff + 1.96 * diff_se,
    row_key    = "h24_naive"
  ) |>
  fselect(score_name, row_key, diff, diff_lower, diff_upper) |>
  as.data.table()

## Part 3: 24-hour horizon, clustered bootstrap -----------------------------
# Defensive handling: boot_h24_raw's exact structure is not documented in the
# scripts we have. We attempt to auto-detect:
#   - If it already contains a per-score difference column (e.g. 'diff_auc'
#     or 'diff') with an 'iter' column, compute 2.5/50/97.5 percentiles.
#   - Otherwise, if it has per-cohort per-iter AUROC columns, pivot and diff.
#   - Otherwise, if it has counts per iter, compute per-iter AUROC via
#     calc_auroc_from_counts and diff across ca_01.
# If none of these match, we fall back to the naive 24h CI with a message.

sf7_h24b = data.table()

compute_boot_h24 = function(boot_dt) {
  
  if (!is.data.table(boot_dt)) boot_dt = as.data.table(boot_dt)
  if (nrow(boot_dt) == 0) return(data.table())
  
  # Defensive: boot_h24_raw is not cleaned by clean_score_names() in 00_load.R,
  # so names may carry the "_total" suffix (e.g., "sirs_total"). Strip it now
  # so downstream format_score() matches the clean levels ("sirs", "qsofa", ...).
  if ("score_name" %in% names(boot_dt)) {
    boot_dt[, score_name := str_remove(score_name, "_total$")]
  }
  
  nms = names(boot_dt)
  iter_col = intersect(c("iter", "boot_iter", "bootstrap_iter", "b"), nms)[1]
  
  if (is.na(iter_col)) {
    message("    boot_h24_raw has no recognizable iteration column; skipping bootstrap row.")
    return(data.table())
  }
  
  # Case A: per-iter diff already computed (a column named diff* that varies
  # across iterations within score).
  diff_col = intersect(c("diff_auc", "diff", "auroc_diff", "auc_diff"), nms)[1]
  if (!is.na(diff_col)) {
    out = boot_dt[, .(
      diff       = median(get(diff_col), na.rm = TRUE),
      diff_lower = quantile(get(diff_col), 0.025, na.rm = TRUE, names = FALSE),
      diff_upper = quantile(get(diff_col), 0.975, na.rm = TRUE, names = FALSE)
    ), by = score_name]
    out[, row_key := "h24_bootstrap"]
    return(out[, .(score_name, row_key, diff, diff_lower, diff_upper)])
  }
  
  # Case B: per-iter per-cohort AUROC column.
  auroc_col = intersect(c("auroc", "auc"), nms)[1]
  if (!is.na(auroc_col) && "ca_01" %in% nms) {
    wide = dcast(
      boot_dt,
      as.formula(paste(iter_col, "+ score_name ~ ca_01")),
      value.var = auroc_col
    )
    setnames(wide, c("0", "1"), c("auroc_nc", "auroc_ca"), skip_absent = TRUE)
    if (all(c("auroc_nc", "auroc_ca") %in% names(wide))) {
      wide[, diff := auroc_ca - auroc_nc]
      out = wide[, .(
        diff       = median(diff, na.rm = TRUE),
        diff_lower = quantile(diff, 0.025, na.rm = TRUE, names = FALSE),
        diff_upper = quantile(diff, 0.975, na.rm = TRUE, names = FALSE)
      ), by = score_name]
      out[, row_key := "h24_bootstrap"]
      return(out[, .(score_name, row_key, diff, diff_lower, diff_upper)])
    }
  }
  
  # Case C: per-iter counts (score_name, ca_01, value, outcome, n, iter).
  if (all(c("score_name", "ca_01", "value", "outcome", "n") %in% nms)) {
    iters = unique(boot_dt[[iter_col]])
    per_iter = vector("list", length(iters))
    for (k in seq_along(iters)) {
      sub = boot_dt[get(iter_col) == iters[k]]
      sub_agg = sub[, .(n = sum(n)), by = .(score_name, ca_01, value, outcome)]
      combos  = unique(sub_agg[, .(score_name, ca_01)])
      lst = vector("list", nrow(combos))
      for (j in seq_len(nrow(combos))) {
        sc = combos$score_name[j]; ca = combos$ca_01[j]
        d  = as.data.frame(sub_agg[score_name == sc & ca_01 == ca])
        auc_res = as.data.frame(calc_auroc_from_counts(d))
        lst[[j]] = data.table(
          score_name = sc, ca_01 = ca, auroc = auc_res$auroc
        )
      }
      pi = rbindlist(lst)
      pi[, boot_id := iters[k]]
      per_iter[[k]] = pi
    }
    big = rbindlist(per_iter)
    wide = dcast(big, boot_id + score_name ~ ca_01, value.var = "auroc")
    setnames(wide, c("0", "1"), c("auroc_nc", "auroc_ca"), skip_absent = TRUE)
    wide[, diff := auroc_ca - auroc_nc]
    out = wide[, .(
      diff       = median(diff, na.rm = TRUE),
      diff_lower = quantile(diff, 0.025, na.rm = TRUE, names = FALSE),
      diff_upper = quantile(diff, 0.975, na.rm = TRUE, names = FALSE)
    ), by = score_name]
    out[, row_key := "h24_bootstrap"]
    return(out[, .(score_name, row_key, diff, diff_lower, diff_upper)])
  }
  
  message("    boot_h24_raw columns: ", paste(nms, collapse = ", "))
  message("    -> structure not recognized; skipping bootstrap row.")
  data.table()
}

if (exists("boot_h24_raw") && nrow(boot_h24_raw) > 0) {
  sf7_h24b = compute_boot_h24(boot_h24_raw)
  if (nrow(sf7_h24b) > 0) {
    message("  24h bootstrap diff computed for ", nrow(sf7_h24b), " scores.")
  }
} else {
  message("  boot_h24_raw not available; skipping bootstrap row.")
}

## Combine all rows and build the figure ------------------------------------

sf7_all = rbindlist(list(sf7_sens, sf7_h24n, sf7_h24b),
                    use.names = TRUE, fill = TRUE)

# Belt-and-suspenders: normalize score_name before format_score() so any
# variant that slipped through with a "_total" suffix (or casing quirk)
# still matches the canonical score levels. Without this, unmatched names
# produce an NA facet on the plot.
sf7_all[, score_name := tolower(str_remove(score_name, "_total$"))]

# Row ordering (bottom to top after coord_flip): sensitivity variants, then
# 24h naive, then 24h bootstrap. Factor levels set so the main row is at
# the bottom and bootstrap at the top.
sf7_row_order = c(
  names(sf7_variant_labels),
  "h24_naive",
  "h24_bootstrap"
)
sf7_row_labels = c(
  sf7_variant_labels,
  h24_naive     = sf7_h24_naive_label,
  h24_bootstrap = sf7_h24_boot_label
)

sf7_all[, `:=`(
  score_label = format_score(score_name),
  row_label   = factor(row_key, levels = sf7_row_order, labels = sf7_row_labels[sf7_row_order]),
  row_group   = fcase(
    row_key %in% names(sf7_variant_labels), "Sensitivity analyses",
    row_key == "h24_naive",                 "24-hour horizon",
    row_key == "h24_bootstrap",             "24-hour horizon",
    default = NA_character_
  )
)]

# Guard against any remaining unrecognized score_name values that would
# otherwise create an NA facet. Report them explicitly instead of silently
# letting them fall out.
sf7_bad = sf7_all[is.na(score_label)]
if (nrow(sf7_bad) > 0) {
  message("  WARNING: dropping ", nrow(sf7_bad),
          " sf7 rows with unrecognized score_name: ",
          paste(unique(sf7_bad$score_name), collapse = ", "))
  sf7_all = sf7_all[!is.na(score_label)]
}

sf7 = ggplot(sf7_all,
             aes(x = row_label, y = diff,
                 ymin = diff_lower, ymax = diff_upper)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Error bars: suppress for the naive 24h row (assumes independence of
  # repeated observations and therefore produces misleadingly-tight CIs).
  # The bootstrap row is the honest uncertainty estimate for that horizon.
  geom_errorbar(
    data = function(d) d[row_key != "h24_naive"],
    width = 0.25, linewidth = 0.5, color = "gray25"
  ) +
  geom_point(size = 2.5, color = "gray25") +
  facet_wrap(~score_label, nrow = 1) +
  scale_y_continuous(
    limits = c(-0.10, 0.06),
    breaks = seq(-0.10, 0.06, 0.04),
    labels = function(x) sprintf("%.2f", x)
  ) +
  coord_flip() +
  labs(
    x = NULL, y = "AUROC Difference (Cancer − Non-Cancer)",
    caption = "Positive values indicate higher discrimination in cancer patients.\nThe 24-hour naive row omits a CI because treating repeated observations as independent produces misleadingly tight uncertainty; the clustered bootstrap row reports the honest CI for that horizon."
  ) +
  theme_ews() +
  theme(
    panel.grid.major.y = element_blank(),
    plot.caption       = element_text(hjust = 0.5, size = 9,
                                      face = "italic", color = "gray40")
  )

ggsave(here("output", "figures", "figure_s07_auroc_gap_robustness.pdf"),
       sf7, width = 14, height = 6)

# ==============================================================================
# SUPP FIGURE 3: Score co-positivity (mirrored UpSet-style)
# ==============================================================================

message("\n== Building Supp Figure 3 (co-positivity UpSet) ==")

score_cols = c("sirs", "qsofa", "mews", "mews_sf", "news")

upset_agg = upset_final |>
  fgroup_by(ca_01, sirs, qsofa, mews, mews_sf, news) |>
  fsummarise(n = fsum(n)) |>
  fungroup() |>
  fmutate(
    combo       = paste(sirs, qsofa, mews, mews_sf, news, sep = "-"),
    combo_label = pmap_chr(
      list(sirs, qsofa, mews, mews_sf, news),
      function(s, q, m, ms, nw) {
        scores = c("SIRS", "QSOFA", "MEWS", "MEWS-SF", "NEWS")[c(s, q, m, ms, nw) == 1]
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

sf3_pal = build_cohort_palette(levels(upset_plot_data$cohort_label))

## bar plot -------------------------------------------------------------------

p_bars = ggplot(upset_plot_data, aes(x = combo_label, y = n_mirror, fill = cohort_label)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = sf3_pal, name = "Cohort") +
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

## matrix plot ----------------------------------------------------------------

matrix_data = upset_plot_data |>
  fselect(combo_label, sirs, qsofa, mews, mews_sf, news) |>
  funique() |>
  pivot_longer(cols = c(sirs, qsofa, mews, mews_sf, news),
               names_to = "score", values_to = "positive") |>
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

sf3 = p_matrix + p_bars + plot_layout(widths = c(1, 3))

ggsave(here("output", "figures", "figure_s03_copositivity.pdf"),
       sf3, width = 10, height = 8)

# ==============================================================================
# SUPP FIGURE 4: Component contributions to score positivity (dumbbell)
# ==============================================================================
#
# For each score, the percentage of first-positive encounters in which a given
# component was contributing at the time of first positivity, by cohort.
# Dumbbell layout: one row per component, two points per row (cohort), segment
# length = percentage-point difference. Components share a single y-axis order
# across all five score panels (physiologic grouping: cardiovascular →
# respiratory → neurologic → other) so the same component sits at the same
# vertical position in every panel. Rows left blank where a score does not
# include the component. SpO2 (NEWS) and SpO2/FiO2 (MEWS-SF) are collapsed
# into a single "Oxygenation" row since no score uses both.

message("\n== Building Supp Figure 4 (component contributions) ==")

## pooled percentages across outcomes -----------------------------------------
# upset_comp_final is stratified by o_primary_01; pool by summing numerators
# and denominators (not averaging stratum percentages).

sfC_raw = upset_comp_final |>
  fgroup_by(ca_01, score, component) |>
  fsummarise(n = fsum(n), n_encs = fsum(n_encs)) |>
  fungroup()

## collapse oxygenation components --------------------------------------------
# spo2 (NEWS) and sf (SpO2/FiO2, MEWS-SF) measure the same physiologic domain
# and no score uses both; collapse into a single "oxygenation" row so the
# y-axis reads physiologically rather than score-specifically.

sfC_raw = sfC_raw |>
  fmutate(component = fifelse(component %in% c("spo2", "sf"), "oxygenation", component)) |>
  fgroup_by(ca_01, score, component) |>
  fsummarise(n = fsum(n), n_encs = fsum(n_encs)) |>
  fungroup()

sfC_data = sfC_raw |>
  fmutate(pct = 100 * n / n_encs)

## wide form + per-row diff ---------------------------------------------------

sfC_wide = sfC_data |>
  fselect(score, component, ca_01, pct) |>
  pivot_wider(names_from = ca_01, values_from = pct, names_prefix = "pct_") |>
  fmutate(diff = pct_1 - pct_0) |>
  as_tidytable()

## component labels and ordering ----------------------------------------------
# Order components physiologically (top to bottom on plot):
# cardiovascular → respiratory → neurologic → other.

component_labels = c(
  hr          = "Heart rate",
  sbp         = "Systolic BP",
  rr          = "Respiratory rate",
  oxygenation = "Oxygenation",
  gcs         = "Mental status",
  temp        = "Temperature",
  wbc         = "WBC"
)

component_order = names(component_labels)

sfC_wide = sfC_wide |>
  fmutate(
    component_lab = component_labels[component],
    component_lab = fifelse(is.na(component_lab), component, component_lab)
  )

## merge diff + labels back onto long form ------------------------------------

sfC_plot = sfC_data |>
  join(
    sfC_wide |> fselect(score, component, diff, component_lab),
    on = c("score", "component"),
    how = "left"
  ) |>
  fmutate(
    score_lab     = format_score(score),
    cohort_lab    = format_cohort(ca_01, COHORT_N),
    component_lab = factor(
      component_lab,
      levels = rev(component_labels[component_order])
    )
  )

# one label row per (score, component) — not per cohort
sfC_labels = sfC_wide |>
  fmutate(
    score_lab     = format_score(score),
    component_lab = factor(
      component_lab,
      levels = rev(component_labels[component_order])
    ),
    diff_txt      = sprintf("%+0.1f", diff)
  )

sfC_pal = build_cohort_palette(levels(sfC_plot$cohort_lab))

## single-row plot with shared y-axis -----------------------------------------
# By using facet_wrap with a single row and NOT dropping unused factor levels,
# every panel shows the same set of components in the same order. Panels for
# scores that don't include a given component will have an empty row at that
# position — this is the "extra white space" that enables visual alignment
# across panels.

sfC = ggplot(sfC_plot, aes(x = pct, y = component_lab)) +
  # segment connecting the two cohort points
  geom_segment(
    data = sfC_wide |>
      fmutate(
        score_lab     = format_score(score),
        component_lab = factor(
          component_lab,
          levels = rev(component_labels[component_order])
        )
      ),
    aes(x = pct_0, xend = pct_1, y = component_lab, yend = component_lab),
    inherit.aes = FALSE,
    color = "gray60", linewidth = 0.6
  ) +
  geom_point(aes(color = cohort_lab), size = 2.8) +
  geom_text(
    data = sfC_labels,
    aes(x = 105, y = component_lab, label = diff_txt),
    inherit.aes = FALSE,
    hjust = 0, size = 2.8, color = "gray30", family = "mono"
  ) +
  facet_wrap(~score_lab, nrow = 1, drop = FALSE) +
  scale_color_manual(values = sfC_pal, name = "Cohort") +
  scale_x_continuous(
    limits = c(0, 115), breaks = seq(0, 100, 25),
    labels = function(x) ifelse(x <= 100, paste0(x, "%"), "")
  ) +
  scale_y_discrete(drop = FALSE) +
  theme_ews() +
  theme(
    legend.position    = "top",
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 9),
    axis.title.x       = element_text(size = 10)
  ) +
  labs(x = "% of score-positive encounters", y = NULL)

ggsave(here("output", "figures", "figure_s04_components.pdf"),
       sfC, width = 18, height = 5)

# ==============================================================================
# SUPP FIGURE 5: Threshold performance characteristics
# ==============================================================================
#
# 5 x 4 grid: scores as rows, performance measures (sens, spec, PPV, NPV) as
# columns. Each small panel has two lines (cancer vs non-cancer) across the
# integer threshold range for that score. The grid structure absorbs the
# complexity: row scans show how one score's four measures trade off across
# thresholds; column scans show how the cancer/non-cancer gap behaves for
# the same measure across all five scores (the NPV column, in particular,
# makes the "negative screens mean less in oncology" point visually).
#
# Data source: sweep_final from 04_threshold_equivalence.R, which already
# carries sens/spec/ppv/alert_rate (and sufficient components for NPV) across
# every integer threshold, with Wilson CIs. For NEWS, thresholds here reflect
# the aggregate-score-only rule; the field-conventional aggregate-plus-any3
# rule used for Table 2 cannot be swept as an integer threshold and is
# noted in the caption.

message("\n== Building Supp Figure 5 (threshold performance grid) ==")

## Data prep: derive NPV from sweep counts, reshape to long ------------------
# sweep_final has tp/fp/tn/fn at every integer threshold per score per cohort,
# which lets us compute NPV alongside the sens/spec/ppv already present.

sf5_data = sweep_final |>
  fmutate(
    npv = fifelse((tn + fn) > 0, tn / (tn + fn), NA_real_)
  ) |>
  fselect(score_name, ca_01, threshold, sens, spec, ppv, npv,
          sens_lo, sens_hi, spec_lo, spec_hi, ppv_lo, ppv_hi) |>
  as.data.table()

# Wilson CI for NPV (not carried in sweep_final by default)
sf5_npv_ci = sweep_final[, {
  ci = wilson_ci(tn, tn + fn)
  .(threshold, score_name, ca_01, npv_lo = ci$lo, npv_hi = ci$hi)
}]
# wilson_ci is defined in 04_threshold_equivalence.R; if not in scope, a local
# fallback is added below.
if (!exists("wilson_ci")) {
  wilson_ci = function(x, n, conf = 0.95) {
    z = qnorm(1 - (1 - conf) / 2)
    p = fifelse(n > 0, x / n, NA_real_)
    denom  = 1 + z^2 / n
    center = (p + z^2 / (2 * n)) / denom
    halfw  = (z / denom) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
    list(lo = pmax(0, center - halfw), hi = pmin(1, center + halfw))
  }
  sf5_npv_ci = sweep_final[, {
    ci = wilson_ci(tn, tn + fn)
    .(threshold, score_name, ca_01, npv_lo = ci$lo, npv_hi = ci$hi)
  }]
}

sf5_data = merge(
  sf5_data, sf5_npv_ci,
  by = c("threshold", "score_name", "ca_01"),
  all.x = TRUE
)

# Pivot to long: one row per (score, cohort, threshold, measure).
sf5_long = sf5_data |>
  pivot_longer(
    cols = c(sens, spec, ppv, npv),
    names_to = "measure", values_to = "value"
  ) |>
  as.data.table()

sf5_ci_long = rbindlist(list(
  sf5_data[, .(score_name, ca_01, threshold,
               measure = "sens", lo = sens_lo, hi = sens_hi)],
  sf5_data[, .(score_name, ca_01, threshold,
               measure = "spec", lo = spec_lo, hi = spec_hi)],
  sf5_data[, .(score_name, ca_01, threshold,
               measure = "ppv",  lo = ppv_lo,  hi = ppv_hi)],
  sf5_data[, .(score_name, ca_01, threshold,
               measure = "npv",  lo = npv_lo,  hi = npv_hi)]
))

sf5_long = merge(sf5_long, sf5_ci_long,
                 by = c("score_name", "ca_01", "threshold", "measure"),
                 all.x = TRUE)

# Drop degenerate endpoint rows where the measure is undefined (e.g., PPV
# when no encounters are threshold-positive). Keep threshold = 0 as the
# pure-maxima/minima anchor: sens = 100%, spec = 0%, PPV = prevalence,
# NPV = undefined (imputed to 1 below). These anchor the left edge of
# each curve honestly.
sf5_long = sf5_long[!is.na(value) & is.finite(value) |
                      (threshold == 0L & measure == "npv")]

# NPV is mathematically undefined at threshold = 0 (no encounters
# classified as negative means 0/0); impute to 1 so the line connects
# to the left edge rather than leaving a gap. The confidence interval
# is set to NA since it cannot be computed from zero counts.
sf5_long[measure == "npv" & threshold == 0L & is.na(value),
         `:=`(value = 1, lo = NA_real_, hi = NA_real_)]

# Restrict to thresholds that actually occur in the data, AND where there
# are enough encounters above the threshold in BOTH cohorts to produce
# stable estimates. Without this, the tail of NEWS (thresholds 15-17+)
# renders with 3-10 encounters per cohort and Wilson CIs balloon; the
# resulting "non-monotonicity" in sensitivity is actually sampling noise
# from vanishing denominators, not a real signal.
#
# The filter uses n_pos (threshold-positive count) rather than raw n
# because that is the denominator that matters for sens/PPV stability.
# We require n_pos >= 50 in each cohort AT each threshold; thresholds
# that fail this in either cohort are trimmed from both so the curves
# remain comparable.

MIN_N_POS_PER_COHORT = 50L

sf5_n_pos = sweep_final[, .(score_name, ca_01, threshold, n_pos)]
sf5_n_pos_wide = dcast(sf5_n_pos, score_name + threshold ~ ca_01,
                       value.var = "n_pos")
setnames(sf5_n_pos_wide, c("0", "1"), c("n_pos_nc", "n_pos_ca"),
         skip_absent = TRUE)
sf5_n_pos_wide[, keep := !is.na(n_pos_nc) & !is.na(n_pos_ca) &
                 (n_pos_nc >= MIN_N_POS_PER_COHORT) &
                 (n_pos_ca >= MIN_N_POS_PER_COHORT)]

# Always keep threshold = 0 even if n_pos is zero-by-definition there.
sf5_n_pos_wide[threshold == 0L, keep := TRUE]

sf5_keep_map = sf5_n_pos_wide[keep == TRUE, .(score_name, threshold)]
sf5_long = merge(sf5_long, sf5_keep_map,
                 by = c("score_name", "threshold"))

# Drop NPV points where the "negative" class is too small to produce a
# meaningful estimate. At low thresholds, specificity near 0 means
# almost every encounter is still classified as positive, leaving a
# sliver of encounters in the negative class -- NPV is then dominated
# by whether those few happen to include a deteriorator, and can
# produce misleading dips (e.g., MEWS at threshold 1) that are sampling
# artifacts rather than real predictive behavior. Require specificity
# >= 1% to include an NPV point.
#
# Threshold = 0 is exempted because NPV there is imputed to 1 and
# anchors the left edge of the curve by convention.
sf5_long = sf5_long[!(measure == "npv" &
                        threshold > 0L &
                        !is.na(spec_hi) & spec_hi < 0.01)]

# Diagnostic: report PPV at threshold 0 vs observed outcome prevalence.
# These should match (PPV at t=0 is prevalence by construction). If they
# don't, something upstream is wrong.
sf5_ppv_at_0 = sf5_long[measure == "ppv" & threshold == 0L,
                        .(score_name, ca_01, ppv_at_0 = value)]
sf5_prevalence = sweep_final[threshold == 0L,
                             .(score_name, ca_01,
                               prev = n_events / n_total)]
sf5_check = merge(sf5_ppv_at_0, sf5_prevalence,
                  by = c("score_name", "ca_01"))
message("  PPV at threshold 0 vs cohort prevalence (should match):")
print(sf5_check[, .(score_name, ca_01,
                    ppv_at_0 = round(ppv_at_0, 3),
                    prev     = round(prev,     3))])

# Labels and ordering
sf5_long[, `:=`(
  score_lab   = format_score(score_name),
  cohort_lab  = format_cohort(ca_01, COHORT_N),
  measure_lab = factor(measure,
                       levels = c("sens", "spec", "ppv", "npv"),
                       labels = c("Sensitivity", "Specificity",
                                  "Positive predictive value",
                                  "Negative predictive value"))
)]

# Belt-and-suspenders: drop any rows with unmatched score_name (shouldn't
# happen since sweep_final comes from cleaned maxscores_ca_raw, but keep
# the defensive pattern for consistency with sf7).
sf5_long = sf5_long[!is.na(score_lab)]

sf5_pal = build_cohort_palette(levels(sf5_long$cohort_lab))

## Plot: facet_grid with measures as rows, scores as columns ----------------
# Scores as columns (left-to-right) mirrors the layout used throughout the
# rest of the paper's figures. Measures as rows (top-to-bottom) means each
# row uses its own native score ranges rather than a shared 0-14 x-axis
# that would squash SIRS (0-4) into a fifth of its panel. y-axis is 0-1 for
# all measures since they're all proportions.

sf5 = ggplot(sf5_long,
             aes(x = threshold, y = value,
                 color = cohort_lab, fill = cohort_lab)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 0.9) +
  facet_grid(measure_lab ~ score_lab, scales = "free_x", switch = "y") +
  scale_color_manual(values = sf5_pal, name = "Cohort") +
  scale_fill_manual(values = sf5_pal,  guide = "none") +
  scale_y_continuous(labels = label_percent(), limits = c(0, 1),
                     breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(breaks = function(x) {
    # Every integer when the range is small (SIRS, QSOFA); every other
    # integer when the range is larger (MEWS, MEWS-SF, NEWS) to avoid
    # label crowding within each score's panel width.
    rng = ceiling(x[2]) - floor(x[1])
    step = if (rng <= 6L) 1L else 2L
    seq(floor(x[1]), ceiling(x[2]), by = step)
  }) +
  labs(
    x = "Threshold (positive if score ≥ value)",
    y = NULL,
    caption = paste(
      "Shaded bands are 95% Wilson confidence intervals.",
      "For NEWS, the sweep uses the aggregate-score threshold only;",
      "the field-conventional \"any single parameter ≥3\" rule cannot",
      "be varied as an integer threshold and is not reflected here."
    )
  ) +
  theme_ews() +
  theme(
    legend.position  = "top",
    strip.text.x     = element_text(face = "bold", size = 10),
    strip.text.y     = element_text(face = "bold", size = 10, angle = 0),
    strip.placement  = "outside",
    panel.spacing.x  = unit(0.7, "lines"),
    panel.spacing.y  = unit(0.5, "lines"),
    plot.caption     = element_text(hjust = 0, size = 9,
                                    face = "italic", color = "gray40")
  )

ggsave(here("output", "figures", "figure_s05_threshold_performance.pdf"),
       sf5, width = 14, height = 10)

# ==============================================================================
# DONE
# ==============================================================================

message("\nAll figures saved to: ", here("output", "figures"))
