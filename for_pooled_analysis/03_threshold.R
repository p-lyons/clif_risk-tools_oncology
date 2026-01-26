# ==============================================================================
# 03_threshold.R
# Threshold-based analyses: sensitivity/specificity, time to positivity,
# cumulative incidence, component analysis (main cohort only)
# ==============================================================================

# setup ------------------------------------------------------------------------

library(ggplot2)

source(here::here("code_pooled", "00_load.R"))

# SENSITIVITY / SPECIFICITY AT STANDARD THRESHOLDS -----------------------------

message("\n== Sensitivity/Specificity at standard thresholds ==")

## aggregate across sites ------------------------------------------------------

sesp = sesp_raw[, .(
  n_total      = sum(n_total),
  n_outcome    = sum(n_outcome),
  n_no_outcome = sum(n_no_outcome),
  n_pos        = sum(n_pos),
  n_neg        = sum(n_neg),
  tp           = sum(tp),
  fp           = sum(fp),
  tn           = sum(tn),
  fn           = sum(fn)
), by = .(score_name, ca_01)]

## calculate pooled metrics ----------------------------------------------------

sesp[, `:=`(
  sensitivity = tp / (tp + fn),
  specificity = tn / (tn + fp),
  ppv         = tp / (tp + fp),
  npv         = tn / (tn + fn),
  prevalence  = n_outcome / n_total,
  positivity  = n_pos / n_total
)]

## format for table ------------------------------------------------------------

sesp_table = sesp[, .(
  score_name,
  ca_01,
  n           = format_n(n_total),
  outcome_n   = format_n(n_outcome),
  outcome_pct = sprintf("%.1f%%", 100 * n_outcome / n_total),
  sens        = sprintf("%.1f%%", 100 * sensitivity),
  spec        = sprintf("%.1f%%", 100 * specificity),
  ppv         = sprintf("%.1f%%", 100 * ppv),
  npv         = sprintf("%.1f%%", 100 * npv),
  positive_n  = format_n(n_pos),
  positive_pct = sprintf("%.1f%%", 100 * positivity)
)]

sesp_table[, `:=`(
  ca_lab    = fifelse(ca_01 == 1, "Cancer", "Non-cancer"),
  score_lab = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
)]

setorder(sesp_table, score_name, ca_01)

message("  Sensitivity/specificity calculated for ", nrow(sesp_table), " score × cancer combinations")

# EVER POSITIVE ANALYSIS -------------------------------------------------------

message("\n== Ever positive analysis ==")

## aggregate across sites ------------------------------------------------------

ever_pos = ever_positive_raw[, .(n = sum(n)), by = .(score_name, ca_01, ever_positive, deadhospice_01)]

## calculate rates by group ----------------------------------------------------

ever_pos_summary = ever_pos[, .(
  n_total        = sum(n),
  n_ever_pos     = sum(n[ever_positive == 1]),
  n_with_outcome = sum(n[deadhospice_01 == 1]),
  n_pos_outcome  = sum(n[ever_positive == 1 & deadhospice_01 == 1])
), by = .(score_name, ca_01)]

ever_pos_summary[, `:=`(
  pct_ever_pos      = 100 * n_ever_pos / n_total,
  pct_outcome       = 100 * n_with_outcome / n_total,
  pct_pos_if_event  = 100 * n_pos_outcome / n_with_outcome
)]

ever_pos_summary[, `:=`(
  ca_lab    = fifelse(ca_01 == 1, "Cancer", "Non-cancer"),
  score_lab = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
)]

message("  Ever positive rates calculated")

# TIME TO FIRST POSITIVITY -----------------------------------------------------

message("\n== Time to first positivity ==")

## aggregate across sites ------------------------------------------------------

first_pos = first_pos_raw[, .(
  n                  = sum(n),
  h_from_admit_sum   = sum(h_from_admit_sum),
  h_from_admit_sumsq = sum(h_from_admit_sumsq),
  n_pos_day0         = sum(n_pos_day0),
  n_pos_day1         = sum(n_pos_day1)
), by = .(score, ca_01, o_primary_01)]

## calculate mean and SD -------------------------------------------------------

first_pos[, `:=`(
  mean_h = h_from_admit_sum / n,
  sd_h   = calculate_sd_from_sums(h_from_admit_sum, h_from_admit_sumsq, n)
)]

## format for presentation -----------------------------------------------------

first_pos_summary = first_pos[, .(
  score,
  ca_01,
  outcome      = o_primary_01,
  n            = format_n(n),
  mean_hours   = round(mean_h, 1),
  sd_hours     = round(sd_h, 1),
  median_days  = round(mean_h / 24, 1),
  pct_day0     = sprintf("%.1f%%", 100 * n_pos_day0 / n),
  pct_day0_1   = sprintf("%.1f%%", 100 * (n_pos_day0 + n_pos_day1) / n)
)]

first_pos_summary[, `:=`(
  ca_lab    = fifelse(ca_01 == 1, "Cancer", "Non-cancer"),
  score_lab = factor(score, levels = c("sirs", "qsofa", "mews", "news", "mews_sf"),
                     labels = c("SIRS", "qSOFA", "MEWS", "NEWS", "MEWS-SF"))
)]

message("  Time to first positivity calculated")

# CUMULATIVE INCIDENCE OF SCORE POSITIVITY -------------------------------------

message("\n== Cumulative incidence of score positivity ==")

## aggregate across sites ------------------------------------------------------

cuminc = cuminc_raw[, .(
  n_at_risk    = sum(n_at_risk),
  n_became_pos = sum(n_became_pos),
  n_pos_in_bin = sum(n_pos_in_bin)
), by = .(score, ca_01, o_primary_01, time_bin_start)]

cuminc[, cum_inc := n_became_pos / n_at_risk]

cuminc[, `:=`(
  ca_lab      = fifelse(ca_01 == 1, "Cancer", "Non-cancer"),
  outcome_lab = fifelse(o_primary_01 == 1, "With deterioration", "Without deterioration"),
  score_lab   = factor(score, levels = c("sirs", "qsofa", "mews", "news", "mews_sf"),
                       labels = c("SIRS", "qSOFA", "MEWS", "NEWS", "MEWS-SF"))
)]

message("  Cumulative incidence calculated for ", uniqueN(cuminc$time_bin_start), " time bins")

# COMPONENT ANALYSIS (UPSET) ---------------------------------------------------

message("\n== Component analysis ==")

## score combination patterns --------------------------------------------------

upset = upset_raw[, .(n = sum(n)), by = .(ca_01, outcome, sirs, qsofa, mews, news, mews_sf)]
setorder(upset, -n)

## top patterns by cancer status -----------------------------------------------

top_patterns_ca = upset[ca_01 == 1][1:20]
top_patterns_no = upset[ca_01 == 0][1:20]

## component contributions -----------------------------------------------------

upset_comp = upset_comp_raw[, .(
  n      = sum(n),
  n_encs = sum(n_encs)
), by = .(ca_01, o_primary_01, score, component)]

upset_comp[, pct := 100 * n / n_encs]

upset_comp[, `:=`(
  ca_lab      = fifelse(ca_01 == 1, "Cancer", "Non-cancer"),
  outcome_lab = fifelse(o_primary_01 == 1, "Deteriorated", "No deterioration")
)]

message("  Component patterns analyzed")

# SUMMARY TABLES ---------------------------------------------------------------

message("\n== Creating summary tables ==")

## Table: Sensitivity/Specificity by score and cancer status -------------------

threshold_table = sesp_table[, .(
  score_lab, ca_lab, n, outcome_pct, positive_pct, sens, spec, ppv, npv
)]

threshold_table_wide = dcast(
  threshold_table,
  score_lab ~ ca_lab,
  value.var = c("sens", "spec", "ppv", "npv", "positive_pct")
)

## Table: Time to positivity ---------------------------------------------------

time_pos_table = first_pos_summary[outcome == 1, .(
  score_lab, ca_lab, n, mean_hours, sd_hours, pct_day0, pct_day0_1
)]

# exports ----------------------------------------------------------------------

message("\n== Threshold analysis complete ==")

# Objects for figures
sesp_final         = sesp
sesp_table_final   = sesp_table
ever_pos_final     = ever_pos_summary
first_pos_final    = first_pos_summary
cuminc_final       = cuminc
upset_final        = upset
upset_comp_final   = upset_comp
threshold_table_final = threshold_table
