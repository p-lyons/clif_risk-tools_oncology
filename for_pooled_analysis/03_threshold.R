# ==============================================================================
# 03_threshold.R
# Threshold-based analyses: sensitivity/specificity, time to positivity,
# cumulative incidence, component analysis (main cohort only)
#
# Round-two change: OUTCOME is a parameter. The site pipeline emits sesp, ever,
# first, cuminc, upset, and upset components once per outcome; round one saw
# only the composite, because legacy_view() filters to outcome_key ==
# "composite" before anything downstream reads the data.
#
# Every block below now aggregates from the all-outcome staging object and
# carries outcome_key through, and each round-one object name is re-bound to
# its composite slice at the bottom so 06_figures.R is untouched. Per-outcome
# results travel under *_by_outcome names.
#
# Discharges the threshold half of R27 (do the operating characteristics differ
# once hospice discharge is removed from the composite) and supplies the
# component-outcome operating characteristics R28 promises.
# ==============================================================================

# setup ------------------------------------------------------------------------

library(ggplot2)

## which outcomes to run -------------------------------------------------------
# Composite first so it is always present in the pooled objects even if a later
# outcome is missing from a partial site collection.

THRESHOLD_OUTCOMES = c("composite", "nohospice", "wardicu", "warddeath", "hospicedc")

#' Restrict an all-outcome staging object to the outcomes actually present.
#' Returns an empty table (not an error) when the object is absent, so a partial
#' collection still runs.
outcome_slice = function(nm, label) {

  if (!exists(nm)) {
    message("  ", label, ": ", nm, " not found; skipping")
    return(data.table())
  }

  d = as.data.table(get(nm))

  if (nrow(d) == 0) return(data.table())

  if (!("outcome_key" %in% names(d))) {
    message("  ", label, ": no outcome_key column; treating as composite")
    d[, outcome_key := "composite"]
  }

  keep = intersect(THRESHOLD_OUTCOMES, unique(d$outcome_key))
  miss = setdiff(THRESHOLD_OUTCOMES, keep)

  if (length(miss) > 0) {
    message("  ", label, ": no rows for ", paste(miss, collapse = ", "))
  }

  d[outcome_key %in% keep]
}

#' Attach the outcome display label used by every table below.
add_outcome_lab = function(dt) {
  if (nrow(dt) == 0) return(dt)
  dt[, outcome_lab := factor(outcome_key,
                             levels = names(OUTCOME_LABS),
                             labels = OUTCOME_LABS)]
  dt[]
}

# SENSITIVITY / SPECIFICITY AT STANDARD THRESHOLDS -----------------------------

message("\n== Sensitivity/Specificity at standard thresholds ==")

## aggregate across sites, within outcome --------------------------------------

sesp_src = outcome_slice("sesp_all_raw", "sesp")

sesp = sesp_src[, .(
  n_total      = sum(n_total),
  n_outcome    = sum(n_outcome),
  n_no_outcome = sum(n_no_outcome),
  n_pos        = sum(n_pos),
  n_neg        = sum(n_neg),
  tp           = sum(tp),
  fp           = sum(fp),
  tn           = sum(tn),
  fn           = sum(fn)
), by = .(score_name, ca_01, outcome_key)]

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
  outcome_key,
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

add_outcome_lab(sesp)
add_outcome_lab(sesp_table)

setorder(sesp_table, outcome_key, score_name, ca_01)

message("  Sensitivity/specificity calculated for ", nrow(sesp_table),
        " score × cohort × outcome combinations across ",
        uniqueN(sesp_table$outcome_key), " outcome(s)")

# EVER POSITIVE ANALYSIS -------------------------------------------------------

message("\n== Ever positive analysis ==")

## aggregate across sites ------------------------------------------------------

ever_src = outcome_slice("ever_all_raw", "ever positive")

if (nrow(ever_src) > 0 && "o_out" %in% names(ever_src)) {
  setnames(ever_src, "o_out", "o_primary_01")
}

ever_pos = ever_src[, .(n = sum(n)),
                    by = .(score_name, ca_01, outcome_key, ever_positive, o_primary_01)]

## calculate rates by group ----------------------------------------------------

ever_pos_summary = ever_pos[, .(
  n_total        = sum(n),
  n_ever_pos     = sum(n[ever_positive == 1]),
  n_with_outcome = sum(n[o_primary_01 == 1]),
  n_pos_outcome  = sum(n[ever_positive == 1 & o_primary_01 == 1])
), by = .(score_name, ca_01, outcome_key)]

ever_pos_summary[, `:=`(
  pct_ever_pos      = 100 * n_ever_pos / n_total,
  pct_outcome       = 100 * n_with_outcome / n_total,
  pct_pos_if_event  = 100 * n_pos_outcome / n_with_outcome
)]

ever_pos_summary[, `:=`(
  ca_lab    = fifelse(ca_01 == 1, "Cancer", "Non-cancer"),
  score_lab = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
)]

add_outcome_lab(ever_pos_summary)

message("  Ever positive rates calculated for ",
        uniqueN(ever_pos_summary$outcome_key), " outcome(s)")

# TIME TO FIRST POSITIVITY -----------------------------------------------------

message("\n== Time to first positivity ==")

## aggregate across sites ------------------------------------------------------

first_src = outcome_slice("first_all_raw", "first positivity")

if (nrow(first_src) > 0 && "o_out" %in% names(first_src)) {
  setnames(first_src, "o_out", "o_primary_01")
}

first_pos = first_src[, .(
  n                  = sum(n),
  h_from_admit_sum   = sum(h_from_admit_sum),
  h_from_admit_sumsq = sum(h_from_admit_sumsq),
  n_pos_day0         = sum(n_pos_day0),
  n_pos_day1         = sum(n_pos_day1)
), by = .(score, ca_01, outcome_key, o_primary_01)]

## calculate mean and SD -------------------------------------------------------

first_pos[, `:=`(
  mean_h = h_from_admit_sum / n,
  sd_h   = calculate_sd_from_sums(h_from_admit_sum, h_from_admit_sumsq, n)
)]

## format for presentation -----------------------------------------------------

first_pos_summary = first_pos[, .(
  score,
  ca_01,
  outcome_key,
  outcome      = o_primary_01,
  n            = format_n(n),
  mean_hours   = round(mean_h, 1),
  sd_hours     = round(sd_h, 1),
  mean_days    = round(mean_h / 24, 1),
  pct_day0     = sprintf("%.1f%%", 100 * n_pos_day0 / n),
  pct_day0_1   = sprintf("%.1f%%", 100 * (n_pos_day0 + n_pos_day1) / n)
)]

first_pos_summary[, `:=`(
  ca_lab    = fifelse(ca_01 == 1, "Cancer", "Non-cancer"),
  score_lab = factor(score, levels = c("sirs", "qsofa", "mews", "news", "mews_sf"),
                     labels = c("SIRS", "qSOFA", "MEWS", "NEWS", "MEWS-SF"))
)]

add_outcome_lab(first_pos_summary)

message("  Time to first positivity calculated for ",
        uniqueN(first_pos_summary$outcome_key), " outcome(s)")

# CUMULATIVE INCIDENCE OF SCORE POSITIVITY -------------------------------------

message("\n== Cumulative incidence of score positivity ==")

## aggregate across sites ------------------------------------------------------

## Round-two note: the site artifact (03a_artifacts.R, cuminc_summary) emits
## n_at_risk, n_became_pos, and a per-site cum_inc only. The round-one column
## n_pos_in_bin no longer exists, so it is not aggregated here; cum_inc is
## recomputed from the pooled counts immediately below.

cuminc_src = outcome_slice("cuminc_all_raw", "cumulative incidence")

if (nrow(cuminc_src) > 0 && "o_out" %in% names(cuminc_src)) {
  setnames(cuminc_src, "o_out", "o_primary_01")
}

cuminc = cuminc_src[, .(
  n_at_risk    = sum(n_at_risk),
  n_became_pos = sum(n_became_pos)
), by = .(score, ca_01, outcome_key, o_primary_01, time_bin_start)]

cuminc[, cum_inc := n_became_pos / n_at_risk]

# outcome_lab here names the EVENT-STATUS split within a panel, not the outcome
# definition; the outcome definition is carried by outcome_key and labelled as
# outcome_def_lab so the two cannot be confused downstream.
cuminc[, `:=`(
  ca_lab      = fifelse(ca_01 == 1, "Cancer", "Non-cancer"),
  outcome_lab = fifelse(o_primary_01 == 1, "With deterioration", "Without deterioration"),
  outcome_def_lab = factor(outcome_key, levels = names(OUTCOME_LABS), labels = OUTCOME_LABS),
  score_lab   = factor(score, levels = c("sirs", "qsofa", "mews", "news", "mews_sf"),
                       labels = c("SIRS", "qSOFA", "MEWS", "NEWS", "MEWS-SF"))
)]

message("  Cumulative incidence calculated for ", uniqueN(cuminc$time_bin_start),
        " time bins across ", uniqueN(cuminc$outcome_key), " outcome(s)")

# COMPONENT ANALYSIS (UPSET) ---------------------------------------------------

message("\n== Component analysis ==")

## score combination patterns --------------------------------------------------

upset_src = outcome_slice("upset_all_raw", "co-positivity")

upset = upset_src[, .(n = sum(n)),
                  by = .(ca_01, outcome_key, outcome, sirs, qsofa, mews, news, mews_sf)]
setorder(upset, outcome_key, -n)

add_outcome_lab(upset)

## top patterns by cancer status (composite outcome) ---------------------------
# head() rather than [1:20], which pads with NA rows when fewer patterns exist.

top_patterns_ca = head(upset[ca_01 == 1 & outcome_key == "composite"], 20L)
top_patterns_no = head(upset[ca_01 == 0 & outcome_key == "composite"], 20L)

## component contributions -----------------------------------------------------

upset_comp_src = outcome_slice("upset_comp_all_raw", "components")

if (nrow(upset_comp_src) > 0 && "o_out" %in% names(upset_comp_src)) {
  setnames(upset_comp_src, "o_out", "o_primary_01")
}

upset_comp = upset_comp_src[, .(
  n      = sum(n),
  n_encs = sum(n_encs)
), by = .(ca_01, outcome_key, o_primary_01, score, component)]

upset_comp[, pct := 100 * n / n_encs]

upset_comp[, `:=`(
  ca_lab          = fifelse(ca_01 == 1, "Cancer", "Non-cancer"),
  outcome_lab     = fifelse(o_primary_01 == 1, "Deteriorated", "No deterioration"),
  outcome_def_lab = factor(outcome_key, levels = names(OUTCOME_LABS), labels = OUTCOME_LABS)
)]

comp_pooled =
  upset_comp_src[outcome_key == "composite"] |>
  fgroup_by(ca_01, score, component) |>
  fsummarize(n = fsum(n), n_encs = fsum(n_encs)) |>
  fmutate(pct = round(100 * n / n_encs, 1)) |>
  roworder(score, ca_01, -pct)

# differences and rank comparisons
comp_diff =
  comp_pooled |>
  pivot_wider(
    names_from  = ca_01,
    values_from = c(n, n_encs, pct),
    names_sep   = "_"
  ) |>
  fmutate(
    diff     = pct_1 - pct_0,
    abs_diff = abs(pct_1 - pct_0)
  ) |>
  roworder(score, -abs_diff)

# rank comparison
comp_ranks =
  comp_pooled |>
  fgroup_by(ca_01, score) |>
  fmutate(rank = frank(-pct)) |>
  fselect(ca_01, score, component, rank) |>
  pivot_wider(
    names_from  = ca_01,
    values_from = rank,
    names_prefix = "rank_"
  ) |>
  fmutate(rank_shift = rank_1 - rank_0)

# join them
comp_summary = join(comp_diff, comp_ranks, on = c("score", "component"))

# print nicely
comp_summary |>
  fselect(score, component, pct_0, pct_1, diff, rank_0, rank_1, rank_shift) |>
  print(n = 50)

message("  Component patterns analyzed")

# SUMMARY TABLES ---------------------------------------------------------------

message("\n== Creating summary tables ==")

## Table: Sensitivity/Specificity by score and cancer status -------------------

threshold_table = sesp_table[outcome_key == "composite", .(
  score_lab, ca_lab, n, outcome_pct, positive_pct, sens, spec, ppv, npv
)]

threshold_table_wide = dcast(
  threshold_table,
  score_lab ~ ca_lab,
  value.var = c("sens", "spec", "ppv", "npv", "positive_pct")
)

## Table: standard-threshold operating characteristics by outcome --------------
# R27 asks whether the cancer gap in operating characteristics survives the
# removal of hospice discharge; R28 asks for the components separately. Both
# read off this table.

threshold_by_outcome = sesp_table[, .(
  Outcome       = outcome_lab,
  Score         = score_lab,
  Cohort        = ca_lab,
  N             = n,
  `Outcome rate` = outcome_pct,
  `Alert rate`  = positive_pct,
  Sensitivity   = sens,
  Specificity   = spec,
  PPV           = ppv,
  NPV           = npv
)]

setorder(threshold_by_outcome, Outcome, Score, Cohort)

if (requireNamespace("flextable", quietly = TRUE) && nrow(threshold_by_outcome) > 0) {

  thr_out_path = if (exists("today")) {
    here("output", "tables", paste0("table_threshold_by_outcome_", today, ".docx"))
  } else {
    here("output", "tables", "table_threshold_by_outcome.docx")
  }

  ft_thr_out = flextable::flextable(as.data.frame(threshold_by_outcome)) |>
    flextable::autofit() |>
    flextable::bold(part = "header")

  flextable::save_as_docx(ft_thr_out, path = thr_out_path)
  message("  Saved: ", thr_out_path, " (", nrow(threshold_by_outcome), " rows)")
}

## Table: Time to positivity ---------------------------------------------------

time_pos_table = first_pos_summary[outcome == 1 & outcome_key == "composite", .(
  score_lab, ca_lab, n, mean_hours, sd_hours, pct_day0, pct_day0_1
)]

# exports ----------------------------------------------------------------------

message("\n== Threshold analysis complete ==")

# Note: COHORT_N, VARIANT_N, and SITE_N are available from 00_load.R
# These should be used for figure labels rather than calculating from data
message("  Using COHORT_N for figure labels: ", 
        format_n(COHORT_N["0"]), " non-cancer, ",
        format_n(COHORT_N["1"]), " cancer")

## composite slice under the round-one names ----------------------------------
# 06_figures.R consumes these and assumes one outcome definition per row.
# Handing it the full per-outcome table would multiply every panel, so the
# legacy names stay composite-only.

sesp_final            = sesp[outcome_key == "composite"]
sesp_table_final      = sesp_table[outcome_key == "composite"]
ever_pos_final        = ever_pos_summary[outcome_key == "composite"]
first_pos_final       = first_pos_summary[outcome_key == "composite"]
cuminc_final          = cuminc[outcome_key == "composite"]
upset_final           = upset[outcome_key == "composite"]
upset_comp_final      = upset_comp[outcome_key == "composite"]
threshold_table_final = threshold_table

## per-outcome objects ---------------------------------------------------------

sesp_by_outcome            = sesp
sesp_table_by_outcome      = sesp_table
ever_pos_by_outcome        = ever_pos_summary
first_pos_by_outcome       = first_pos_summary
cuminc_by_outcome          = cuminc
upset_by_outcome           = upset
upset_comp_by_outcome      = upset_comp
threshold_by_outcome_final = threshold_by_outcome

message("  Composite objects bound under their round-one names; ",
        uniqueN(sesp$outcome_key), " outcome(s) available in the *_by_outcome objects")
