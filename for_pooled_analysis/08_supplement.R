# ==============================================================================
# 08_supplement.R
# Pooled analyses for artifact families that 00_load.R reads and that no other
# script consumed. Each block discharges a specific reviewer response.
#
#   S1  first-ward-observation score distributions          R06
#   S2  hospital type and hospital count                    R08, R18
#   S3  admission diagnoses by chapter and code stem        R09
#   S4  number needed to evaluate, lead time, crossing      R25
#   S5  monitoring intensity and component missingness      R30
#   S6  between-site descriptive heterogeneity              R31
#
# Every block guards on its source object, reports what it found, and writes a
# Word table to output/tables. A missing family produces a message and an empty
# object rather than an error, so a partial site collection still runs.
# ==============================================================================

# setup ------------------------------------------------------------------------

library(flextable)
library(officer)
library(here)

message("\n== 08_supplement.R ==")

sup_path = function(stem) {
  if (exists("today")) {
    here("output", "tables", paste0(stem, "_", today, ".docx"))
  } else {
    here("output", "tables", paste0(stem, ".docx"))
  }
}

save_supplement = function(dt, stem, label) {

  if (is.null(dt) || nrow(dt) == 0) {
    message("  ", label, ": nothing to write")
    return(invisible(NULL))
  }

  ft = flextable(as.data.frame(dt)) |>
    autofit() |>
    bold(part = "header")

  p = sup_path(stem)
  save_as_docx(ft, path = p)
  message("  Saved: ", p, " (", nrow(dt), " rows)")
  invisible(p)
}

have = function(nm) exists(nm) && nrow(get(nm)) > 0

## weighted summary helpers ----------------------------------------------------
# The site artifacts are count tables, so every distributional summary has to be
# computed from (value, n) pairs rather than from observations.

weighted_mean_from_counts = function(value, n) {
  sum(as.numeric(value) * n) / sum(n)
}

weighted_sd_from_counts = function(value, n) {
  m  = weighted_mean_from_counts(value, n)
  nn = sum(n)
  if (nn <= 1) return(NA_real_)
  sqrt(sum(n * (as.numeric(value) - m)^2) / (nn - 1))
}

weighted_quantile_from_counts = function(value, n, probs) {
  o   = order(value)
  v   = as.numeric(value)[o]
  w   = as.numeric(n)[o]
  cum = cumsum(w) / sum(w)
  vapply(probs, function(p) v[which(cum >= p)[1]], numeric(1))
}

fmt_pct = function(x, digits = 1) {
  fifelse(is.na(x), "—", sprintf(paste0("%.", digits, "f%%"), 100 * x))
}

# ==============================================================================
# S1. FIRST-WARD-OBSERVATION SCORE DISTRIBUTIONS  (R06)
# ==============================================================================
# Reviewer 1 asked for a subanalysis by emergency-department triage category.
# CLIF carries no triage acuity, so the response offers the score at the first
# ward observation as a measure of presenting physiology at the moment the
# monitoring period begins. firstscore_raw holds (score_name, ca_01, value, n).

message("\n-- S1: first-ward-observation score distributions (R06) --")

firstscore_summary = data.table()

if (have("firstscore_raw")) {

  fs = as.data.table(firstscore_raw)[!is.na(value)]

  fs_pooled = fs[, .(n = sum(n)), by = .(score_name, ca_01, value)]

  firstscore_summary = fs_pooled[, {
    q = weighted_quantile_from_counts(value, n, c(0.25, 0.50, 0.75))
    .(
      n_encounters = sum(n),
      mean_score   = weighted_mean_from_counts(value, n),
      sd_score     = weighted_sd_from_counts(value, n),
      median_score = q[2],
      q25          = q[1],
      q75          = q[3]
    )
  }, by = .(score_name, ca_01)]

  firstscore_summary[, `:=`(
    score_lab = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS),
    ca_lab    = fifelse(ca_01 == 1L, "Cancer", "Non-cancer")
  )]

  setorder(firstscore_summary, score_lab, -ca_01)

  firstscore_table = firstscore_summary[, .(
    Score      = score_lab,
    Cohort     = ca_lab,
    N          = format_n(n_encounters),
    `Mean (SD)`     = sprintf("%.2f (%.2f)", mean_score, sd_score),
    `Median (IQR)`  = sprintf("%.0f (%.0f-%.0f)", median_score, q25, q75)
  )]

  save_supplement(firstscore_table, "table_s_firstscore", "First-ward score")

} else {
  message("  firstscore_raw unavailable; skipping")
  firstscore_table = data.table()
}

# ==============================================================================
# S2. HOSPITAL TYPE AND HOSPITAL COUNT  (R08, R18)
# ==============================================================================
# hospital_types_raw carries one row per (site, hospital_type) with n_hospitals,
# resolved at the site by majority vote over the ADT hospital_type field. This
# supplies the hospital-composition columns of eTable 1 and the corrected total
# that R18 concedes. Bed counts, IRB numbers, and date ranges are not pipeline
# outputs and are joined by hand from site_metadata.csv if that file exists.

message("\n-- S2: hospital types and count (R08, R18) --")

hospital_type_table = data.table()

if (have("hospital_types_raw")) {

  ht = as.data.table(hospital_types_raw)[, .(n_hospitals = sum(n_hospitals)),
                                         by = .(site, hospital_type)]

  ht_wide = dcast(ht, site ~ hospital_type, value.var = "n_hospitals", fill = 0L)

  type_cols = setdiff(names(ht_wide), "site")
  ht_wide[, hospitals_total := rowSums(.SD), .SDcols = type_cols]

  total_hospitals = sum(ht_wide$hospitals_total)

  message("  Pooled hospital count: ", total_hospitals, " across ",
          nrow(ht_wide), " health system(s)")
  message("  By type: ",
          paste(type_cols, "=", colSums(ht_wide[, ..type_cols]), collapse = ", "))

  ## optional hand-maintained metadata (beds, IRB, date range) ------------------

  meta_file = here("site_metadata.csv")

  if (file.exists(meta_file)) {
    site_meta = fread(meta_file)
    ht_wide   = merge(ht_wide, site_meta, by = "site", all.x = TRUE)
    message("  Joined site_metadata.csv")
  } else {
    message("  site_metadata.csv not found; bed counts, IRB numbers, and date ",
            "ranges left for manual entry")
    ht_wide[, `:=`(
      adult_inpatient_beds = NA_character_,
      adult_icu_beds       = NA_character_,
      date_range           = NA_character_,
      irb                  = NA_character_
    )]
  }

  ht_wide[, site_lab := toupper(site)]
  setcolorder(ht_wide, c("site_lab", type_cols, "hospitals_total"))
  ht_wide[, site := NULL]

  hospital_type_table = ht_wide

  save_supplement(hospital_type_table, "table_s_hospital_types", "Hospital types")

} else {
  message("  hospital_types_raw unavailable; skipping")
  total_hospitals = NA_integer_
}

# ==============================================================================
# S3. ADMISSION DIAGNOSES  (R09)
# ==============================================================================
# Reviewer 2 asked why patients were hospitalized and whether the reasons differ
# by cancer status. CLIF carries no harmonized chief complaint, so the response
# offers primary hospital diagnoses by ICD-10-CM chapter and by the most
# frequent three-character stems. Site tallies below six encounters are already
# suppressed at export (02_scores.R), so pooled percentages are computed over
# the retained rows and the denominator is stated.

message("\n-- S3: admission diagnoses (R09) --")

adm_chapter_table = data.table()
adm_stem_table    = data.table()

if (have("adm_dx_chapter_raw")) {

  ch = as.data.table(adm_dx_chapter_raw)[, .(n = sum(n)), by = .(ca_01, chapter)]
  ch[, pct := n / sum(n), by = ca_01]

  adm_chapter_table = dcast(
    ch[, .(chapter, ca_01, cell = paste0(format_n(n), " (", fmt_pct(pct), ")"))],
    chapter ~ ca_01,
    value.var = "cell",
    fill = "—"
  )

  setnames(adm_chapter_table, c("0", "1"), c("No Cancer", "Cancer"), skip_absent = TRUE)
  setnames(adm_chapter_table, "chapter", "ICD-10-CM chapter")

  save_supplement(adm_chapter_table, "table_s_admission_dx_chapter", "Admission dx chapter")

} else {
  message("  adm_dx_chapter_raw unavailable; skipping chapter table")
}

if (have("adm_dx_stem_raw")) {

  st = as.data.table(adm_dx_stem_raw)[, .(n = sum(n)), by = .(ca_01, code_stem)]
  st[, pct := n / sum(n), by = ca_01]

  # top 15 stems within each cohort, unioned so both columns are comparable
  top_stems = unique(c(
    st[ca_01 == 1][order(-n)][1:min(15L, .N)]$code_stem,
    st[ca_01 == 0][order(-n)][1:min(15L, .N)]$code_stem
  ))
  top_stems = top_stems[!is.na(top_stems)]

  adm_stem_table = dcast(
    st[code_stem %in% top_stems,
       .(code_stem, ca_01, cell = paste0(format_n(n), " (", fmt_pct(pct), ")"))],
    code_stem ~ ca_01,
    value.var = "cell",
    fill = "—"
  )

  setnames(adm_stem_table, c("0", "1"), c("No Cancer", "Cancer"), skip_absent = TRUE)
  setnames(adm_stem_table, "code_stem", "ICD-10-CM code")

  save_supplement(adm_stem_table, "table_s_admission_dx_stem", "Admission dx stem")

} else {
  message("  adm_dx_stem_raw unavailable; skipping stem table")
}

# ==============================================================================
# S4. NUMBER NEEDED TO EVALUATE, LEAD TIME, AND CROSSING  (R25)
# ==============================================================================
# Reviewer 4 asked for number needed to evaluate and for time from threshold
# crossing to deterioration. Number needed to evaluate is 1/PPV at each score's
# standard threshold and derives from the sesp cells already pooled.
#
# Lead time arrives in two artifacts. leadtime holds cumulative counts at fixed
# thresholds, which pool by summation and give the encounter-weighted "within x
# hours" proportions directly. leadtime_median holds exact medians computed at
# each site for the health system and for each hospital_id; those are NOT
# pooled, and the range across units is the reported statistic. Both carry
# crossing_def: `first` is the conventional first-alert-of-the-encounter
# definition, and `final_onset` is the onset of the last uninterrupted positive
# run, which is the warning time a clinician means. crossclass gives the 2x2 of
# ever-crossed against ever-had-the-event.
#
# No mean is reported. The distribution is bounded below at zero with a heavy
# right tail, and the round-one mean of 182 hours with an SD of 255 described it
# badly. The sufficient statistics are still exported, so a mean can be
# recovered if a reviewer asks.

message("\n-- S4: NNE, lead time, threshold crossing (R25) --")

nne_table = data.table()

if (have("sesp_raw")) {

  sesp_pool = as.data.table(sesp_raw)[, .(
    n_total = sum(n_total),
    tp      = sum(tp),
    fp      = sum(fp),
    tn      = sum(tn),
    fn      = sum(fn),
    n_pos   = sum(n_pos)
  ), by = .(score_name, ca_01)]

  sesp_pool[, `:=`(
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    ppv         = tp / (tp + fp),
    positivity  = n_pos / n_total
  )]

  # Number needed to evaluate: alerts that must be worked up per true event.
  sesp_pool[, nne := fifelse(ppv > 0, 1 / ppv, NA_real_)]

  sesp_pool[, `:=`(
    score_lab = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS),
    ca_lab    = fifelse(ca_01 == 1L, "Cancer", "Non-cancer")
  )]

  setorder(sesp_pool, score_lab, -ca_01)

  nne_table = sesp_pool[, .(
    Score           = score_lab,
    Cohort          = ca_lab,
    `Alert rate`    = fmt_pct(positivity),
    Sensitivity     = fmt_pct(sensitivity),
    PPV             = fmt_pct(ppv),
    `NNE (1/PPV)`   = fifelse(is.na(nne), "—", sprintf("%.1f", nne))
  )]

  save_supplement(nne_table, "table_s_nne", "Number needed to evaluate")

} else {
  message("  sesp_raw unavailable; skipping NNE table")
}

leadtime_table        = data.table()
leadtime_median_table = data.table()

## S4a. cumulative lead-time distribution (poolable) ---------------------------
# leadtime_raw holds n_at_or_below and n_total at fixed thresholds, so pooling is
# summation with no distributional assumption. "Within x" comes straight from the
# ratio; "beyond x" is its complement. Reported alongside the alert rate, because
# a long warning time from a score that flags most encounters is a different
# claim from the same warning time at low positivity.

if (have("leadtime_raw")) {

  lt = as.data.table(leadtime_raw)[, .(
    n_at_or_below = sum(n_at_or_below),
    n_total       = sum(n_total)
  ), by = .(score_name, ca_01, outcome_key, crossing_def, threshold_h)]

  lt[, pct_within := n_at_or_below / n_total]

  lt[, `:=`(
    score_lab   = factor(score_name,  levels = names(SCORE_LABS),   labels = SCORE_LABS),
    ca_lab      = fifelse(ca_01 == 1L, "Cancer", "Non-cancer"),
    outcome_lab = factor(outcome_key, levels = names(OUTCOME_LABS), labels = OUTCOME_LABS),
    def_lab     = fifelse(crossing_def == "first",
                          "First crossing", "Final positive run")
  )]

  setorder(lt, outcome_lab, def_lab, score_lab, -ca_01, threshold_h)

  leadtime_table = dcast(
    lt[, .(outcome_lab, def_lab, score_lab, ca_lab,
           th = paste0("<= ", threshold_h, " h"),
           cell = fmt_pct(pct_within))],
    outcome_lab + def_lab + score_lab + ca_lab ~ factor(th, levels = unique(th)),
    value.var = "cell",
    fill = "—"
  )

  setnames(leadtime_table,
           c("outcome_lab", "def_lab", "score_lab", "ca_lab"),
           c("Outcome", "Crossing", "Score", "Cohort"))

  save_supplement(leadtime_table, "table_s_leadtime", "Lead-time distribution")

} else {
  message("  leadtime_raw unavailable; skipping lead-time distribution table")
}

## S4b. per-unit medians, reported as a range ----------------------------------
# Each median is exact, computed on line-level data at its own site. Sites are
# not pooled; the reported statistic is the range across units. Hospital-level
# units are filtered by LEADTIME_EVENT_FLOOR first, because the range is the
# statistic most sensitive to the smallest contributors and some health systems
# report nine or more hospitals.

LEADTIME_EVENT_FLOOR = 20L

if (have("leadtime_median_raw")) {

  lm = as.data.table(leadtime_median_raw)

  n_dropped = lm[unit == "hospital" & n_events < LEADTIME_EVENT_FLOOR, .N]
  if (n_dropped > 0L) {
    message("  Dropped ", n_dropped, " hospital row(s) below ",
            LEADTIME_EVENT_FLOOR, " events")
  }

  lm = lm[unit == "health_system" | n_events >= LEADTIME_EVENT_FLOOR]

  lm_range = lm[, .(
    k_units  = .N,
    n_events = sum(n_events),
    md_lo    = min(median_h),
    md_md    = median(median_h),
    md_hi    = max(median_h)
  ), by = .(score_name, ca_01, outcome_key, crossing_def, unit)]

  lm_range[, `:=`(
    score_lab   = factor(score_name,  levels = names(SCORE_LABS),   labels = SCORE_LABS),
    ca_lab      = fifelse(ca_01 == 1L, "Cancer", "Non-cancer"),
    outcome_lab = factor(outcome_key, levels = names(OUTCOME_LABS), labels = OUTCOME_LABS),
    def_lab     = fifelse(crossing_def == "first",
                          "First crossing", "Final positive run"),
    unit_lab    = fifelse(unit == "health_system", "Health system", "Hospital")
  )]

  setorder(lm_range, outcome_lab, def_lab, unit_lab, score_lab, -ca_01)

  leadtime_median_table = lm_range[, .(
    Outcome                     = outcome_lab,
    Crossing                    = def_lab,
    Unit                        = unit_lab,
    Score                       = score_lab,
    Cohort                      = ca_lab,
    Units                       = k_units,
    Events                      = format_n(n_events),
    `Median of unit medians, h` = sprintf("%.1f", md_md),
    `Range across units, h`     = sprintf("%.1f-%.1f", md_lo, md_hi)
  )]

  save_supplement(leadtime_median_table, "table_s_leadtime_median", "Lead-time medians")

} else {
  message("  leadtime_median_raw unavailable; skipping lead-time median table")
}

crossclass_table = data.table()

if (have("crossclass_raw")) {

  cc = as.data.table(crossclass_raw)[, .(n = sum(n)),
                                     by = .(score_name, ca_01, outcome_key, crossed, event)]

  cc_wide = cc[, .(
    n_total       = sum(n),
    crossed_event = sum(n[crossed == 1L & event == 1L]),
    crossed_no    = sum(n[crossed == 1L & event == 0L]),
    nocross_event = sum(n[crossed == 0L & event == 1L]),
    nocross_no    = sum(n[crossed == 0L & event == 0L])
  ), by = .(score_name, ca_01, outcome_key)]

  cc_wide[, `:=`(
    pct_crossed        = (crossed_event + crossed_no) / n_total,
    pct_event_crossed  = crossed_event / (crossed_event + nocross_event),
    pct_missed         = nocross_event / (crossed_event + nocross_event)
  )]

  cc_wide[, `:=`(
    score_lab   = factor(score_name,  levels = names(SCORE_LABS),   labels = SCORE_LABS),
    ca_lab      = fifelse(ca_01 == 1L, "Cancer", "Non-cancer"),
    outcome_lab = factor(outcome_key, levels = names(OUTCOME_LABS), labels = OUTCOME_LABS)
  )]

  setorder(cc_wide, outcome_lab, score_lab, -ca_01)

  crossclass_table = cc_wide[, .(
    Outcome                      = outcome_lab,
    Score                        = score_lab,
    Cohort                       = ca_lab,
    `Ever crossed threshold`     = fmt_pct(pct_crossed),
    `Events preceded by crossing` = fmt_pct(pct_event_crossed),
    `Events never flagged`       = fmt_pct(pct_missed)
  )]

  save_supplement(crossclass_table, "table_s_crossclass", "Threshold crossing")

} else {
  message("  crossclass_raw unavailable; skipping crossing table")
}

# ==============================================================================
# S5. MONITORING INTENSITY AND COMPONENT MISSINGNESS  (R30)
# ==============================================================================
# Reviewer 4 observed that carry-forward plus treat-missing-as-normal interacts
# with monitoring intensity, which may differ by cancer status. The monitoring
# artifact reports, per encounter and measure, the number of distinct
# measurement times within the at-risk ward window and the implied rate per 24
# hours; missing_vlab reports the share of encounters with no value at all.

message("\n-- S5: monitoring intensity and missingness (R30) --")

monitoring_table = data.table()

if (have("monitoring_raw")) {

  mon = as.data.table(monitoring_raw)[, .(
    n_enc         = sum(n_enc),
    sum_obs       = sum(sum_obs),
    sumsq_obs     = sum(sumsq_obs),
    sum_ward_h    = sum(sum_ward_hours),
    sum_rate      = sum(sum_rate_per24h),
    sumsq_rate    = sum(sumsq_rate_per24h)
  ), by = .(ca_01, measure)]

  mon[, `:=`(
    mean_obs   = sum_obs  / n_enc,
    sd_obs     = calculate_sd_from_sums(sum_obs,  sumsq_obs,  n_enc),
    mean_rate  = sum_rate / n_enc,
    sd_rate    = calculate_sd_from_sums(sum_rate, sumsq_rate, n_enc)
  )]

  mon[, ca_lab := fifelse(ca_01 == 1L, "Cancer", "Non-cancer")]
  setorder(mon, measure, -ca_01)

  monitoring_table = mon[, .(
    Measure                       = measure,
    Cohort                        = ca_lab,
    Encounters                    = format_n(n_enc),
    `Measurements, mean (SD)`     = sprintf("%.1f (%.1f)", mean_obs,  sd_obs),
    `Per 24 h, mean (SD)`         = sprintf("%.2f (%.2f)", mean_rate, sd_rate)
  )]

  save_supplement(monitoring_table, "table_s_monitoring", "Monitoring intensity")

} else {
  message("  monitoring_raw unavailable; skipping monitoring table")
}

monitoring_bins_table = data.table()

if (have("monitoring_bins_raw")) {

  mb = as.data.table(monitoring_bins_raw)[, .(n_enc = sum(n_enc)),
                                          by = .(ca_01, measure, rate_bin)]
  mb[, pct := n_enc / sum(n_enc), by = .(ca_01, measure)]

  monitoring_bins_table = dcast(
    mb[, .(measure, rate_bin, ca_01, cell = fmt_pct(pct))],
    measure + rate_bin ~ ca_01,
    value.var = "cell",
    fill = "—"
  )

  setnames(monitoring_bins_table, c("0", "1"), c("No Cancer", "Cancer"), skip_absent = TRUE)
  setnames(monitoring_bins_table, c("measure", "rate_bin"),
           c("Measure", "Measurements per 24 h"))

  save_supplement(monitoring_bins_table, "table_s_monitoring_bins", "Monitoring bins")

} else {
  message("  monitoring_bins_raw unavailable; skipping monitoring-bin table")
}

missingness_table = data.table()

if (have("missing_vlab_raw")) {

  mv = as.data.table(missing_vlab_raw)[, .(
    n_total   = sum(n_total),
    n_missing = sum(n_missing)
  ), by = .(variable, ca_01)]

  mv[, pct_missing := n_missing / n_total]
  mv[, ca_lab := fifelse(ca_01 == 1L, "Cancer", "Non-cancer")]

  missingness_table = dcast(
    mv[, .(variable, ca_lab,
           cell = paste0(format_n(n_missing), " (", fmt_pct(pct_missing), ")"))],
    variable ~ ca_lab,
    value.var = "cell",
    fill = "—"
  )

  setnames(missingness_table, "variable", "Component")

  save_supplement(missingness_table, "table_s_missingness", "Component missingness")

} else {
  message("  missing_vlab_raw unavailable; skipping missingness table")
}

# ==============================================================================
# S6. BETWEEN-SITE DESCRIPTIVE HETEROGENEITY  (R31)
# ==============================================================================
# R31 promises, alongside I-squared and tau-squared (produced in
# 02_discrimination.R), the range across sites of the primary outcome rate and
# of standard-threshold alert positivity. Both come from the per-site sesp
# cells: the outcome rate is invariant across scores within a site, so one score
# supplies it.

message("\n-- S6: between-site descriptive heterogeneity (R31) --")

site_range_table = data.table()

if (have("sesp_raw")) {

  ss = as.data.table(sesp_raw)

  ## outcome rate by site (one score as reference) ------------------------------

  ref_score = ss$score_name[1]

  rate_by_site = ss[score_name == ref_score, .(
    outcome_rate = sum(n_outcome) / sum(n_total)
  ), by = .(site, ca_01)]

  rate_range = rate_by_site[, .(
    metric = "Primary outcome rate",
    score  = "—",
    k      = .N,
    lo     = min(outcome_rate),
    md     = median(outcome_rate),
    hi     = max(outcome_rate)
  ), by = ca_01]

  ## alert positivity by site and score ----------------------------------------

  pos_by_site = ss[, .(positivity = sum(n_pos) / sum(n_total)),
                   by = .(site, ca_01, score_name)]

  pos_range = pos_by_site[, .(
    metric = "Alert positivity at standard threshold",
    k      = .N,
    lo     = min(positivity),
    md     = median(positivity),
    hi     = max(positivity)
  ), by = .(ca_01, score_name)]

  setnames(pos_range, "score_name", "score")

  site_range = rbindlist(list(rate_range, pos_range), use.names = TRUE, fill = TRUE)

  site_range[, `:=`(
    ca_lab    = fifelse(ca_01 == 1L, "Cancer", "Non-cancer"),
    score_lab = fifelse(score == "—", "—",
                        as.character(factor(score,
                                            levels = names(SCORE_LABS),
                                            labels = SCORE_LABS)))
  )]

  setorder(site_range, metric, score_lab, -ca_01)

  site_range_table = site_range[, .(
    Measure               = metric,
    Score                 = score_lab,
    Cohort                = ca_lab,
    Sites                 = k,
    `Median (range)`      = paste0(fmt_pct(md), " (", fmt_pct(lo), "-", fmt_pct(hi), ")")
  )]

  save_supplement(site_range_table, "table_s_site_range", "Between-site ranges")

} else {
  message("  sesp_raw unavailable; skipping between-site range table")
}

# ==============================================================================
# EXPORTS
# ==============================================================================

message("\n== 08_supplement.R complete ==")

firstscore_table_final      = firstscore_table
hospital_type_table_final   = hospital_type_table
adm_chapter_table_final     = adm_chapter_table
adm_stem_table_final        = adm_stem_table
nne_table_final             = nne_table
leadtime_table_final        = leadtime_table
leadtime_median_table_final = leadtime_median_table
crossclass_table_final      = crossclass_table
monitoring_table_final      = monitoring_table
monitoring_bins_table_final = monitoring_bins_table
missingness_table_final     = missingness_table
site_range_table_final      = site_range_table

gc()
