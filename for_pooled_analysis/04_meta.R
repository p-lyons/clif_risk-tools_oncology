# ==============================================================================
# 04_threshold_equivalence.R
#
# Replaces the former 04_meta.R (score x cancer interaction meta-analysis).
# The interaction approach produced a per-SD / per-point OR that (a) was
# guaranteed to be "significant" at our N regardless of clinical meaning, and
# (b) did not directly answer the operational question: if we simply move the
# threshold for each score in cancer patients, can we recover equivalent
# performance?
#
# Three complementary analyses, all built from existing federated count tables
# (maxscores_ca_raw, sesp_raw) without requiring new site-level extractions:
#
#   1. Risk-difference table: observed deterioration rate (cancer - noncancer)
#      at each integer score value, with Newcombe CIs on the difference and a
#      weighted-regression test for trend in the risk differences across score
#      values (parallel-shift test).
#
#   2. Threshold-equivalence table: at the standard threshold in noncancer,
#      report sens/spec/PPV/alert rate; then find the integer threshold in
#      cancer that matches on sensitivity, and separately on alert rate,
#      reporting the resulting tradeoffs.
#
#   3. Efficiency curves: sensitivity (x) vs proportion alerted (y), one curve
#      per cohort per score, built by sweeping the integer threshold. Based on
#      Romero-Brufau et al. (JAMIA Open 2020; PMC7237982).
#
# Flextables for the threshold-equivalence analysis are constructed at the
# end of this script and written to output/tables/:
#   - table_threshold_eq_wide.docx   : all brackets, grouped headers (suppl)
#   - table_threshold_eq_sens.docx   : sens-matched only (main text candidate)
#   - table_threshold_eq_alert.docx  : alert-rate-matched only
#
# Primary analyses use crude pooled counts (summed across sites) because the
# clinical question is operational ("what sens/alert-rate would these 8 sites
# see at this threshold"), not inferential ("what would a hypothetical new
# site see"). Random-effects meta-analysis would inflate uncertainty for this
# descriptive question. However, pooled counts are vulnerable to Simpson's
# paradox if site-level cancer mix and baseline deterioration rates are
# heterogeneous. To address this, we also compute the threshold sweep
# PER SITE (site_sweep_final) and summarize between-site variation at the
# standard threshold (het_diag_final). If any score shows max deviation from
# the pooled estimate exceeding HET_FLAG_THRESHOLD (10 pp), it is flagged,
# and we will present a meta-analyzed version for that score as a sensitivity
# check in the supplement.
#
# Wilson CIs for proportions; Newcombe method 10 for risk differences.
#
# ----- Positivity conventions -----
#
# For NEWS, the field's "any single parameter scores 3" rule triggers
# positivity in addition to aggregate NEWS >= 5. The site-side sesp artifact
# incorporates both rules. The maxscores artifact, however, captures only the
# aggregate news_total score (the news_any3 flag is stored separately and not
# aggregated into the count table). Because the threshold-equivalence and
# efficiency-curve analyses depend on sweeping an integer threshold across
# the score, and the any-3 rule cannot be varied as a threshold, these
# analyses use the aggregate-score-only rule for all five scores including
# NEWS. The per-site supplementary check (site_ops) which is sourced from
# sesp_raw therefore differs from the main sweep for NEWS at the standard
# threshold; we report this explicitly rather than trying to paper over it.
# The standard-positivity numbers quoted in the main methods (e.g., Table 2)
# continue to use the field's conventional NEWS rule (aggregate-plus-any3)
# and come from sesp_raw directly.
#
# ----- Note for future federated extraction -----
#
# The primary outcome is the composite of ward-ICU transfer, ward death, and
# hospice discharge. A reviewer may reasonably ask whether threshold-
# equivalence results change when hospice is excluded, since hospice
# discharge often reflects goals-of-care decisions rather than physiologic
# deterioration that a vital-signs score is designed to detect. If component
# outcomes are exported separately in a future run, we can rerun these
# analyses for the ICU-transfer-or-death subcomposite as a planned
# sensitivity analysis.
# ==============================================================================

# setup ------------------------------------------------------------------------

library(ggplot2)

## standard thresholds for aggregate-score positivity -------------------------
# Applied to cleaned score_name values (post "_total" strip in 00_load.R).

STD_THRESHOLDS = c(
  sirs    = 2L,
  qsofa   = 2L,
  mews    = 5L,
  mews_sf = 7L,
  news    = 5L
)

## minimum cell size for presentation-table inclusion -------------------------
# Keep all rows in underlying _final objects; apply this only to the
# formatted rd_table to avoid tiny-N cells dominating visual scan.

MIN_N_PER_COHORT = 50L

# HELPERS ---------------------------------------------------------------------

## Wilson score interval for a proportion (vectorized) ------------------------

wilson_ci = function(x, n, conf = 0.95) {
  z = qnorm(1 - (1 - conf) / 2)
  p = fifelse(n > 0, x / n, NA_real_)
  denom = 1 + z^2 / n
  center = (p + z^2 / (2 * n)) / denom
  halfw  = (z / denom) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  list(
    lo = pmax(0, center - halfw),
    hi = pmin(1, center + halfw)
  )
}

## Newcombe method 10 CI for a risk difference (vectorized) -------------------
# RD = p1 - p2; caller passes cancer first to get (cancer - noncancer).
# Reference: Newcombe RG. Stat Med 1998;17:873-890.

newcombe_rd_ci = function(x1, n1, x2, n2, conf = 0.95) {
  ci1 = wilson_ci(x1, n1, conf)
  ci2 = wilson_ci(x2, n2, conf)
  p1 = fifelse(n1 > 0, x1 / n1, NA_real_)
  p2 = fifelse(n2 > 0, x2 / n2, NA_real_)
  delta = sqrt((p1 - ci1$lo)^2 + (ci2$hi - p2)^2)
  epsil = sqrt((ci1$hi - p1)^2 + (p2 - ci2$lo)^2)
  list(
    rd    = p1 - p2,
    rd_lo = (p1 - p2) - delta,
    rd_hi = (p1 - p2) + epsil
  )
}

# RISK DIFFERENCE BY SCORE VALUE ----------------------------------------------

message("\n== Risk difference by score value ==")

## pooled counts at each (score, cohort, score_value) -------------------------
# maxscores_ca_raw with analysis == "main" is already ED-filtered upstream
# (see materialize_variant_max() in site-side 03_analysis.R). score_name is
# already cleaned of the "_total" suffix by clean_score_names() in 00_load.R.

rd_counts = maxscores_ca_raw |>
  fsubset(tolower(analysis) == "main") |>
  fmutate(n_events = n * outcome) |>
  fgroup_by(score_name, ca_01, max_value) |>
  fsummarise(n = fsum(n), events = fsum(n_events)) |>
  fungroup() |>
  as.data.table()

## per-cell event rates with Wilson CIs ---------------------------------------

rd_counts[, c("p", "p_lo", "p_hi") := {
  ci = wilson_ci(events, n)
  list(fifelse(n > 0, events / n, NA_real_), ci$lo, ci$hi)
}]

## wide: one row per (score, max_value), cohort in columns --------------------

rd_wide = dcast(
  rd_counts,
  score_name + max_value ~ ca_01,
  value.var = c("n", "events", "p", "p_lo", "p_hi")
)

setnames(
  rd_wide,
  c("n_0", "n_1", "events_0", "events_1", "p_0", "p_1",
    "p_lo_0", "p_lo_1", "p_hi_0", "p_hi_1"),
  c("n_nc", "n_ca", "ev_nc", "ev_ca", "p_nc", "p_ca",
    "p_nc_lo", "p_ca_lo", "p_nc_hi", "p_ca_hi")
)

## Newcombe RD CIs (cancer minus noncancer) -----------------------------------
# dcast may leave NA in n_nc or n_ca for score values present in only one
# cohort. Newcombe helpers propagate NA cleanly.

rd_wide[, c("rd", "rd_lo", "rd_hi") := {
  ci = newcombe_rd_ci(ev_ca, n_ca, ev_nc, n_nc)
  list(ci$rd, ci$rd_lo, ci$rd_hi)
}]

## test for trend in RD across score values -----------------------------------
# Is the gap roughly constant (parallel shift), or does it fan/narrow?
# Fit weighted linear regression of RD on score value per score; weights are
# inverse approximate variance of the per-cell RD (Wald approx; appropriate
# for this trend test even though Newcombe is used for the CIs themselves).
# Slope near zero => parallel shift (threshold fix plausible).
#
# Guard against tail volatility: RD is bounded in [-1, 1] and the tails of
# most scores have tiny Ns where Wilson intervals widen, small-cell noise
# dominates, and the linear-trend assumption breaks down. Require a minimum
# cell size in BOTH cohorts (TREND_MIN_N) to contribute to the slope fit.

TREND_MIN_N = 100L

rd_trend = rd_wide[, {
  
  var_nc = p_nc * (1 - p_nc) / n_nc
  var_ca = p_ca * (1 - p_ca) / n_ca
  w      = 1 / (var_nc + var_ca)
  ok     = is.finite(w) & w > 0 & is.finite(rd) &
    !is.na(n_nc) & !is.na(n_ca) &
    n_nc >= TREND_MIN_N & n_ca >= TREND_MIN_N
  
  if (sum(ok) >= 3) {
    dt_fit = data.table(rd = rd[ok], max_value = max_value[ok], w = w[ok])
    fit    = lm(rd ~ max_value, data = dt_fit, weights = w)
    cf     = summary(fit)$coefficients
    list(
      slope_rd_per_pt = cf["max_value", "Estimate"],
      slope_se        = cf["max_value", "Std. Error"],
      slope_p         = cf["max_value", "Pr(>|t|)"],
      n_score_values  = sum(ok),
      trend_min_n     = TREND_MIN_N
    )
  } else {
    list(
      slope_rd_per_pt = NA_real_, slope_se = NA_real_,
      slope_p = NA_real_, n_score_values = sum(ok),
      trend_min_n = TREND_MIN_N
    )
  }
}, by = score_name]

## presentation-ready RD table (filter to adequately-populated cells) ---------

rd_wide[, has_both := !is.na(n_nc) & !is.na(n_ca) &
          n_nc >= MIN_N_PER_COHORT & n_ca >= MIN_N_PER_COHORT]

rd_table = rd_wide[has_both == TRUE, .(
  score_name,
  max_value,
  n_nc    = format_n(n_nc),
  rate_nc = sprintf("%.1f%%", 100 * p_nc),
  n_ca    = format_n(n_ca),
  rate_ca = sprintf("%.1f%%", 100 * p_ca),
  rd_pct  = sprintf("%.1f",   100 * rd),
  rd_ci   = sprintf("(%.1f, %.1f)", 100 * rd_lo, 100 * rd_hi)
)]

rd_table[, score_lab := factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)]
setorder(rd_table, score_lab, max_value)

message("  Risk differences: ", nrow(rd_wide), " cells across ",
        uniqueN(rd_wide$score_name), " scores; ", nrow(rd_table),
        " cells pass MIN_N_PER_COHORT = ", MIN_N_PER_COHORT)

# THRESHOLD EQUIVALENCE -------------------------------------------------------

message("\n== Threshold equivalence ==")

## operating characteristics at every integer threshold -----------------------
# Sweep t from 0 to max(max_value)+1; positivity if max_value >= t.

op_chars = copy(rd_counts)[, .(score_name, ca_01, max_value, n, events)]

op_totals = op_chars[, .(
  N_total  = sum(n),
  N_events = sum(events)
), by = .(score_name, ca_01)]

sweep_list = list()
for (sn in unique(op_chars$score_name)) {
  for (ca in c(0L, 1L)) {
    
    sub = op_chars[score_name == sn & ca_01 == ca]
    if (nrow(sub) == 0) next
    tot = op_totals[score_name == sn & ca_01 == ca]
    if (tot$N_total == 0) next
    
    thresholds = seq(0L, max(sub$max_value, na.rm = TRUE) + 1L)
    
    for (t in thresholds) {
      pos   = sub[max_value >= t]
      n_pos = sum(pos$n)
      tp    = sum(pos$events)
      fp    = n_pos - tp
      fn    = tot$N_events - tp
      tn    = tot$N_total - tot$N_events - fp
      
      sens = if (tot$N_events > 0) tp / tot$N_events else NA_real_
      spec = if ((tot$N_total - tot$N_events) > 0) tn / (tot$N_total - tot$N_events) else NA_real_
      ppv  = if (n_pos > 0) tp / n_pos else NA_real_
      alert_rate = n_pos / tot$N_total
      
      sweep_list[[length(sweep_list) + 1L]] = data.table(
        score_name = sn, ca_01 = ca, threshold = t,
        tp = tp, fp = fp, tn = tn, fn = fn,
        n_pos = n_pos, n_total = tot$N_total, n_events = tot$N_events,
        sens = sens, spec = spec, ppv = ppv, alert_rate = alert_rate
      )
    }
  }
}
sweep = rbindlist(sweep_list)

## Wilson CIs on sens/spec/ppv/alert_rate -------------------------------------

sweep[, c("sens_lo",  "sens_hi")  := wilson_ci(tp, tp + fn)]
sweep[, c("spec_lo",  "spec_hi")  := wilson_ci(tn, tn + fp)]
sweep[, c("ppv_lo",   "ppv_hi")   := wilson_ci(tp, tp + fp)]
sweep[, c("alert_lo", "alert_hi") := wilson_ci(n_pos, n_total)]

## flag the standard-threshold point per score -------------------------------

sweep[, is_std := threshold == STD_THRESHOLDS[score_name]]

message("  Threshold sweep: ", nrow(sweep), " (score x cohort x threshold) rows")

## match-by-sensitivity and match-by-alert-rate -------------------------------
# Find integer thresholds in the cancer cohort that bracket the noncancer
# standard-threshold operating point. Since scores are integer, exact matches
# are rare; report both the inclusive bracket (lower t, metric at or above
# target => more positives) and the conservative bracket (higher t, metric
# at or below target => fewer positives).

find_bracketing_thresholds = function(sweep_cohort, target,
                                      metric = c("sens", "alert_rate")) {
  
  metric = match.arg(metric)
  v = sweep_cohort[[metric]]
  t = sweep_cohort$threshold
  ord = order(t); v = v[ord]; t = t[ord]
  
  # sens and alert_rate are both nonincreasing in threshold
  at_or_above = which(v >= target)
  at_or_below = which(v <= target)
  
  if (length(at_or_above) == 0 || length(at_or_below) == 0) {
    return(list(t_inclusive = NA_integer_, t_conservative = NA_integer_))
  }
  
  list(
    t_inclusive    = t[max(at_or_above)],  # highest t where metric still >= target
    t_conservative = t[min(at_or_below)]   # lowest  t where metric has fallen <= target
  )
}

pick_closer_threshold = function(t_incl, t_cons, target, ca_sweep, metric) {
  
  if (is.na(t_incl) && is.na(t_cons)) return(NA_integer_)
  if (is.na(t_incl)) return(t_cons)
  if (is.na(t_cons)) return(t_incl)
  
  r_incl = ca_sweep[threshold == t_incl]
  r_cons = ca_sweep[threshold == t_cons]
  if (nrow(r_incl) == 0) return(t_cons)
  if (nrow(r_cons) == 0) return(t_incl)
  
  d_incl = abs(r_incl[[metric]] - target)
  d_cons = abs(r_cons[[metric]] - target)
  
  # tie -> more conservative (fewer alerts)
  if (d_incl < d_cons) t_incl else t_cons
}

match_rows = list()
for (sn in unique(sweep$score_name)) {
  
  nc = sweep[score_name == sn & ca_01 == 0L]
  ca = sweep[score_name == sn & ca_01 == 1L]
  if (nrow(nc) == 0 || nrow(ca) == 0) next
  
  t_std  = STD_THRESHOLDS[sn]
  std_nc = nc[threshold == t_std]
  std_ca = ca[threshold == t_std]
  if (nrow(std_nc) == 0) next
  
  br_sens  = find_bracketing_thresholds(ca, std_nc$sens,       "sens")
  br_alert = find_bracketing_thresholds(ca, std_nc$alert_rate, "alert_rate")
  
  t_pick_sens  = pick_closer_threshold(br_sens$t_inclusive,  br_sens$t_conservative,
                                       std_nc$sens,  ca, "sens")
  t_pick_alert = pick_closer_threshold(br_alert$t_inclusive, br_alert$t_conservative,
                                       std_nc$alert_rate, ca, "alert_rate")
  
  # operating characteristics at each bracket so the presentation table can
  # show both sides of the bracket (readers pick their own operating point)
  pull_row = function(sweep_cohort, t) {
    if (is.na(t)) return(sweep_cohort[FALSE])
    sweep_cohort[threshold == t]
  }
  
  at_sens_incl = pull_row(ca, br_sens$t_inclusive)
  at_sens_cons = pull_row(ca, br_sens$t_conservative)
  at_alert_incl = pull_row(ca, br_alert$t_inclusive)
  at_alert_cons = pull_row(ca, br_alert$t_conservative)
  
  at_sens  = pull_row(ca, t_pick_sens)
  at_alert = pull_row(ca, t_pick_alert)
  
  ca_at_std = function(col) if (nrow(std_ca)) std_ca[[col]] else NA_real_
  grab      = function(dt, col) if (nrow(dt))  dt[[col]] else NA_real_
  
  match_rows[[length(match_rows) + 1L]] = data.table(
    score_name        = sn,
    
    # noncancer standard operating point (the target)
    thr_std           = t_std,
    sens_target       = std_nc$sens,
    spec_target       = std_nc$spec,
    ppv_target        = std_nc$ppv,
    alert_target      = std_nc$alert_rate,
    
    # cancer at the same threshold (baseline)
    sens_ca_std       = ca_at_std("sens"),
    spec_ca_std       = ca_at_std("spec"),
    ppv_ca_std        = ca_at_std("ppv"),
    alert_ca_std      = ca_at_std("alert_rate"),
    
    # sens-matching: bracket thresholds and operating characteristics ---------
    t_sens_inclusive        = br_sens$t_inclusive,
    sens_at_sens_incl       = grab(at_sens_incl, "sens"),
    spec_at_sens_incl       = grab(at_sens_incl, "spec"),
    alert_at_sens_incl      = grab(at_sens_incl, "alert_rate"),
    
    t_sens_conservative     = br_sens$t_conservative,
    sens_at_sens_cons       = grab(at_sens_cons, "sens"),
    spec_at_sens_cons       = grab(at_sens_cons, "spec"),
    alert_at_sens_cons      = grab(at_sens_cons, "alert_rate"),
    
    # closer-of-two (preserved for backward compat with any figure code)
    thr_match_sens          = t_pick_sens,
    sens_match_sens         = grab(at_sens, "sens"),
    spec_match_sens         = grab(at_sens, "spec"),
    ppv_match_sens          = grab(at_sens, "ppv"),
    alert_match_sens        = grab(at_sens, "alert_rate"),
    
    # alert-rate-matching: bracket thresholds and operating characteristics ---
    t_alert_inclusive       = br_alert$t_inclusive,
    sens_at_alert_incl      = grab(at_alert_incl, "sens"),
    spec_at_alert_incl      = grab(at_alert_incl, "spec"),
    alert_at_alert_incl     = grab(at_alert_incl, "alert_rate"),
    
    t_alert_conservative    = br_alert$t_conservative,
    sens_at_alert_cons      = grab(at_alert_cons, "sens"),
    spec_at_alert_cons      = grab(at_alert_cons, "spec"),
    alert_at_alert_cons     = grab(at_alert_cons, "alert_rate"),
    
    # closer-of-two (preserved for backward compat with any figure code)
    thr_match_alert         = t_pick_alert,
    sens_match_alert        = grab(at_alert, "sens"),
    spec_match_alert        = grab(at_alert, "spec"),
    ppv_match_alert         = grab(at_alert, "ppv"),
    alert_match_alert       = grab(at_alert, "alert_rate")
  )
}
thr_eq = rbindlist(match_rows)

## presentation table ---------------------------------------------------------
# Per reviewer feedback: show BOTH bracketing thresholds (inclusive and
# conservative) rather than collapsing to the closer-of-two. Lets the reader
# pick their own operating point and see the bracket width. The "closer"
# values are still available in thr_eq for downstream figures.

pct = function(x) fifelse(is.na(x), "—", sprintf("%.1f%%", 100 * x))
int = function(x) fifelse(is.na(x), "—", as.character(x))

thr_eq_table = thr_eq[, .(
  score_lab = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS),
  thr_std,
  
  # noncancer target
  `NC sens`       = pct(sens_target),
  `NC spec`       = pct(spec_target),
  `NC alert rate` = pct(alert_target),
  
  # cancer at same (unshifted) threshold
  `CA sens @ std`  = pct(sens_ca_std),
  `CA spec @ std`  = pct(spec_ca_std),
  `CA alert @ std` = pct(alert_ca_std),
  
  # sens-matching: both brackets -------------------------------------------
  `Sens-match thr (incl)`   = int(t_sens_inclusive),
  `CA sens @ S-incl`        = pct(sens_at_sens_incl),
  `CA spec @ S-incl`        = pct(spec_at_sens_incl),
  `CA alert @ S-incl`       = pct(alert_at_sens_incl),
  
  `Sens-match thr (cons)`   = int(t_sens_conservative),
  `CA sens @ S-cons`        = pct(sens_at_sens_cons),
  `CA spec @ S-cons`        = pct(spec_at_sens_cons),
  `CA alert @ S-cons`       = pct(alert_at_sens_cons),
  
  # alert-matching: both brackets ------------------------------------------
  `Alert-match thr (incl)`  = int(t_alert_inclusive),
  `CA sens @ A-incl`        = pct(sens_at_alert_incl),
  `CA spec @ A-incl`        = pct(spec_at_alert_incl),
  `CA alert @ A-incl`       = pct(alert_at_alert_incl),
  
  `Alert-match thr (cons)`  = int(t_alert_conservative),
  `CA sens @ A-cons`        = pct(sens_at_alert_cons),
  `CA spec @ A-cons`        = pct(spec_at_alert_cons),
  `CA alert @ A-cons`       = pct(alert_at_alert_cons)
)]

setorder(thr_eq_table, score_lab)
message("  Threshold-equivalence table built for ", nrow(thr_eq_table), " scores")

# EFFICIENCY CURVES -----------------------------------------------------------

message("\n== Efficiency curves ==")

## For each (score, cohort), every integer threshold yields one (sens, alert)
## point. Default orientation: alert rate on x, sensitivity on y -- lower curve
## is better (more sens per alert), consistent with Romero-Brufau et al. The
## plotting layer decides orientation; we export both variables.
##
## Drop degenerate endpoints (t=0 with sens=1, or t>max with sens=0) for
## visual cleanliness, but retain the standard-threshold point regardless.

eff_data = copy(sweep)[, .(
  score_name, ca_01, threshold, is_std,
  sens, sens_lo, sens_hi,
  alert_rate, alert_lo, alert_hi,
  n_pos, n_events
)]

eff_data = eff_data[(sens > 0 & sens < 1) | is_std == TRUE]

eff_data[, `:=`(
  cohort_lab = fifelse(ca_01 == 1L, "Cancer", "Non-cancer"),
  score_lab  = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
)]

message("  Efficiency curve data: ", nrow(eff_data), " points across ",
        uniqueN(eff_data$score_name), " scores")

# SUPPLEMENT: PER-SITE HETEROGENEITY CHECK ------------------------------------

message("\n== Per-site sens/alert-rate check (supplement) ==")

## sesp_raw carries tp/fp/tn/fn at each site at the standard positivity
## definition. NOTE: this includes the NEWS any-3 rule for NEWS, unlike the
## aggregate-score-only sweep above. Flagged on the output.

if (exists("sesp_raw") && nrow(sesp_raw) > 0) {
  
  site_ops = copy(as.data.table(sesp_raw))
  
  # strip "_total" suffix to match cleaned score_name elsewhere in the pipeline.
  # clean_score_names() in 00_load.R does not touch sesp_raw, so we do it here.
  if ("score_name" %in% names(site_ops)) {
    site_ops[, score_name := str_remove(score_name, "_total$")]
  }
  
  site_ops = site_ops[, .(
    site, score_name, ca_01, tp, fp, tn, fn,
    n_total    = tp + fp + tn + fn,
    n_pos      = tp + fp,
    sens       = fifelse((tp + fn) > 0, tp / (tp + fn), NA_real_),
    spec       = fifelse((tn + fp) > 0, tn / (tn + fp), NA_real_),
    alert_rate = fifelse((tp + fp + tn + fn) > 0,
                         (tp + fp) / (tp + fp + tn + fn), NA_real_)
  )]
  
  site_ops[, c("sens_lo",  "sens_hi")  := wilson_ci(tp, tp + fn)]
  site_ops[, c("alert_lo", "alert_hi") := wilson_ci(n_pos, n_total)]
  
  site_ops[, `:=`(
    cohort_lab = fifelse(ca_01 == 1L, "Cancer", "Non-cancer"),
    score_lab  = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS),
    # NEWS only: flag that site-level sens/alert_rate reflects
    # aggregate-plus-any3, unlike the sweep-based main analysis
    news_any3_warning = score_name == "news"
  )]
  
  message("  Site-level ops characteristics: ", nrow(site_ops), " rows")
  
} else {
  message("  sesp_raw not available; skipping per-site check")
  site_ops = data.table()
}

# NEWS FIELD-CONVENTIONAL OPERATING POINT (for Figure 3) ----------------------
#
# The aggregate-score-only sweep above cannot represent the field-conventional
# NEWS rule ("aggregate >= 5 OR any single parameter >= 3") because the any-3
# escape clause is not an integer threshold. But Table 2 (and the manuscript's
# main Results text) report NEWS sensitivity/PPV/NPV using this rule, which
# is carried by sesp_raw -> site_ops above.
#
# Figure 3 needs a single pooled operating point per cohort (non-cancer,
# cancer) computed under the field-conventional rule, so it can plot that
# point off the aggregate-only curve. We pool tp/fp/tn/fn across sites for
# NEWS and recompute sens and alert_rate from the pooled counts, matching the
# pooled-counts convention used by the rest of this script.

if (nrow(site_ops) > 0 && any(site_ops$score_name == "news")) {
  
  news_fc_operating = site_ops[score_name == "news", .(
    tp = sum(tp, na.rm = TRUE),
    fp = sum(fp, na.rm = TRUE),
    tn = sum(tn, na.rm = TRUE),
    fn = sum(fn, na.rm = TRUE)
  ), by = ca_01]
  
  news_fc_operating[, `:=`(
    n_total    = tp + fp + tn + fn,
    sens       = fifelse((tp + fn) > 0, tp / (tp + fn), NA_real_),
    alert_rate = fifelse((tp + fp + tn + fn) > 0,
                         (tp + fp) / (tp + fp + tn + fn), NA_real_)
  )]
  
  # Retain only the columns Figure 3 needs; drop the counts now that
  # derived quantities are computed.
  news_fc_operating = news_fc_operating[, .(ca_01, sens, alert_rate)]
  setorder(news_fc_operating, ca_01)
  
  message("  NEWS field-conventional operating point:")
  message(sprintf("    non-cancer: sens=%.3f, alert_rate=%.3f",
                  news_fc_operating[ca_01 == 0, sens],
                  news_fc_operating[ca_01 == 0, alert_rate]))
  message(sprintf("    cancer:     sens=%.3f, alert_rate=%.3f",
                  news_fc_operating[ca_01 == 1, sens],
                  news_fc_operating[ca_01 == 1, alert_rate]))
  
} else {
  message("  NEWS rows not present in site_ops; news_fc_operating not built")
  news_fc_operating = data.table(
    ca_01 = integer(), sens = numeric(), alert_rate = numeric()
  )
}

# PER-SITE THRESHOLD SWEEP (Simpson's paradox check) --------------------------
#
# Pooled counts answer "what happens at our 8 sites if we move the threshold."
# A legitimate concern is that site-level heterogeneity in cancer mix and in
# baseline deterioration rate could make the pooled curve a poor summary of
# what any individual site would see. Rather than switching primary analysis
# to meta-analyzed logits (which adds variance for an operational question),
# we compute the threshold sweep PER SITE and report heterogeneity alongside
# the pooled result. If any score shows large site-level variation we will
# present the meta-analyzed version for that score as a contingency; by
# default the pooled sweep is primary.
#
# Built from maxscores_ca_raw (which carries `site` from the federated load)
# so it is internally consistent with the main sweep.

message("\n== Per-site threshold sweep ==")

## per-site counts at each (score, cohort, score_value) -----------------------

site_counts = maxscores_ca_raw |>
  fsubset(tolower(analysis) == "main") |>
  fmutate(n_events = n * outcome) |>
  fgroup_by(site, score_name, ca_01, max_value) |>
  fsummarise(n = fsum(n), events = fsum(n_events)) |>
  fungroup() |>
  as.data.table()

site_totals = site_counts[, .(
  N_total  = sum(n),
  N_events = sum(events)
), by = .(site, score_name, ca_01)]

## per-site threshold sweep ---------------------------------------------------
# Mirror the pooled sweep logic exactly, but one additional stratification.

site_sweep_list = list()
for (st in unique(site_counts$site)) {
  for (sn in unique(site_counts$score_name)) {
    for (ca in c(0L, 1L)) {
      
      sub = site_counts[site == st & score_name == sn & ca_01 == ca]
      if (nrow(sub) == 0) next
      tot = site_totals[site == st & score_name == sn & ca_01 == ca]
      if (nrow(tot) == 0 || tot$N_total == 0) next
      
      thresholds = seq(0L, max(sub$max_value, na.rm = TRUE) + 1L)
      
      for (t in thresholds) {
        pos   = sub[max_value >= t]
        n_pos = sum(pos$n)
        tp    = sum(pos$events)
        fp    = n_pos - tp
        fn    = tot$N_events - tp
        tn    = tot$N_total - tot$N_events - fp
        
        sens       = if (tot$N_events > 0) tp / tot$N_events else NA_real_
        spec       = if ((tot$N_total - tot$N_events) > 0) tn / (tot$N_total - tot$N_events) else NA_real_
        alert_rate = n_pos / tot$N_total
        
        site_sweep_list[[length(site_sweep_list) + 1L]] = data.table(
          site = st, score_name = sn, ca_01 = ca, threshold = t,
          tp = tp, fp = fp, tn = tn, fn = fn,
          n_pos = n_pos, n_total = tot$N_total, n_events = tot$N_events,
          sens = sens, spec = spec, alert_rate = alert_rate
        )
      }
    }
  }
}
site_sweep = rbindlist(site_sweep_list)

site_sweep[, c("sens_lo",  "sens_hi")  := wilson_ci(tp, tp + fn)]
site_sweep[, c("alert_lo", "alert_hi") := wilson_ci(n_pos, n_total)]
site_sweep[, is_std := threshold == STD_THRESHOLDS[score_name]]

site_sweep[, `:=`(
  cohort_lab = fifelse(ca_01 == 1L, "Cancer", "Non-cancer"),
  score_lab  = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
)]

message("  Per-site sweep: ", nrow(site_sweep), " (site x score x cohort x threshold) rows")

# HETEROGENEITY DIAGNOSTIC AT STANDARD THRESHOLD -------------------------------
#
# For each score x cohort, summarize how much the site-level operating
# characteristics at the standard threshold vary around the pooled estimate.
# Flag scores with wide site variation for potential meta-analyzed reporting.
#
# Small-site guard: a site with very few encounters in a cohort will produce
# extreme sens/alert_rate estimates from small-cell noise, not true
# heterogeneity. Exclude site x cohort cells below HET_MIN_N_PER_SITE from
# the heterogeneity summary; retain them in site_sweep for completeness.

HET_MIN_N_PER_SITE = 500L

het_sweep_std = site_sweep[is_std == TRUE & n_total >= HET_MIN_N_PER_SITE]

het_diag = het_sweep_std[, .(
  k_sites        = .N,
  sens_min       = min(sens, na.rm = TRUE),
  sens_median    = median(sens, na.rm = TRUE),
  sens_max       = max(sens, na.rm = TRUE),
  sens_range     = max(sens, na.rm = TRUE) - min(sens, na.rm = TRUE),
  alert_min      = min(alert_rate, na.rm = TRUE),
  alert_median   = median(alert_rate, na.rm = TRUE),
  alert_max      = max(alert_rate, na.rm = TRUE),
  alert_range    = max(alert_rate, na.rm = TRUE) - min(alert_rate, na.rm = TRUE)
), by = .(score_name, ca_01)]

# compare site-level point estimates to the pooled sweep at std threshold
pooled_std = sweep[is_std == TRUE, .(
  score_name, ca_01,
  sens_pooled  = sens,
  alert_pooled = alert_rate
)]

het_diag = merge(het_diag, pooled_std, by = c("score_name", "ca_01"), all.x = TRUE)

het_diag[, `:=`(
  sens_max_dev_from_pooled  = pmax(abs(sens_min  - sens_pooled),
                                   abs(sens_max  - sens_pooled)),
  alert_max_dev_from_pooled = pmax(abs(alert_min - alert_pooled),
                                   abs(alert_max - alert_pooled))
)]

## thresholds for flagging concerning heterogeneity ---------------------------
# Deliberately conservative: a score is "flagged" if the worst site differs
# from the pooled estimate by more than 10 percentage points on either metric.
# Flag does NOT mean results are invalid -- it means the meta-analyzed version
# should be reported for that score as a sensitivity check.

HET_FLAG_THRESHOLD = 0.10

het_diag[, flag_sens  := sens_max_dev_from_pooled  > HET_FLAG_THRESHOLD]
het_diag[, flag_alert := alert_max_dev_from_pooled > HET_FLAG_THRESHOLD]
het_diag[, flag_any   := flag_sens | flag_alert]

het_diag[, `:=`(
  cohort_lab = fifelse(ca_01 == 1L, "Cancer", "Non-cancer"),
  score_lab  = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
)]

setorder(het_diag, score_lab, ca_01)

n_flagged = het_diag[flag_any == TRUE, .N]
if (n_flagged > 0) {
  message("  HETEROGENEITY FLAG: ", n_flagged, " (score x cohort) cells exceed ",
          HET_FLAG_THRESHOLD, " pp max deviation from pooled. See het_diag_final.")
  print(het_diag[flag_any == TRUE,
                 .(score_lab, cohort_lab, k_sites,
                   sens_max_dev = round(sens_max_dev_from_pooled, 3),
                   alert_max_dev = round(alert_max_dev_from_pooled, 3))])
} else {
  message("  No (score x cohort) cells exceed ", HET_FLAG_THRESHOLD,
          " pp max deviation from pooled at the standard threshold.")
}

# PRESENTATION TABLES (flextable) ---------------------------------------------
#
# Three variants of the threshold-equivalence table:
#   (A) Wide: all brackets in one 24-column table, grouped headers, for
#       supplement.
#   (B) Sens-matching subtable: candidate for main text. 12 columns, shows
#       both inclusive and conservative brackets for sensitivity matching.
#   (C) Alert-matching subtable: same structure as (B) but for alert-rate
#       matching.

library(flextable)
library(officer)

if (!dir.exists(here("output", "tables"))) {
  dir.create(here("output", "tables"), recursive = TRUE)
}

## formatting helpers --------------------------------------------------------

pct_fmt = function(x) {
  fifelse(is.na(x), "\u2014", sprintf("%.1f%%", 100 * x))
}
int_fmt = function(x) {
  fifelse(is.na(x), "\u2014", as.character(x))
}

## (A) wide table with grouped headers ---------------------------------------

message("\n== Building threshold-equivalence wide table (grouped headers) ==")

thr_wide = thr_eq[, .(
  
  score_lab = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS),
  thr_std,
  
  # non-cancer at standard threshold
  sens_nc    = pct_fmt(sens_target),
  spec_nc    = pct_fmt(spec_target),
  alert_nc   = pct_fmt(alert_target),
  
  # cancer at standard threshold (unshifted)
  sens_ca_s  = pct_fmt(sens_ca_std),
  spec_ca_s  = pct_fmt(spec_ca_std),
  alert_ca_s = pct_fmt(alert_ca_std),
  
  # sens-match inclusive bracket
  t_si       = int_fmt(t_sens_inclusive),
  sens_si    = pct_fmt(sens_at_sens_incl),
  spec_si    = pct_fmt(spec_at_sens_incl),
  alert_si   = pct_fmt(alert_at_sens_incl),
  
  # sens-match conservative bracket
  t_sc       = int_fmt(t_sens_conservative),
  sens_sc    = pct_fmt(sens_at_sens_cons),
  spec_sc    = pct_fmt(spec_at_sens_cons),
  alert_sc   = pct_fmt(alert_at_sens_cons),
  
  # alert-match inclusive bracket
  t_ai       = int_fmt(t_alert_inclusive),
  sens_ai    = pct_fmt(sens_at_alert_incl),
  spec_ai    = pct_fmt(spec_at_alert_incl),
  alert_ai   = pct_fmt(alert_at_alert_incl),
  
  # alert-match conservative bracket
  t_ac       = int_fmt(t_alert_conservative),
  sens_ac    = pct_fmt(sens_at_alert_cons),
  spec_ac    = pct_fmt(spec_at_alert_cons),
  alert_ac   = pct_fmt(alert_at_alert_cons)
)]

setorder(thr_wide, score_lab)

thr_wide_headers = data.frame(
  col_keys = names(thr_wide),
  top = c(
    "", "",
    "Non-cancer (target)", "Non-cancer (target)", "Non-cancer (target)",
    "Cancer @ std threshold", "Cancer @ std threshold", "Cancer @ std threshold",
    "Cancer: sens-match (inclusive)", "Cancer: sens-match (inclusive)",
    "Cancer: sens-match (inclusive)", "Cancer: sens-match (inclusive)",
    "Cancer: sens-match (conservative)", "Cancer: sens-match (conservative)",
    "Cancer: sens-match (conservative)", "Cancer: sens-match (conservative)",
    "Cancer: alert-match (inclusive)", "Cancer: alert-match (inclusive)",
    "Cancer: alert-match (inclusive)", "Cancer: alert-match (inclusive)",
    "Cancer: alert-match (conservative)", "Cancer: alert-match (conservative)",
    "Cancer: alert-match (conservative)", "Cancer: alert-match (conservative)"
  ),
  bottom = c(
    "Score", "Std thr",
    "Sens", "Spec", "Alert",
    "Sens", "Spec", "Alert",
    "Thr", "Sens", "Spec", "Alert",
    "Thr", "Sens", "Spec", "Alert",
    "Thr", "Sens", "Spec", "Alert",
    "Thr", "Sens", "Spec", "Alert"
  ),
  stringsAsFactors = FALSE
)

ft_wide = flextable(thr_wide) |>
  set_header_df(mapping = thr_wide_headers, key = "col_keys") |>
  merge_h(part = "header") |>
  merge_v(part = "header") |>
  theme_booktabs() |>
  align(align = "center", part = "all") |>
  align(j = 1, align = "left", part = "body") |>
  fontsize(size = 8, part = "all") |>
  padding(padding = 2, part = "all") |>
  autofit()

save_as_docx(
  "Threshold equivalence (full, supplement)" = ft_wide,
  path = here("output", "tables", "table_threshold_eq_wide.docx")
)

## (B) sens-matching subtable (main text candidate) --------------------------

message("\n== Building sens-matching subtable ==")

thr_sens_tbl = thr_eq[, .(
  score_lab  = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS),
  thr_std,
  sens_nc    = pct_fmt(sens_target),
  alert_nc   = pct_fmt(alert_target),
  
  sens_ca_s  = pct_fmt(sens_ca_std),
  alert_ca_s = pct_fmt(alert_ca_std),
  
  t_si       = int_fmt(t_sens_inclusive),
  sens_si    = pct_fmt(sens_at_sens_incl),
  alert_si   = pct_fmt(alert_at_sens_incl),
  
  t_sc       = int_fmt(t_sens_conservative),
  sens_sc    = pct_fmt(sens_at_sens_cons),
  alert_sc   = pct_fmt(alert_at_sens_cons)
)]

setorder(thr_sens_tbl, score_lab)

thr_sens_headers = data.frame(
  col_keys = names(thr_sens_tbl),
  top = c(
    "", "",
    "Non-cancer (target)", "Non-cancer (target)",
    "Cancer @ std threshold", "Cancer @ std threshold",
    "Cancer: matched on sensitivity, inclusive bracket",
    "Cancer: matched on sensitivity, inclusive bracket",
    "Cancer: matched on sensitivity, inclusive bracket",
    "Cancer: matched on sensitivity, conservative bracket",
    "Cancer: matched on sensitivity, conservative bracket",
    "Cancer: matched on sensitivity, conservative bracket"
  ),
  bottom = c(
    "Score", "Std thr",
    "Sens", "Alert rate",
    "Sens", "Alert rate",
    "Thr", "Sens", "Alert rate",
    "Thr", "Sens", "Alert rate"
  ),
  stringsAsFactors = FALSE
)

ft_sens = flextable(thr_sens_tbl) |>
  set_header_df(mapping = thr_sens_headers, key = "col_keys") |>
  merge_h(part = "header") |>
  merge_v(part = "header") |>
  theme_booktabs() |>
  align(align = "center", part = "all") |>
  align(j = 1, align = "left", part = "body") |>
  fontsize(size = 9, part = "all") |>
  padding(padding = 3, part = "all") |>
  autofit()

save_as_docx(
  "Threshold equivalence: sensitivity-matched" = ft_sens,
  path = here("output", "tables", "table_threshold_eq_sens.docx")
)

## (C) alert-matching subtable (main text candidate) -------------------------

message("\n== Building alert-matching subtable ==")

thr_alert_tbl = thr_eq[, .(
  score_lab  = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS),
  thr_std,
  sens_nc    = pct_fmt(sens_target),
  alert_nc   = pct_fmt(alert_target),
  
  sens_ca_s  = pct_fmt(sens_ca_std),
  alert_ca_s = pct_fmt(alert_ca_std),
  
  t_ai       = int_fmt(t_alert_inclusive),
  sens_ai    = pct_fmt(sens_at_alert_incl),
  alert_ai   = pct_fmt(alert_at_alert_incl),
  
  t_ac       = int_fmt(t_alert_conservative),
  sens_ac    = pct_fmt(sens_at_alert_cons),
  alert_ac   = pct_fmt(alert_at_alert_cons)
)]

setorder(thr_alert_tbl, score_lab)

thr_alert_headers = data.frame(
  col_keys = names(thr_alert_tbl),
  top = c(
    "", "",
    "Non-cancer (target)", "Non-cancer (target)",
    "Cancer @ std threshold", "Cancer @ std threshold",
    "Cancer: matched on alert rate, inclusive bracket",
    "Cancer: matched on alert rate, inclusive bracket",
    "Cancer: matched on alert rate, inclusive bracket",
    "Cancer: matched on alert rate, conservative bracket",
    "Cancer: matched on alert rate, conservative bracket",
    "Cancer: matched on alert rate, conservative bracket"
  ),
  bottom = c(
    "Score", "Std thr",
    "Sens", "Alert rate",
    "Sens", "Alert rate",
    "Thr", "Sens", "Alert rate",
    "Thr", "Sens", "Alert rate"
  ),
  stringsAsFactors = FALSE
)

ft_alert = flextable(thr_alert_tbl) |>
  set_header_df(mapping = thr_alert_headers, key = "col_keys") |>
  merge_h(part = "header") |>
  merge_v(part = "header") |>
  theme_booktabs() |>
  align(align = "center", part = "all") |>
  align(j = 1, align = "left", part = "body") |>
  fontsize(size = 9, part = "all") |>
  padding(padding = 3, part = "all") |>
  autofit()

save_as_docx(
  "Threshold equivalence: alert-rate-matched" = ft_alert,
  path = here("output", "tables", "table_threshold_eq_alert.docx")
)

message("  Threshold-equivalence tables written to output/tables/")

# EXPORTS ---------------------------------------------------------------------

message("\n== Threshold-equivalence analyses complete ==")

rd_counts_final  = rd_counts
rd_wide_final    = rd_wide
rd_trend_final   = rd_trend
rd_table_final   = rd_table

sweep_final        = sweep
thr_eq_final       = thr_eq
thr_eq_table_final = thr_eq_table

eff_data_final   = eff_data
site_ops_final   = site_ops      # per-site ops at standard threshold (sesp_raw-based)
news_fc_operating_final = news_fc_operating  # pooled NEWS field-conventional operating point (Figure 3)

# Simpson's-paradox / site heterogeneity exports -----------------------------
site_sweep_final = site_sweep    # per-site sweep across all integer thresholds
het_diag_final   = het_diag      # heterogeneity diagnostic at standard threshold

# Flextable objects (also saved to output/tables/ above) ---------------------
ft_thr_wide_final  = ft_wide
ft_thr_sens_final  = ft_sens
ft_thr_alert_final = ft_alert

gc()
