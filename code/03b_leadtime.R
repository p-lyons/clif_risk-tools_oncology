# 03b_leadtime.R — lead time from threshold crossing to deterioration (P6)
#
# Reviewer 4 asked whether an alert affords an actionable window. For each
# encounter, score, and outcome this reports the hours from a threshold-positive
# score to the event, exported as poolable exceedance counts and as exact
# per-unit medians (never individual intervals), together with a four-cell
# crossed × event classification of every encounter.
#
#   Inputs (proj_tables/):  scores_full.parquet, cohort.parquet, outcome_times.parquet
#   Reads (upload_to_box_v2/threshold/): the P5 sesp artifact, for a QC cross-check
#   Outcomes: composite and nohospice only.
#
# Positivity uses the standard THRESHOLDS with the NEWS any-3 rule, exactly as in
# 03a's ever_positive block. scores_full is truncated before the event, so every
# crossing precedes it and lead times are non-negative (asserted, not filtered).
# The crossclass denominator is the full cohort ED-admit set (from outcome_times),
# matching 03a's sesp, so the QC ties out.
#
# ------------------------------------------------------------------------------
# Round-two redesign
#
# Round one exported nine coarse bins plus sum/sumsq, and the coordinating center
# reported a mean. Both choices failed on the observed distribution: 91% of SIRS
# crossings in the cancer cohort exceed 24 hours, so the median falls inside a
# 96-hour-wide bin (or the open top bin) where interpolation is undefined, and a
# mean of 182 hours with an SD of 255 describes a heavy right tail badly.
#
# Two artifacts replace the single binned file.
#
#   leadtime         Cumulative counts at fixed thresholds (6, 12, 24, 48, 72,
#                    96, 168, 336 hours) with the denominator. These pool by
#                    summation with no distributional assumption, and yield both
#                    "within x" and "beyond x" by subtraction at the
#                    coordinating center.
#
#   leadtime_median  Exact median and quartiles computed on line-level data at
#                    the site, reported per health system AND per hospital_id,
#                    each with its own n_events so the coordinating center can
#                    apply an event floor before reporting a range. Sites are
#                    NOT pooled into a single median; the range across units is
#                    the reported statistic, and the exceedance curve above
#                    supplies the encounter-weighted counterpart independently.
#
# Both carry crossing_def, taking two values.
#
#   first        Time from the FIRST threshold-positive score of the encounter
#                to the event. The conventional definition, and what a reader
#                expects. At the alert rates these scores carry (NEWS is
#                ever-positive in 71% of encounters) it fires early in the stay
#                and behaves partly as a length-of-stay measure.
#
#   final_onset  Time from the onset of the LAST uninterrupted positive run to
#                the event: the alert that came on and stayed on going into the
#                deterioration. Defined whenever the encounter ever crossed, so
#                the denominator matches `first` exactly and the two are
#                directly comparable. The gap between them measures how long a
#                score sat positive without an event, which is the alert-burden
#                argument in R13 and R14.

if (!exists("BOX_DIR")) {
  stop("BOX_DIR not found. Did you run 00_setup first?", call. = FALSE)
}
if (!exists("site_lowercase")) {
  site_lowercase = as.character(read_parquet(here("proj_tables", "site_lowercase.parquet"))$site_lowercase)
}

# constants --------------------------------------------------------------------

THRESHOLDS = data.table(
  score_name = c("sirs_total", "qsofa_total", "mews_total", "news_total", "mews_sf_total"),
  threshold  = c(2L, 2L, 5L, 5L, 7L)
)

LEADTIME_OUTCOMES = c("composite", "nohospice")
OUT_FLAG  = c(composite = "o_composite_01",      nohospice = "o_nohospice_01")
OUT_EVENT = c(composite = "event_composite_dttm", nohospice = "event_nohospice_dttm")

# Cumulative thresholds, in hours. The short end (6, 12, 24) carries
# final_onset, which is expected to concentrate under a day; the long end
# (72 through 336) carries first, whose median sits somewhere past two days.
# 48 is retained because two days is the landmark a clinician reads without
# conversion and because it halves the bracket the pooled curve places around
# the first-crossing median.
LEAD_THRESHOLDS = c(6, 12, 24, 48, 72, 96, 168, 336)

CROSSING_DEFS = c("first", "final_onset")

# helpers ----------------------------------------------------------------------

#' Threshold-positive crossing times per (encounter, score), ED-admit only,
#' under both crossing definitions.
#'
#' Positivity is the union of the standard threshold rule and, for NEWS, the
#' single-parameter any-3 rule. Round one computed those two as separate tables
#' and took the earlier; folding the any-3 rule into a single `pos` flag gives
#' the identical first-crossing time and additionally makes the positivity
#' SERIES available, which final_onset requires.
#'
#' Runs are identified with the standard cumsum-of-negatives trick: within an
#' encounter and score ordered by time, every positive row in the same
#' uninterrupted run shares a run_grp value, and later runs carry larger values.
#' The last run is therefore the one whose run_grp is maximal among positive
#' rows, and its onset is that run's earliest time.
compute_crossings = function(scores, thresholds) {

  sc_ed = as.data.table(scores)[ed_admit_01 == 1L]

  long = melt(sc_ed,
    id.vars       = c("joined_hosp_id", "time", "news_any3"),
    measure.vars  = thresholds$score_name,
    variable.name = "score_name", value.name = "value")
  long[, score_name := as.character(score_name)]
  long = merge(long, thresholds, by = "score_name")

  long[, pos := as.integer(!is.na(value) & value >= threshold)]
  long[score_name == "news_total" & !is.na(news_any3) & news_any3 == 1L, pos := 1L]

  setorder(long, joined_hosp_id, score_name, time)
  long[, run_grp := cumsum(pos == 0L), by = .(joined_hosp_id, score_name)]

  pos_rows = long[pos == 1L]
  pos_rows[, run_max := max(run_grp), by = .(joined_hosp_id, score_name)]

  first_cross = pos_rows[
    , .(cross_time = time[1L]), by = .(joined_hosp_id, score_name)
  ][, crossing_def := "first"]

  final_cross = pos_rows[
    run_grp == run_max, .(cross_time = time[1L]), by = .(joined_hosp_id, score_name)
  ][, crossing_def := "final_onset"]

  # The two definitions must describe the identical set of crossed encounters,
  # since final_onset is the onset of a run that first_cross also sees.
  if (nrow(first_cross) != nrow(final_cross)) {
    stop("Crossing definitions disagree on which encounters ever crossed.", call. = FALSE)
  }

  out = rbindlist(list(first_cross, final_cross), use.names = TRUE)

  # final_onset can never precede first: it is the onset of a later or identical run
  chk = dcast(out, joined_hosp_id + score_name ~ crossing_def, value.var = "cross_time")
  if (any(chk$final_onset < chk$first)) {
    stop("final_onset precedes first crossing; run identification is wrong.", call. = FALSE)
  }

  rm(long, pos_rows, chk)
  out[]
}

#' Lead-time + crossclass artifacts for one outcome.
leadtime_for_outcome = function(outcome_key, scores, cohort, outcome_tt, crossings,
                                site_lowercase) {

  flag = OUT_FLAG[[outcome_key]]
  evc  = OUT_EVENT[[outcome_key]]

  scores     = as.data.table(scores)
  cohort     = as.data.table(cohort)
  outcome_tt = as.data.table(outcome_tt)

  # full cohort ED-admit denominator with this outcome's flag (matches 03a sesp)
  ae = merge(
    cohort[ed_admit_01 == 1L, .(joined_hosp_id, ca_01, hospital_id)],
    outcome_tt[, .(joined_hosp_id, o_out = get(flag))],
    by = "joined_hosp_id", all.x = TRUE
  )
  ae[is.na(o_out), o_out := 0L]

  # per-encounter event time (for lead), from scores
  enc_ev = scores[ed_admit_01 == 1L, .(event_dttm = get(evc)[1L]), by = joined_hosp_id]

  # encounter × score grid. crossclass is a property of the encounter, not of the
  # crossing definition, so it is built from the `first` rows only; the crossed
  # set is identical either way (asserted in compute_crossings).
  cross_first = crossings[crossing_def == "first", .(joined_hosp_id, score_name, cross_time)]

  grid = ae[, .(score_name = THRESHOLDS$score_name),
            by = .(joined_hosp_id, ca_01, hospital_id, o_out)]
  grid = merge(grid, cross_first, by = c("joined_hosp_id", "score_name"), all.x = TRUE)
  grid[, `:=`(crossed = as.integer(!is.na(cross_time)), event = o_out)]

  crossclass = grid[
    , .(n = .N), by = .(score_name, ca_01, crossed, event)
  ][, `:=`(outcome_key = outcome_key, site = site_lowercase)]
  setcolorder(crossclass, c("score_name", "ca_01", "outcome_key", "crossed", "event", "n", "site"))

  # lead time for encounters that crossed AND had the event, under both definitions
  eligible = grid[crossed == 1L & event == 1L,
                  .(joined_hosp_id, score_name, ca_01, hospital_id)]

  lead = merge(eligible, crossings, by = c("joined_hosp_id", "score_name"),
               allow.cartesian = TRUE)
  lead = merge(lead, enc_ev, by = "joined_hosp_id", all.x = TRUE)
  lead[, lead_h := as.numeric(difftime(event_dttm, cross_time, units = "hours"))]

  if (nrow(lead) > 0L && any(lead$lead_h < 0, na.rm = TRUE)) {
    stop(sprintf("Negative lead time for outcome '%s' (%d rows).",
                 outcome_key, sum(lead$lead_h < 0, na.rm = TRUE)), call. = FALSE)
  }
  if (anyNA(lead$lead_h)) {
    stop(sprintf("Missing lead time for outcome '%s' (%d rows).",
                 outcome_key, sum(is.na(lead$lead_h))), call. = FALSE)
  }

  ## artifact 1: cumulative counts at fixed thresholds -------------------------
  # One row per (score, cohort, crossing_def, threshold). n_at_or_below and
  # n_total both pool by summation, so the coordinating center can report
  # "within x hours" directly and "beyond x hours" as n_total - n_at_or_below.

  lead_cum = lead[
    , .(
      threshold_h   = LEAD_THRESHOLDS,
      n_at_or_below = vapply(LEAD_THRESHOLDS, function(t) sum(lead_h <= t), numeric(1)),
      n_total       = .N
    ),
    by = .(score_name, ca_01, crossing_def)
  ][, `:=`(outcome_key = outcome_key, site = site_lowercase)]

  setcolorder(lead_cum, c("score_name", "ca_01", "outcome_key", "crossing_def",
                          "threshold_h", "n_at_or_below", "n_total", "site"))

  ## artifact 2: exact medians per unit ----------------------------------------
  # Computed on line-level lead_h at the site, so no interpolation is involved.
  # Emitted for the health system as a whole and for each hospital_id, each with
  # its own n_events. Sites are deliberately NOT pooled: the coordinating center
  # reports the range across units, and applies an event floor using n_events
  # rather than having that choice made here.

  # Expects dt to carry a unit_id column already; returns one row per
  # (score, cohort, crossing_def, unit_id).
  median_block = function(dt, unit_lab) {
    dt[
      , {
        q = as.numeric(quantile(lead_h, probs = c(0.25, 0.50, 0.75),
                                type = 7, names = FALSE))
        .(
          unit        = unit_lab,
          n_events    = .N,
          median_h    = q[2],
          q25_h       = q[1],
          q75_h       = q[3],
          mean_h      = mean(lead_h),
          sd_h        = if (.N > 1L) sd(lead_h) else NA_real_,
          sum_hours   = sum(lead_h),
          sumsq_hours = sum(lead_h^2)
        )
      },
      by = .(score_name, ca_01, crossing_def, unit_id)
    ]
  }

  lead_sys = copy(lead)[, unit_id := site_lowercase]
  lead_hos = copy(lead)[, unit_id := as.character(hospital_id)]

  lead_median = rbindlist(
    list(
      median_block(lead_sys, "health_system"),
      median_block(lead_hos, "hospital")
    ),
    use.names = TRUE
  )[, `:=`(outcome_key = outcome_key, site = site_lowercase)]

  setcolorder(lead_median, c("score_name", "ca_01", "outcome_key", "crossing_def",
                             "unit", "unit_id", "n_events", "median_h",
                             "q25_h", "q75_h", "mean_h", "sd_h",
                             "sum_hours", "sumsq_hours", "site"))

  ## QC ------------------------------------------------------------------------
  # The health-system rows must account for exactly the same events as the
  # threshold artifact's denominator, and the hospital rows must sum to them.

  denom_cum = unique(lead_cum[, .(score_name, ca_01, crossing_def, n_total)])
  denom_sys = lead_median[unit == "health_system",
                          .(score_name, ca_01, crossing_def, n_sys = n_events)]
  denom_hos = lead_median[unit == "hospital",
                          .(n_hos = sum(n_events)), by = .(score_name, ca_01, crossing_def)]

  chk = merge(denom_cum, denom_sys, by = c("score_name", "ca_01", "crossing_def"))
  chk = merge(chk,      denom_hos, by = c("score_name", "ca_01", "crossing_def"))

  if (any(chk$n_total != chk$n_sys) || any(chk$n_total != chk$n_hos)) {
    stop(sprintf("Lead-time denominators disagree across artifacts for outcome '%s'.",
                 outcome_key), call. = FALSE)
  }

  # The cumulative counts must be monotone non-decreasing in threshold_h.
  mono = lead_cum[order(threshold_h),
                  .(ok = all(diff(n_at_or_below) >= 0)),
                  by = .(score_name, ca_01, crossing_def)]
  if (any(!mono$ok)) {
    stop("Cumulative lead-time counts are not monotone in threshold.", call. = FALSE)
  }

  # final_onset can never exceed first, so its cumulative count at any threshold
  # must be at least as large.
  ord = dcast(lead_cum, score_name + ca_01 + threshold_h ~ crossing_def,
              value.var = "n_at_or_below")
  if (any(ord$final_onset < ord$first)) {
    stop("final_onset lead times exceed first-crossing lead times.", call. = FALSE)
  }

  # crossclass events (summed over crossed) reproduce the sesp event counts
  sesp_path = here(BOX_DIR, "threshold", paste0("sesp-ca-", outcome_key, "-", site_lowercase, ".csv"))
  if (file.exists(sesp_path)) {
    sesp   = fread(sesp_path)
    cc_ev  = crossclass[event == 1L, .(n_ev = sum(n)), by = .(score_name, ca_01)]
    chk2   = merge(sesp[, .(score_name, ca_01, n_outcome)], cc_ev, by = c("score_name", "ca_01"))
    if (nrow(chk2) != nrow(cc_ev) || any(chk2$n_outcome != chk2$n_ev)) {
      stop(sprintf("crossclass event counts disagree with sesp for outcome '%s'.", outcome_key),
           call. = FALSE)
    }
  } else {
    message("  (sesp artifact not found for ", outcome_key, "; skipping cross-check)")
  }

  list(leadtime = lead_cum, leadtime_median = lead_median, crossclass = crossclass)
}

# run --------------------------------------------------------------------------

message("\n== 03b: lead time to deterioration ==")

scores     = read_parquet(here("proj_tables", "scores_full.parquet"))
cohort     = read_parquet(here("proj_tables", "cohort.parquet"))
outcome_tt = as.data.table(read_parquet(here("proj_tables", "outcome_times.parquet")))

crossings = compute_crossings(scores, THRESHOLDS)

message("  Crossings: ",
        format(uniqueN(crossings[, .(joined_hosp_id, score_name)]), big.mark = ","),
        " encounter-score pairs, both definitions")

thr_dir = here(BOX_DIR, "threshold")
if (!dir.exists(thr_dir)) dir.create(thr_dir, recursive = TRUE, showWarnings = FALSE)

for (ok in LEADTIME_OUTCOMES) {

  message("  outcome: ", ok)

  res = leadtime_for_outcome(ok, scores, cohort, outcome_tt, crossings, site_lowercase)

  fwrite(res$leadtime,
         file.path(thr_dir, paste0("leadtime-ca-", ok, "-", site_lowercase, ".csv")))
  fwrite(res$leadtime_median,
         file.path(thr_dir, paste0("leadtime_median-ca-", ok, "-", site_lowercase, ".csv")))
  fwrite(res$crossclass,
         file.path(thr_dir, paste0("crossclass-ca-", ok, "-", site_lowercase, ".csv")))

  message(sprintf("    %d threshold rows | %d median rows (%d hospital unit(s))",
                  nrow(res$leadtime), nrow(res$leadtime_median),
                  uniqueN(res$leadtime_median[unit == "hospital"]$unit_id)))
}

message("\n== 03b complete ==")
message("Files written to: ", thr_dir)

# end of site pipeline; run_all.R builds the manifest next

################################################################################
