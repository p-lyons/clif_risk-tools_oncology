# 03b_leadtime.R — lead time from first threshold crossing to deterioration (P6)
#
# Reviewer 4 asked whether an alert affords an actionable window. For each
# encounter, score, and outcome this reports the hours from the first
# threshold-positive score to the event, exported as binned counts plus
# sufficient statistics (never individual intervals), together with a four-cell
# crossed × event classification of every encounter.
#
#   Inputs (proj_tables/):  scores_full.parquet, cohort.parquet, outcome_times.parquet
#   Reads (upload_to_box_v2/threshold/): the P5 sesp artifact, for a QC cross-check
#   Outcomes: composite and nohospice only.
#
# Positivity uses the standard THRESHOLDS with the NEWS any-3 rule, exactly as in
# 03a's ever_positive block. scores_full is truncated before the event, so every
# crossing precedes it and lead times are strictly positive (asserted, not
# filtered). The crossclass denominator is the full cohort ED-admit set (from
# outcome_times), matching 03a's sesp, so the QC ties out.

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

LEAD_BREAKS = c(0, 2, 4, 8, 12, 24, 48, 72, 168, Inf)
LEAD_LABELS = c("(0,2]", "(2,4]", "(4,8]", "(8,12]", "(12,24]", "(24,48]", "(48,72]", "(72,168]", ">168")

# helpers ----------------------------------------------------------------------

#' First threshold-positive time per (encounter, score), ED-admit only.
#' Standard threshold crossings for all five scores, plus the NEWS single-
#' parameter any-3 rule folded into news_total; earliest wins.
compute_crossings = function(scores, thresholds) {

  sc_ed = as.data.table(scores)[ed_admit_01 == 1L]

  long = melt(sc_ed,
    id.vars       = c("joined_hosp_id", "time"),
    measure.vars  = thresholds$score_name,
    variable.name = "score_name", value.name = "value")
  long[, score_name := as.character(score_name)]
  long = merge(long, thresholds, by = "score_name")

  cross_std = long[value >= threshold][
    order(time), .(cross_time = time[1L]), by = .(joined_hosp_id, score_name)]

  cross_n3 = sc_ed[news_any3 == 1L][
    order(time), .(cross_time = time[1L]), by = joined_hosp_id][
    , score_name := "news_total"]

  rbindlist(list(cross_std, cross_n3), use.names = TRUE)[
    order(cross_time), .(cross_time = cross_time[1L]), by = .(joined_hosp_id, score_name)]
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
    cohort[ed_admit_01 == 1L, .(joined_hosp_id, ca_01)],
    outcome_tt[, .(joined_hosp_id, o_out = get(flag))],
    by = "joined_hosp_id", all.x = TRUE
  )
  ae[is.na(o_out), o_out := 0L]

  # per-encounter event time (for lead), from scores
  enc_ev = scores[ed_admit_01 == 1L, .(event_dttm = get(evc)[1L]), by = joined_hosp_id]

  # encounter × score grid, classified by crossed × event
  grid = ae[, .(score_name = THRESHOLDS$score_name), by = .(joined_hosp_id, ca_01, o_out)]
  grid = merge(grid, crossings, by = c("joined_hosp_id", "score_name"), all.x = TRUE)
  grid[, `:=`(crossed = as.integer(!is.na(cross_time)), event = o_out)]

  crossclass = grid[
    , .(n = .N), by = .(score_name, ca_01, crossed, event)
  ][, `:=`(outcome_key = outcome_key, site = site_lowercase)]
  setcolorder(crossclass, c("score_name", "ca_01", "outcome_key", "crossed", "event", "n", "site"))

  # lead time for encounters that crossed AND had the event
  lead = merge(grid[crossed == 1L & event == 1L], enc_ev, by = "joined_hosp_id", all.x = TRUE)
  lead[, lead_h := as.numeric(difftime(event_dttm, cross_time, units = "hours"))]

  if (nrow(lead) > 0L && any(lead$lead_h < 0, na.rm = TRUE)) {
    stop(sprintf("Negative lead time for outcome '%s' (%d rows).",
                 outcome_key, sum(lead$lead_h < 0, na.rm = TRUE)), call. = FALSE)
  }

  lead[, lead_bin := as.character(cut(lead_h, breaks = LEAD_BREAKS, labels = LEAD_LABELS, right = TRUE))]
  leadtime = lead[
    , .(n = .N, sum_hours = sum(lead_h), sumsq_hours = sum(lead_h^2)),
    by = .(score_name, ca_01, lead_bin)
  ][, `:=`(outcome_key = outcome_key, site = site_lowercase)]
  setcolorder(leadtime, c("score_name", "ca_01", "outcome_key", "lead_bin", "n", "sum_hours", "sumsq_hours", "site"))

  # QC: crossclass events (summed over crossed) reproduce the sesp event counts
  sesp_path = here(BOX_DIR, "threshold", paste0("sesp-ca-", outcome_key, "-", site_lowercase, ".csv"))
  if (file.exists(sesp_path)) {
    sesp   = fread(sesp_path)
    cc_ev  = crossclass[event == 1L, .(n_ev = sum(n)), by = .(score_name, ca_01)]
    chk    = merge(sesp[, .(score_name, ca_01, n_outcome)], cc_ev, by = c("score_name", "ca_01"))
    if (nrow(chk) != nrow(cc_ev) || any(chk$n_outcome != chk$n_ev)) {
      stop(sprintf("crossclass event counts disagree with sesp for outcome '%s'.", outcome_key),
           call. = FALSE)
    }
  } else {
    message("  (sesp artifact not found for ", outcome_key, "; skipping cross-check)")
  }

  list(leadtime = leadtime, crossclass = crossclass)
}

# run --------------------------------------------------------------------------

message("\n== 03b: lead time to deterioration ==")

scores     = read_parquet(here("proj_tables", "scores_full.parquet"))
cohort     = read_parquet(here("proj_tables", "cohort.parquet"))
outcome_tt = as.data.table(read_parquet(here("proj_tables", "outcome_times.parquet")))

crossings = compute_crossings(scores, THRESHOLDS)

thr_dir = here(BOX_DIR, "threshold")
if (!dir.exists(thr_dir)) dir.create(thr_dir, recursive = TRUE, showWarnings = FALSE)

for (ok in LEADTIME_OUTCOMES) {
  message("  outcome: ", ok)
  res = leadtime_for_outcome(ok, scores, cohort, outcome_tt, crossings, site_lowercase)
  fwrite(res$leadtime,   file.path(thr_dir, paste0("leadtime-ca-",   ok, "-", site_lowercase, ".csv")))
  fwrite(res$crossclass, file.path(thr_dir, paste0("crossclass-ca-", ok, "-", site_lowercase, ".csv")))
}

message("\n== 03b complete ==")
message("Files written to: ", thr_dir)

# go to 03c

################################################################################
