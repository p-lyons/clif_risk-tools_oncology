# 02b_carryforward.R — carry-forward sensitivity analysis (P3)
#
# Reviewers asked whether the results depend on how long a stale vital sign is
# carried forward. This script re-runs the LOCF-plus-totals stage at vitals
# windows of 2, 6, and 12 hours (the laboratory window is held fixed at 12h,
# since the concern is vitals/oxygenation monitoring intensity) and exports, per
# window, the encounter-level maximum scores and the 24-hour horizon counts.
#
# The extraction and point-assignment stage runs ONCE (in 02_scores.R, which
# persists the pre-carry components); here only the carry-forward and downstream
# stages repeat. Correctness test: the cf6 outputs must equal the main round-two
# outputs exactly. Rather than depend on 03a having run (it runs after this
# script), "main" is computed here directly from scores_full.parquet with the
# same functions used for the cf windows, and cf6 is asserted equal to it.
#
#   Inputs (proj_tables/):
#     scores_components.parquet  point-assigned components, pre-LOCF
#     news_o2_stream.parquet     pre-carry supplemental-oxygen stream
#     ward_times.parquet         ward-stay intervals
#     cohort.parquet             patient/hospital ids, code status, ca/ed flags
#     outcome_times.parquet      composite outcome time for truncation + horizon
#     scores_full.parquet        the round-two main scored data (cf6 == main gate)

if (!exists("BOX_DIR")) {
  stop("BOX_DIR not found. Did you run 00_setup first?", call. = FALSE)
}

if (!exists("site_lowercase")) {
  site_lowercase = as.character(read_parquet(here("proj_tables", "site_lowercase.parquet"))$site_lowercase)
}

# parameters -------------------------------------------------------------------

VITALS_CF_HOURS = c(2L, 6L, 12L)
LABS_CF_HOURS   = 12L
REFERENCE_H     = 6L

# inputs -----------------------------------------------------------------------

components   = as.data.table(read_parquet(here("proj_tables", "scores_components.parquet")))
news_o2_base = as.data.table(read_parquet(here("proj_tables", "news_o2_stream.parquet")))
ward_times   = read_parquet(here("proj_tables", "ward_times.parquet"))
cohort       = read_parquet(here("proj_tables", "cohort.parquet"))
scores_full  = as.data.table(read_parquet(here("proj_tables", "scores_full.parquet")))

outcomes = as.data.table(read_parquet(here("proj_tables", "outcome_times.parquet")))
## a structurally all-missing dttm column reads back as logical; coerce so pmin
## stays in the datetime domain (outcome_dttm is populated in practice)
if (is.logical(outcomes$outcome_dttm)) {
  set(outcomes, j = "outcome_dttm", value = .POSIXct(rep(NA_real_, nrow(outcomes)), tz = "UTC"))
}
outcomes_min = outcomes[, .(joined_hosp_id, outcome_dttm)]

# cohort-derived helpers -------------------------------------------------------

jp          = select(cohort, patient_id, joined_hosp_id, hospital_id)
fc          = fsubset(cohort, tolower(initial_code_status) == "full")$joined_hosp_id
cancer_encs = fsubset(cohort, ca_01 == 1)$joined_hosp_id
ed_encs     = fsubset(cohort, ed_admit_01 == 1)$joined_hosp_id

# helpers ----------------------------------------------------------------------

make_y = function(h_to_event, horizon) {
  as.integer(!is.na(h_to_event) & h_to_event >= 0 & h_to_event <= horizon)
}

## LOCF within N hours for one column (identical to 02_scores.R)
locf <- function(df, col, hours, id = "joined_hosp_id", tcol = "time") {

  cy     <- as.difftime(hours, units = "hours")
  tstamp <- paste0("t__", col)
  is_int <- is.integer(df[[col]])
  na_val <- if (is_int) NA_integer_ else NA_real_

  df |>
    arrange(!!sym(id), !!sym(tcol)) |>
    mutate(!!tstamp := if_else(!is.na(.data[[col]]), .data[[tcol]], as.POSIXct(NA_real_))) |>
    fill(!!sym(tstamp), .direction = "down", .by = !!sym(id)) |>
    mutate(
      !!col := {
        v  <- data.table::nafill(.data[[col]], type = "locf")
        tt <- .data[[tcol]] - .data[[tstamp]]
        if_else(tt <= cy, v, na_val)
      },
      .by = !!sym(id)
    ) |>
    select(-!!sym(tstamp))
}

# build the scored table at a given vitals carry-forward window ----------------
# Reproduces the post-extraction pipeline of 02_scores.R with the vitals LOCF
# window as a parameter. At vitals_h = 6 this reproduces scores_full.parquet on
# every column that feeds the maxscores/counts artifacts (verified by the
# cf6 == main assertion below).

build_scores_core = function(vitals_h) {

  dt = copy(components)

  score_cols = setdiff(names(dt), c("joined_hosp_id", "time"))
  lab_cols   = intersect(c("sirs_wbc", "sirs_co2", "sirs_bands"), score_cols)
  vital_cols = setdiff(score_cols, lab_cols)

  for (cn in vital_cols) dt = locf(dt, cn, hours = vitals_h)
  for (cn in lab_cols)   dt = locf(dt, cn, hours = LABS_CF_HOURS)

  ## news_o2: 6h-style as-of carry at the same vitals window (news_o2 only)
  setDT(dt)
  o2 = news_o2_base[, .(joined_hosp_id, time, news_o2)]
  setkey(o2, joined_hosp_id, time)
  dt = o2[dt, on = .(joined_hosp_id, time), roll = vitals_h * 3600]
  setDT(dt)

  dt = replace_na(dt, value = 0L, cols = NULL, set = F, type = "const")
  setDT(dt)

  ## sirs preprep (identical to 02_scores.R)
  dt$sirs_rr = if_else(dt$sirs_rr == 1 | dt$sirs_co2 == 1, 1L, 0L)
  dt = fselect(dt, -sirs_co2, -spo2)
  dt = rename(dt, sf = mews_sf)

  ## totals
  setDT(dt)
  dt[, sirs_total  := rowSums(.SD, na.rm = T), .SDcols = patterns("^sirs_")]
  dt[, mews_total  := rowSums(.SD, na.rm = T), .SDcols = patterns("^mews_")]
  dt[, news_total  := rowSums(.SD, na.rm = T), .SDcols = patterns("^news_")]
  dt[, qsofa_total := rowSums(.SD, na.rm = T), .SDcols = patterns("^qsofa_")]
  dt[, news_any3 := as.integer(
    news_temp == 3L | news_hr == 3L | news_rr == 3L |
      news_sbp == 3L | news_spo2 == 3L | news_gcs == 3L
  )]

  ## wards only (identical to 02_scores.R)
  dt =
    tidytable(dt) |>
    ftransform(mews_sf_total = mews_total + sf) |>
    join(ward_times, how = "inner", multiple = T) |>
    fsubset(time >= in_dttm & time <= out_dttm) |>
    select(-hospitalization_id) |>
    fgroup_by(joined_hosp_id, in_dttm, out_dttm, time) |>
    fmax() |>
    fgroup_by(joined_hosp_id) |>
    fmutate(in_dttm = fmin(in_dttm), out_dttm = fmax(out_dttm)) |>
    fungroup() |>
    funique()

  ## composite outcome + truncation (row set matches 02_scores.R exactly)
  dt =
    join(dt, outcomes_min, how = "left", multiple = T) |>
    fmutate(o_primary_01 = if_else(!is.na(outcome_dttm), 1L, 0L)) |>
    fmutate(end_dttm     = pmin(out_dttm, outcome_dttm, na.rm = T)) |>
    fsubset(time < end_dttm)

  ## cancer / ED flags (membership, matching 02_scores.R)
  dt = as.data.table(dt)
  dt[, ca_01       := fifelse(joined_hosp_id %in% cancer_encs, 1L, 0L)]
  dt[, ed_admit_01 := fifelse(joined_hosp_id %in% ed_encs,     1L, 0L)]

  dt[]
}

# attach patient/hospital ids, as 03_analysis.R does before the analyses -------

prep = function(scored) {
  join(scored, jp, how = "inner", multiple = FALSE)
}

# encounter-level maximum scores (materialize_variant_max "main" + aggregate) --

make_maxscores = function(scored) {

  dt_scores = fsubset(scored, ed_admit_01 == 1L)

  dt_max =
    fgroup_by(dt_scores, joined_hosp_id) |>
    fsummarize(
      patient_id    = ffirst(patient_id),
      hospital_id   = ffirst(hospital_id),
      ca_01         = ffirst(ca_01),
      ed_admit_01   = ffirst(ed_admit_01),
      fullcode_01   = fifelse(ffirst(joined_hosp_id) %in% fc, 1L, 0L),
      sirs_max      = fmax(sirs_total,    na.rm = TRUE),
      qsofa_max     = fmax(qsofa_total,   na.rm = TRUE),
      mews_max      = fmax(mews_total,    na.rm = TRUE),
      news_max      = fmax(news_total,    na.rm = TRUE),
      mews_sf_max   = fmax(mews_sf_total, na.rm = TRUE),
      enc_news_any3 = fmax(news_any3,     na.rm = TRUE),
      outcome       = ffirst(o_primary_01)
    ) |>
    pivot_longer(cols = ends_with("_max"), names_to = "score_name", values_to = "max_value") |>
    ftransform(score_name = str_remove(score_name, "_max"))

  dt_max = as.data.table(dt_max)
  dt_max[is.infinite(max_value), max_value := NA_integer_]
  dt_max[is.infinite(enc_news_any3), enc_news_any3 := 0L]

  dt_max[, .(n = .N), by = .(score_name, ca_01, max_value, outcome)][, site := site_lowercase][]
}

# 24-hour horizon counts (scores_long_base "main" + run_horizon_counts h24) -----

make_counts_h24 = function(scored) {

  dt = as.data.table(scored)[ed_admit_01 == 1L]

  long = melt(
    dt,
    id.vars       = c("ca_01", "time", "outcome_dttm"),
    measure.vars  = c("sirs_total", "qsofa_total", "mews_total", "news_total", "mews_sf_total"),
    variable.name = "score_name",
    value.name    = "value"
  )

  long[, score_name := as.character(score_name)]
  long[, h_to_event := as.numeric(difftime(outcome_dttm, time, units = "hours"))]
  long[, outcome := make_y(h_to_event, 24L)]

  long[, .(n = .N), by = .(score_name, ca_01, value, outcome)][, `:=`(site = site_lowercase, h = 24L)][]
}

# main (round-two) reference from scores_full.parquet --------------------------

message("\n== Carry-forward sensitivity (P3) ==")
message("  Computing main-variant reference from scores_full.parquet")

main_scored = prep(scores_full)
main_max    = make_maxscores(main_scored)
main_counts = make_counts_h24(main_scored)

# output destination -----------------------------------------------------------

sens_dir = here(BOX_DIR, "sensitivity")
if (!dir.exists(sens_dir)) dir.create(sens_dir, recursive = TRUE, showWarnings = FALSE)

# run each window --------------------------------------------------------------

cf6_max      = NULL
cf6_counts   = NULL
mean_max_all = list()

for (H in VITALS_CF_HOURS) {

  message(sprintf("  Window cf%d: building scored table and artifacts", H))

  scored_H = prep(build_scores_core(H))
  mx = make_maxscores(scored_H)
  ct = make_counts_h24(scored_H)

  fwrite(mx, file.path(sens_dir, paste0("maxscores-ca-cf",     H, "-", site_lowercase, ".csv")))
  fwrite(ct, file.path(sens_dir, paste0("counts-ca-h24-cf",    H, "-", site_lowercase, ".csv")))

  mean_max_all[[as.character(H)]] =
    mx[!is.na(max_value), .(mean_max = sum(max_value * n) / sum(n)), by = score_name][, H := H]

  if (H == REFERENCE_H) {
    cf6_max    = mx
    cf6_counts = ct
  }

  rm(scored_H, mx, ct); gc()
}

# QC: cf6 must equal the main round-two outputs exactly ------------------------

if (!fsetequal(cf6_max, main_max)) {
  stop("cf6 maxscores do not equal the main-variant maxscores; the carry-forward refactor has drifted.",
       call. = FALSE)
}
if (!fsetequal(cf6_counts, main_counts)) {
  stop("cf6 24h counts do not equal the main-variant counts; the carry-forward refactor has drifted.",
       call. = FALSE)
}
message("✅ cf6 == main (maxscores and 24h counts).")

# QC: mean maximum score is monotone non-decreasing in the window length -------

mean_max_tbl = rbindlist(mean_max_all, use.names = TRUE)
setorder(mean_max_tbl, score_name, H)

mono = mean_max_tbl[, .(ok = all(diff(mean_max) >= -1e-9)), by = score_name]
if (!all(mono$ok)) {
  bad = mono[ok == FALSE]$score_name
  stop(sprintf("Mean max score not monotone non-decreasing in window length for: %s",
               paste(bad, collapse = ", ")),
       call. = FALSE)
}
message("✅ Mean max score is monotone non-decreasing across cf2 < cf6 < cf12 for every score.")

message(sprintf("✅ Carry-forward artifacts written to: %s", sens_dir))

# cleanup ----------------------------------------------------------------------

rm(components, news_o2_base, ward_times, cohort, scores_full, outcomes, outcomes_min,
   jp, fc, cancer_encs, ed_encs, main_scored, main_max, main_counts,
   cf6_max, cf6_counts, mean_max_all, mean_max_tbl, mono,
   locf, make_y, build_scores_core, prep, make_maxscores, make_counts_h24)
gc()

# go to 02c

################################################################################
