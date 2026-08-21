# 02c_monitoring.R — monitoring intensity and missingness by cancer status (P4)
#
# Reviewers asked whether monitoring intensity differs between patients with and
# without cancer, since a difference could by itself produce a score-performance
# gap. This script quantifies, per encounter and measure, how often each vital
# and the leukocyte count are actually recorded during the at-risk ward period,
# and re-runs the missingness export stratified by cancer status.
#
#   Inputs (from proj_tables/):
#     vital_lab_extract.parquet  one row per (encounter, measure, measurement
#                                time) for the six vitals and WBC, written by
#                                02_scores.R (spans the whole encounter)
#     scores_full.parquet        supplies, per encounter, the ward-entry time
#                                (in_dttm) and the composite at-risk end
#                                (end_composite_dttm), plus ca_01 / ed_admit_01
#     cohort.parquet             ED-admit scope and ca_01 for the missingness
#                                export (kept identical to the round-one scope)
#
# Note on inputs: the spec lists the extract + cohort.parquet, but ward-hours are
# defined as difftime(end_composite_dttm, in_dttm) and neither field is in
# cohort.parquet, so scores_full.parquet is read for the per-encounter ward
# window. Nothing individual-level is exported; every artifact is aggregated to
# counts and sums by cancer status.

if (!exists("BOX_DIR")) {
  stop("BOX_DIR not found. Did you run 00_setup first?", call. = FALSE)
}

if (!exists("site_lowercase")) {
  site_lowercase = as.character(read_parquet(here("proj_tables", "site_lowercase.parquet"))$site_lowercase)
}

# inputs -----------------------------------------------------------------------

extract     = as.data.table(read_parquet(here("proj_tables", "vital_lab_extract.parquet")))
scores_full = as.data.table(read_parquet(here("proj_tables", "scores_full.parquet")))
cohort      = as.data.table(read_parquet(here("proj_tables", "cohort.parquet")))

VITAL_MEASURES = c("heart_rate", "respiratory_rate", "temp_c", "sbp", "spo2", "gcs", "wbc")

message("\n== Monitoring intensity and missingness (P4) ==")

# per-encounter main frame and ward-hours at risk ------------------------------
# The main-variant encounter set is the ED-admitted encounters that carry score
# rows. in_dttm and end_composite_dttm are per-encounter constants in
# scores_full, so one row per encounter after unique().

enc = unique(scores_full[
  , .(joined_hosp_id, ca_01, ed_admit_01, in_dttm, end_composite_dttm)
])
enc = enc[ed_admit_01 == 1L]

enc[, ward_hours := as.numeric(difftime(end_composite_dttm, in_dttm, units = "hours"))]

## QC: ward-hours must be strictly positive
if (any(enc$ward_hours <= 0)) {
  stop(sprintf("Ward-hours <= 0 for %d encounter(s).", sum(enc$ward_hours <= 0)),
       call. = FALSE)
}

N_main = nrow(enc)

# measurement counts within the at-risk ward window ----------------------------
# The extract spans the whole encounter, so it is intersected with each
# encounter's [in_dttm, end_composite_dttm) window — the same interval the score
# grid occupies — before counting. Extract rows are unique per (encounter,
# measure, time), so a row count is a distinct-measurement-time count.

ext = merge(
  extract,
  enc[, .(joined_hosp_id, in_dttm, end_composite_dttm)],
  by = "joined_hosp_id"
)
ext = ext[time >= in_dttm & time < end_composite_dttm]

obs = ext[, .(n_obs = .N), by = .(joined_hosp_id, measure)]

# encounter × measure grid (every main encounter appears once per measure, with
# n_obs = 0 where a measure was never recorded) --------------------------------

grid = enc[, .(measure = VITAL_MEASURES), by = .(joined_hosp_id, ca_01, ward_hours)]
grid = merge(grid, obs, by = c("joined_hosp_id", "measure"), all.x = TRUE)
grid[is.na(n_obs), n_obs := 0L]
grid[, rate_per24h := n_obs / ward_hours * 24]

# monitoring-ca ----------------------------------------------------------------

monitoring = grid[
  , .(
    n_enc             = .N,
    sum_obs           = sum(n_obs),
    sumsq_obs         = sum(as.numeric(n_obs)^2),
    sum_ward_hours    = sum(ward_hours),
    sum_rate_per24h   = sum(rate_per24h),
    sumsq_rate_per24h = sum(rate_per24h^2)
  ),
  by = .(ca_01, measure)
][, site := site_lowercase][]

# monitoring_bins-ca -----------------------------------------------------------

bin_breaks = c(0, 1, 2, 4, 6, 12, Inf)
bin_labs   = c("[0,1)", "[1,2)", "[2,4)", "[4,6)", "[6,12)", "12+")

grid[, rate_bin := as.character(
  cut(rate_per24h, breaks = bin_breaks, labels = bin_labs, right = FALSE)
)]

monitoring_bins = grid[
  , .(n_enc = .N), by = .(ca_01, measure, rate_bin)
][, site := site_lowercase][]

# QC: encounters accounted for exactly once per measure in each artifact -------

chk_mon = monitoring[, .(tot = sum(n_enc)), by = measure]
if (any(chk_mon$tot != N_main)) {
  stop("monitoring n_enc does not sum to the main-variant encounter count for every measure.",
       call. = FALSE)
}

chk_bin = monitoring_bins[, .(tot = sum(n_enc)), by = measure]
if (any(chk_bin$tot != N_main)) {
  stop("monitoring_bins n_enc does not sum to the main-variant encounter count for every measure.",
       call. = FALSE)
}

message(sprintf("✅ Monitoring QC passed | %s main encounters | %d measures.",
                format(N_main, big.mark = ","), length(VITAL_MEASURES)))

# write monitoring artifacts ---------------------------------------------------

diag_dir = here(BOX_DIR, "diagnostics")
if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(monitoring,      file.path(diag_dir, paste0("monitoring-ca-",      site_lowercase, ".csv")))
fwrite(monitoring_bins, file.path(diag_dir, paste0("monitoring_bins-ca-", site_lowercase, ".csv")))

# missing_vlab-ca --------------------------------------------------------------
# Re-runs the round-one missing_vlab export stratified by cancer status. Scope
# and variable set match round one: ED-admitted cohort encounters (a superset of
# the main analysis set by the 34 encounters that carry no score rows), and the
# six variables the original tracked -- systolic blood pressure was not among
# them. "Missing" means the encounter has no measurement of that type anywhere.

miss_map = c(
  heart_rate = "heart_rate",
  resp_rate  = "respiratory_rate",
  temp       = "temp_c",
  spo2       = "spo2",
  gcs        = "gcs",
  wbc        = "wbc"
)

scope = cohort[ed_admit_01 == 1L, .(joined_hosp_id, ca_01)]

missing_vlab = rbindlist(lapply(seq_along(miss_map), function(i) {

  variable_label = names(miss_map)[i]
  extract_measure = miss_map[[i]]
  has_measure    = unique(extract[measure == extract_measure]$joined_hosp_id)

  scope[
    , .(
      variable  = variable_label,
      n_total   = .N,
      n_missing = sum(!joined_hosp_id %in% has_measure)
    ),
    by = ca_01
  ]
}))

missing_vlab[, pct_missing := round(100 * n_missing / n_total, 2)]
missing_vlab[, site := site_lowercase]
setcolorder(missing_vlab, c("variable", "ca_01", "n_total", "n_missing", "pct_missing", "site"))

fwrite(missing_vlab, file.path(diag_dir, paste0("missing_vlab-ca-", site_lowercase, ".csv")))

message(sprintf("✅ Monitoring and missingness written to: %s", diag_dir))

# cleanup ----------------------------------------------------------------------

rm(extract, scores_full, cohort, enc, ext, obs, grid,
   monitoring, monitoring_bins, chk_mon, chk_bin,
   scope, missing_vlab, miss_map)
gc()

# go to 03a

################################################################################
