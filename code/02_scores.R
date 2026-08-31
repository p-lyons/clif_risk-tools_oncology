# create risk scores by looping over vital sign and lab components
# process one vital/lab at a time, assign all relevant scores, then combine

if (!exists("BOX_DIR")) {
  stop("BOX_DIR not found. Did you run 00_setup first?", call. = FALSE)
}

## allow_sparse_o2 is read and validated once in 00_setup.R, from
## config/config_clif_oncrisk.yaml, and carried through the stage-01 cleanup. It
## governs the sparse-oxygen guard at the end of this file. Nothing in this file
## is site-modifiable; a site that needs the override edits its config, not the
## pipeline.

if (!exists("allow_sparse_o2")) {
  stop("allow_sparse_o2 not found. Did you run 00_setup first?", call. = FALSE)
}

# extract vital signs and assign score points ----------------------------------

## make a list of data frames for each vital sign ------------------------------

get_each_vital <- function(vn) {
  data_list[["vitals"]] |>
    dplyr::filter(hospitalization_id %in% cohort_hids) |>
    dplyr::filter(vital_category == vn) |>
    dplyr::select(hospitalization_id, time = recorded_dttm, vital_value) |>
    dplyr::collect() |>
    join(hid_jid_crosswalk, how = "inner", multiple = T) |>
    fgroup_by(joined_hosp_id, time) |>
    fsummarize(val = fmax(vital_value)) |>
    roworder(joined_hosp_id, time)
}

vitals_list = setNames(map(req_vitals, get_each_vital), req_vitals)

## temp_c (sirs, mews, news) ---------------------------------------------------

vitals_list$temp_c =
  vitals_list[["temp_c"]] |>
  ftransform(
    sirs_temp = if_else(val > 38 | val < 36, 1L, 0L),
    mews_temp = case_when(
      val >= 38.5 ~ 2L,
      val >  35.0 ~ 0L,
      val <= 35.0 ~ 2L,
      TRUE        ~ NA_integer_
    ),
    news_temp = case_when(
      val >= 39.1 ~ 2L,
      val >= 38.1 ~ 1L,
      val >= 36.1 ~ 0L,
      val >= 35.1 ~ 1L,
      val <= 35.0 ~ 3L,
      TRUE        ~ NA_integer_
    )
  ) |>
  fselect(-val)

## heart rate (sirs, mews, news) -----------------------------------------------

vitals_list$heart_rate =
  vitals_list[["heart_rate"]] |>
  ftransform(
    sirs_hr = if_else(val > 90, 1L, 0L),
    mews_hr = case_when(
      val >= 130  ~ 3L,
      val >= 111  ~ 2L,
      val >= 101  ~ 1L,
      val >= 51   ~ 0L,
      val >= 41   ~ 1L,
      val <= 40   ~ 2L,
      TRUE        ~ NA_integer_
    ),
    news_hr = case_when(
      val >= 131  ~ 3L,
      val >= 111  ~ 2L,
      val >= 91   ~ 1L,
      val >= 51   ~ 0L,
      val >= 41   ~ 1L,
      val <= 40   ~ 3L,
      TRUE        ~ NA_integer_
    )
  ) |>
  fselect(-val)

## respiratory rate (sirs, qsofa, mews, news) ----------------------------------

vitals_list$respiratory_rate =
  vitals_list[["respiratory_rate"]] |>
  ftransform(
    sirs_rr  = if_else(val >= 21, 1L, 0L),
    qsofa_rr = if_else(val >= 22, 1L, 0L),
    mews_rr  = case_when(
      val >= 30 ~ 3L,
      val >= 21 ~ 2L,
      val >= 15 ~ 1L,
      val >= 9  ~ 0L,
      val <= 8  ~ 2L,
      TRUE      ~ NA_integer_
    ),
    news_rr = case_when(
      val >= 25 ~ 3L,
      val >= 21 ~ 2L,
      val >= 12 ~ 0L,
      val >= 9  ~ 1L,
      val <= 8  ~ 3L,
      TRUE       ~ NA_integer_
    )
  ) |>
  fselect(-val)

## sbp (qsofa, mews, news) -----------------------------------------------------

vitals_list$sbp =
  vitals_list[["sbp"]] |>
  ftransform(
    qsofa_sbp = if_else(val <= 100, 1L, 0L),
    mews_sbp  = case_when(
      val >= 200 ~ 2L,
      val >= 101 ~ 0L,
      val >= 81  ~ 1L,
      val >= 71  ~ 2L,
      val <= 70  ~ 3L,
      TRUE       ~ NA_integer_
    ),
    news_sbp = case_when(
      val >= 220 ~ 3L,
      val >= 111 ~ 0L,
      val >= 101 ~ 1L,
      val >= 91  ~ 2L,
      val <= 90  ~ 3L,
      TRUE       ~ NA_integer_
    )
  ) |>
  fselect(-val)

## gcs (qsofa) -----------------------------------------------------------------

vitals_list$gcs =
  data_list[["patient_assessments"]] |>
  dplyr::filter(hospitalization_id %in% cohort_hids) |>
  dplyr::filter(assessment_category == "gcs_total") |>
  dplyr::filter(!is.na(numerical_value) & numerical_value <= 15) |>
  dplyr::select(
    hospitalization_id,
    time = recorded_dttm,
    numerical_value
  ) |>
  dplyr::collect() |>
  ftransform(
    qsofa_gcs = if_else(numerical_value < 15, 1L, 0L, NA_integer_),
    mews_gcs  = case_when(
      numerical_value == 15 ~ 0L,
      numerical_value >= 13 ~ 1L,
      numerical_value >= 08 ~ 2L,
      numerical_value >= 03 ~ 3L,
      TRUE                  ~ NA_integer_
    ),
    news_gcs  = if_else(numerical_value < 15, 3L, 0L, NA_integer_)
  ) |>
  fselect(-numerical_value)

vitals_list$gcs =
  join(vitals_list$gcs, hid_jid_crosswalk, how = "inner", multiple = T) |>
  fgroup_by(joined_hosp_id, time) |>
  fsummarize(
    qsofa_gcs = fmax(qsofa_gcs),
    mews_gcs  = fmax(mews_gcs),
    news_gcs  = fmax(news_gcs)
  ) |>
  roworder(joined_hosp_id, time)

## spo2 (news) -----------------------------------------------------------------

vitals_list$spo2 =
  vitals_list[["spo2"]] |>
  ftransform(
    news_spo2 = case_when(
      val >= 96 ~ 0L,
      val >= 94 ~ 1L,
      val >= 92 ~ 2L,
      val <= 91 ~ 3L,
      TRUE      ~ NA_integer_
    )
  ) |>
  rename(spo2 = val)

## de-duplicate ----------------------------------------------------------------

vitals_list$heart_rate =
  vitals_list[["heart_rate"]] |>
  fgroup_by(joined_hosp_id, time) |>
  fmax()

vitals_list$respiratory_rate =
  vitals_list[["respiratory_rate"]] |>
  fgroup_by(joined_hosp_id, time) |>
  fmax()

vitals_list$sbp =
  vitals_list[["sbp"]] |>
  fgroup_by(joined_hosp_id, time) |>
  fmax()

vitals_list$temp_c =
  vitals_list[["temp_c"]] |>
  fgroup_by(joined_hosp_id, time) |>
  fmax()

vitals_list$gcs =
  vitals_list[["gcs"]] |>
  fgroup_by(joined_hosp_id, time) |>
  fmax()

vitals_list$spo2 =
  vitals_list[["spo2"]] |>
  fgroup_by(joined_hosp_id, time) |>
  fmax()

## sf ratio (mews_sf, pmid 32114753) -------------------------------------------

### prepare fio2 from respiratory table ----------------------------------------

resp =
  data_list[["respiratory_support"]] |>
  dplyr::filter(hospitalization_id %in% cohort_hids) |>
  dplyr::select(
    hospitalization_id,
    time = recorded_dttm,
    device_category,
    lpm_set,
    fio2_set
  ) |>
  dplyr::filter(!is.na(lpm_set) | !is.na(fio2_set) | tolower(device_category) == "room air") |>
  dplyr::collect()

resp$fio2_impute = case_when(
  resp$lpm_set == 0 ~ 0.21,
  resp$lpm_set <= 6 ~ 0.21 + (0.04*resp$lpm_set),
  resp$lpm_set >  0 ~ 0.21 + (6*0.04) + ((resp$lpm_set - 6)*0.03),
  TRUE              ~ NA_real_
)

resp$fio2 = case_when(
  resp$fio2_set >= 0.21 & resp$fio2_set <= 1   ~ resp$fio2_set,
  resp$fio2_set >= 21   & resp$fio2_set <= 100 ~ resp$fio2_set*0.01,
  resp$fio2_impute >= 1                        ~ 1,
  tolower(resp$device_category) == "room air"  ~ 0.21,
  TRUE                                         ~ resp$fio2_impute
)

### news supplemental-oxygen item: full-coverage stream ------------------------
# Published NEWS (Smith 2013) includes a 2-point supplemental-oxygen item that
# the round-one build omitted; round two completes the score. The item needs
# coverage of every encounter, so it is derived here from the
# full respiratory_support extract -- before the keep_ids filter below that
# restricts the SpO2/FiO2 (MEWS-SF) stream to encounters with an SpO2 value --
# reusing the fio2 computed just above verbatim.
#
#   news_o2 = 2  when fio2 > 0.21 or lpm_set > 0        (any supplemental O2)
#   news_o2 = 0  when room air, or (fio2 == 0.21 & lpm_set == 0)
#   news_o2 = NA otherwise  (device category alone cannot resolve the item)
#
# o2_src records which branch resolved each measurement -- 1 fio2 > 0.21, 2
# lpm_set > 0, 3 room air -- and is kept only for the resolution diagnostic.
# Grid-fixed handling (confirmed with Brenna): this stream is carried onto the
# existing score rows by a 6h as-of join later, rather than added as its own
# vitals_list element. Keeping it out of the Reduce/full-join grid means no new
# score-row timestamps are introduced, so SIRS, qSOFA, MEWS, and MEWS-SF rows
# stay byte-identical to round one and only NEWS changes.

news_o2_stream =
  ftransform(
    resp,
    news_o2 = case_when(
      fio2 > 0.21 | lpm_set > 0                     ~ 2L,
      tolower(device_category) == "room air"        ~ 0L,
      fio2 == 0.21 & lpm_set == 0                    ~ 0L,
      TRUE                                          ~ NA_integer_
    ),
    o2_src = case_when(
      fio2 > 0.21                                   ~ 1L,
      lpm_set > 0                                    ~ 2L,
      tolower(device_category) == "room air"        ~ 3L,
      fio2 == 0.21 & lpm_set == 0                    ~ 3L,
      TRUE                                          ~ NA_integer_
    )
  ) |>
  fsubset(!is.na(news_o2)) |>
  join(hid_jid_crosswalk, how = "inner", multiple = T) |>
  fgroup_by(joined_hosp_id, time) |>
  fsummarize(
    news_o2 = fmax(news_o2),
    o2_src  = fmin(o2_src)
  ) |>
  roworder(joined_hosp_id, time)

resp =
  join(resp, hid_jid_crosswalk, how = "inner", multiple = T) |>
  fgroup_by(joined_hosp_id, time) |>
  fsummarize(fio2 = fmax(fio2)) |>
  roworder(joined_hosp_id, time)

keep_ids = pull(vitals_list$spo2, joined_hosp_id)
resp     = fsubset(resp, !is.na(fio2) & joined_hosp_id %in% keep_ids)

### sf based on 6h fio2 carryforward -------------------------------------------

vitals_list$sf =
  fsubset(vitals_list$spo2, spo2 >= 0 & spo2 < 97) |>
  join(resp, how = "full", multiple = F) |>
  ftransform(tf = if_else(!is.na(fio2), time, as.POSIXct(NA))) |>
  roworder(joined_hosp_id, time) |>
  fill(tf,   .direction = "down", .by = joined_hosp_id) |>
  fill(fio2, .direction = "down", .by = joined_hosp_id) |>
  ftransform(hdf = as.numeric(difftime(time, tf), units = "hours")) |>
  ftransform(f2  = if_else(hdf <= 6 & !is.na(hdf), fio2, NA_real_)) |>
  ftransform(sf  = spo2/f2)

### compute sf scores ----------------------------------------------------------

vitals_list$sf =
  fsubset(vitals_list$sf, !is.na(spo2)) |>
  ftransform(
    mews_sf = case_when(
      f2 == 0.21 ~ 0L,
      sf <= 235  ~ 3L,
      sf <= 315  ~ 2L,
      TRUE       ~ 0L
    )
  ) |>
  fselect(joined_hosp_id, time, mews_sf)

# extract labs and assign score points -----------------------------------------

req_labs = c("pco2_arterial", "wbc")

labs =
  data_list[["labs"]] |>
  dplyr::filter(hospitalization_id %in% cohort_hids) |>
  dplyr::filter(lab_category %in% req_labs) |>
  dplyr::select(
    hospitalization_id,
    time = lab_result_dttm,
    lab_category,
    lab_value_numeric
  ) |>
  dplyr::collect()

labs =
  join(labs, hid_jid_crosswalk, how = "inner", multiple = T) |>
  fgroup_by(joined_hosp_id, time, lab_category) |>
  fsummarize(val = fmin(lab_value_numeric))

## wbc, pco2 (sirs) ------------------------------------------------------------

labs =
  pivot_wider(labs, names_from = lab_category, values_from = val) |>
  ftransform(sirs_wbc = if_else(wbc < 4 | wbc > 12, 1L, 0L)) |>
  ftransform(sirs_co2 = if_else(pco2_arterial < 32, 1L, 0L)) |>
  select(joined_hosp_id, time, starts_with("sirs")) |>
  fsubset(!is.na(sirs_wbc) | !is.na(sirs_co2))

# NOTE: encounter-level vitals/lab missingness now lives in 02c_monitoring.R,
# which reports it stratified by cancer status (missing_vlab-ca). The round-one
# non-stratified missing_vlab export was removed here to complete that swap.

# combine score components -----------------------------------------------------

scores =
  Reduce(
    function(x, y) join(x, y,
                        on            = c("joined_hosp_id", "time"),
                        how           = "full",
                        multiple      = TRUE,
                        drop.dup.cols = TRUE
    ),
    c(vitals_list, list(labs = labs))
  )

if ("val" %in% names(scores)) scores <- select(scores, -val)

setDT(scores)
score_component_cols = setdiff(names(scores), c("joined_hosp_id", "time"))
scores = scores[, lapply(.SD, max, na.rm = TRUE), by = .(joined_hosp_id, time), .SDcols = score_component_cols]

for (col in score_component_cols) {
  scores[is.infinite(get(col)), (col) := NA_integer_]
}

# persist intermediates for downstream sensitivity/monitoring scripts ----------
# These are written here, before the extracts are dropped, so P3 and P4 do not
# re-extract raw data.
#
#   scores_components.parquet  Point-assigned score components, one row per
#                              encounter-time, BEFORE carry-forward, ward
#                              restriction, and outcome truncation. P3
#                              (02b_carryforward.R) re-runs the LOCF-plus-totals
#                              stage on this at 2/6/12h; the cf6 pass must
#                              reproduce this script's main output exactly.
#   news_o2_stream.parquet     The supplemental-oxygen measurement stream BEFORE
#                              the 6h carry, so P3 can reproduce the NEWS O2
#                              item at each vitals window.
#   vital_lab_extract.parquet  One row per (encounter, measure, measurement time)
#                              for the six vitals and WBC, spanning the whole
#                              encounter (not ward-restricted). P4
#                              (02c_monitoring.R) counts distinct measurement
#                              times per 24 ward-hours from this.
#   ward_times.parquet         The ward intervals, built in 01_cohort.R and
#                              otherwise kept only in memory. P3 and P4 need it to
#                              re-apply the ward restriction after their own
#                              carry-forward / monitoring steps. Written here
#                              exactly as this script consumes it. Columns are
#                              joined_hosp_id, in_dttm, out_dttm: 01 merges
#                              overlapping and touching ADT ward stays into
#                              disjoint blocks per encounter, so every score time
#                              falls inside at most one interval, and drops
#                              hospitalization_id, which a merged block can span.

write_parquet(scores,         here("proj_tables", "scores_components.parquet"))
write_parquet(news_o2_stream, here("proj_tables", "news_o2_stream.parquet"))
write_parquet(ward_times,     here("proj_tables", "ward_times.parquet"))

vital_lab_extract =
  rowbind(
    ftransform(fselect(vitals_list$heart_rate,       joined_hosp_id, time), measure = "heart_rate"),
    ftransform(fselect(vitals_list$respiratory_rate, joined_hosp_id, time), measure = "respiratory_rate"),
    ftransform(fselect(vitals_list$temp_c,           joined_hosp_id, time), measure = "temp_c"),
    ftransform(fselect(vitals_list$sbp,              joined_hosp_id, time), measure = "sbp"),
    ftransform(fselect(vitals_list$spo2,             joined_hosp_id, time), measure = "spo2"),
    ftransform(fselect(vitals_list$gcs,              joined_hosp_id, time), measure = "gcs"),
    ftransform(fselect(fsubset(labs, !is.na(sirs_wbc)), joined_hosp_id, time), measure = "wbc")
  )

write_parquet(vital_lab_extract, here("proj_tables", "vital_lab_extract.parquet"))

message(
  sprintf("✅ Intermediates written | components: %s | o2 stream: %s | extract: %s | ward_times: %s",
          format(nrow(scores),            big.mark = ","),
          format(nrow(news_o2_stream),    big.mark = ","),
          format(nrow(vital_lab_extract), big.mark = ","),
          format(nrow(ward_times),        big.mark = ","))
)

rm(vital_lab_extract)

rm(resp, labs, vitals_list, get_each_vital); gc()

# carryforward (vs 4h, labs 12h) -----------------------------------------------

## set up locf parameters ------------------------------------------------------

score_cols = setdiff(names(scores), c("joined_hosp_id","time"))
lab_cols   = intersect(c("sirs_wbc","sirs_co2","sirs_bands"), score_cols)
vital_cols = setdiff(score_cols, lab_cols)

## function for LOCF within N hours for one column -----------------------------

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

## run locf: vitals 6h, labs 12h -----------------------------------------------

for (cn in vital_cols) scores = locf(scores, cn, hours = 6)
for (cn in lab_cols)   scores = locf(scores, cn, hours = 12)

## news O2 item: 6h as-of carry onto the fixed score grid ----------------------
# Grid-fixed decision (confirmed with Brenna): news_o2 is deliberately NOT in
# vital_cols above, so the score grid built by the Reduce/full-join is identical
# to round one. Here the O2 stream is carried onto each existing score row by a
# rolling last-observation-carried-forward join within the same 6h window used
# for vitals; this roll join *is* the O2 item's LOCF. Because the stream spans
# every respiratory_support timestamp (not just those that happen to coincide
# with a score row), this also picks up O2 measured between score rows, which a
# grid-only LOCF would miss. o2_meas_time records which measurement supplied the
# value so the resolution diagnostic can separate a direct hit (measured at this
# instant) from a carried-forward value.

setDT(scores)
setDT(news_o2_stream)
news_o2_stream[, o2_meas_time := time]
setkey(news_o2_stream, joined_hosp_id, time)

roll_secs = 6 * 3600  # 6h, matching the vital-sign LOCF window

scores = news_o2_stream[scores, on = .(joined_hosp_id, time), roll = roll_secs]
setDT(scores)

scores[, o2_branch := fcase(
  !is.na(o2_meas_time) & o2_meas_time == time & o2_src == 1L, 1L,  # fio2_gt_21
  !is.na(o2_meas_time) & o2_meas_time == time & o2_src == 2L, 2L,  # lpm_gt_0
  !is.na(o2_meas_time) & o2_meas_time == time & o2_src == 3L, 3L,  # room_air
  !is.na(o2_meas_time),                                       4L,  # carried_forward
  default = 5L                                                     # defaulted_zero
)]

scores[, c("o2_src", "o2_meas_time") := NULL]

rm(news_o2_stream); gc()

scores = replace_na(scores, value = 0L, cols = NULL, set = F, type = "const")

rm(score_cols, lab_cols, vital_cols, locf); gc()

# total scores -----------------------------------------------------------------

## sirs preprep ----------------------------------------------------------------

scores$sirs_rr = if_else(scores$sirs_rr == 1 | scores$sirs_co2 == 1, 1L, 0L)
scores = fselect(scores, -sirs_co2, -spo2)
scores = rename(scores, sf = mews_sf)

## all scores ------------------------------------------------------------------

setDT(scores)

scores[, sirs_total  := rowSums(.SD, na.rm = T), .SDcols = patterns("^sirs_")]
scores[, mews_total  := rowSums(.SD, na.rm = T), .SDcols = patterns("^mews_")]
scores[, news_total  := rowSums(.SD, na.rm = T), .SDcols = patterns("^news_")]
scores[, qsofa_total := rowSums(.SD, na.rm = T), .SDcols = patterns("^qsofa_")]

# NEWS single-parameter rule: positive if any component scores 3.
# news_o2 is intentionally excluded: it scores 0 or 2 and can never satisfy a
# single-parameter-3 rule, so news_any3 is unchanged from round one.
scores[, news_any3 := as.integer(
  news_temp == 3L | news_hr == 3L | news_rr == 3L |
    news_sbp == 3L | news_spo2 == 3L | news_gcs == 3L
)]

## QC: NEWS supplemental-oxygen build ------------------------------------------
# news_o2 is 0/2 only; the new news_total is exactly the round-one NEWS total
# plus the O2 item, so it can only be >= the round-one value; news_any3 stays
# 0/1 and does not reference the O2 item. (The definitive byte-identical check
# for SIRS/qSOFA/MEWS/MEWS-SF is the external comparison against commit 37897cf.)

news_total_orig =
  scores$news_temp + scores$news_hr + scores$news_rr +
  scores$news_sbp  + scores$news_spo2 + scores$news_gcs

if (!all(scores$news_o2 %in% c(0L, 2L))) {
  stop("news_o2 took a value other than 0 or 2.", call. = FALSE)
}
if (!all(scores$news_total == news_total_orig + scores$news_o2)) {
  stop("news_total is not the round-one NEWS total plus the O2 item.", call. = FALSE)
}
if (!all(scores$news_total >= news_total_orig)) {
  stop("news_total fell below the round-one NEWS total.", call. = FALSE)
}
if (!all(scores$news_any3 %in% c(0L, 1L))) {
  stop("news_any3 is not 0/1.", call. = FALSE)
}

message("✅ NEWS O2-item QC passed (O2 item additive; news_any3 unchanged).")

rm(news_total_orig)

# scores on the wards only -----------------------------------------------------

scores =
  tidytable(scores) |>
  ftransform(mews_sf_total = mews_total + sf) |>
  join(ward_times, how = "inner", multiple = T) |>
  fsubset(time >= in_dttm & time <= out_dttm) |>
  fgroup_by(joined_hosp_id, in_dttm, out_dttm, time) |>
  fmax() |>
  fgroup_by(joined_hosp_id) |>
  fmutate(
    in_dttm  = fmin(in_dttm),
    out_dttm = fmax(out_dttm)
  ) |>
  fungroup() |>
  funique()

## QC: one score row per encounter and timestamp -------------------------------
# Guaranteed by the disjoint ward blocks built in 01_cohort.R. Asserted rather
# than assumed, because before that merge the funique() above was the only thing
# removing a score row duplicated across two overlapping or touching intervals,
# and it removed it silently.

n_score_keys = fnunique(fselect(scores, joined_hosp_id, time))

if (n_score_keys != nrow(scores)) {
  stop(
    sprintf("Ward restriction produced %d rows for %d distinct encounter-timestamps.",
            nrow(scores), n_score_keys),
    call. = FALSE
  )
}

rm(n_score_keys)

# add outcomes -----------------------------------------------------------------
# P1 rewrote outcome_times to carry, for five outcomes over the competing-
# risk universe {icu, death, hospice}, an indicator (o_k_01), an event time
# (event_k_dttm), and a censoring time (censor_k_dttm). This left join brings
# those through to scores_full.parquet, and for each outcome computes
#   end_k_dttm = pmin(out_dttm, event_k_dttm, censor_k_dttm)
# the last score time that outcome admits; 03a will truncate each outcome at its
# own end_k_dttm. The composite carries no censoring, so end_composite_dttm
# equals the round-one end_dttm exactly and we still truncate on the composite
# here, leaving the exported row set unchanged. The round-one o_primary_01 /
# o_nohospc_01 / end_dttm fields are retained so the existing 03_analysis.R keeps
# running until P5 lands.

## Parquet preserves column types, so every dttm column round-trips as POSIXct --
## including structurally all-missing ones such as censor_composite_dttm. (Under
## the old CSV path those read back as logical and needed a coercion pass; parquet
## removes that hazard entirely, which is the point of the migration.)
outcomes = as.data.table(read_parquet(here("proj_tables", "outcome_times.parquet")))

scores =
  join(scores, outcomes, how = "left", multiple = T) |>
  fmutate(o_primary_01 = if_else(!is.na(outcome_dttm),         1L, 0L)) |>
  fmutate(o_nohospc_01 = if_else(!is.na(outcome_nohospc_dttm), 1L, 0L)) |>
  fmutate(
    end_composite_dttm = pmin(out_dttm, event_composite_dttm, censor_composite_dttm, na.rm = T),
    end_nohospice_dttm = pmin(out_dttm, event_nohospice_dttm, censor_nohospice_dttm, na.rm = T),
    end_wardicu_dttm   = pmin(out_dttm, event_wardicu_dttm,   censor_wardicu_dttm,   na.rm = T),
    end_warddeath_dttm = pmin(out_dttm, event_warddeath_dttm, censor_warddeath_dttm, na.rm = T),
    end_hospicedc_dttm = pmin(out_dttm, event_hospicedc_dttm, censor_hospicedc_dttm, na.rm = T)
  ) |>
  fmutate(end_dttm = pmin(out_dttm, outcome_dttm, na.rm = T)) |>
  fsubset(time < end_dttm)

## QC: outcome carry-through ---------------------------------------------------
# The composite reproduces the round-one truncation exactly and its indicator
# matches the round-one primary indicator.

if (!identical(scores$o_composite_01, scores$o_primary_01)) {
  stop("o_composite_01 (from outcome_times) disagrees with o_primary_01.", call. = FALSE)
}
if (!all(scores$end_composite_dttm == scores$end_dttm)) {
  stop("end_composite_dttm does not equal the round-one end_dttm.", call. = FALSE)
}

message("✅ Outcome carry-through QC passed (composite matches round one).")

rm(outcomes)

# add cancer dx to score df ----------------------------------------------------

cancer =
  fsubset(cohort, ca_01 == 1) |>
  pull(joined_hosp_id)

scores$ca_01 = if_else(scores$joined_hosp_id %in% cancer, 1L, 0L)

# add ed admission to score df -------------------------------------------------

ed =
  fsubset(cohort, ed_admit_01 == 1) |>
  pull(joined_hosp_id)

scores$ed_admit_01 = if_else(scores$joined_hosp_id %in% ed, 1L, 0L)

# news O2 resolution diagnostic ------------------------------------------------
# Counts of score-observation rows by how the O2 item was resolved, by cancer
# status. Under the grid-fixed carry (see above), a concurrent resolution
# (fio2_gt_21 / lpm_gt_0 / room_air) requires an O2 measurement whose timestamp
# coincides exactly with a score row. A value carried from a measurement within
# the preceding 6h is carried_forward: those rows are resolved, merely not
# concurrent, and they account for the large majority of rows (roughly 85% at
# the reference site). A row with no O2 measurement within 6h is defaulted_zero:
# news_o2 is zero-filled and the patient is treated as being on room air. Only
# defaulted_zero represents an absent measurement, and the guard below acts on
# its share of rows.

o2_branch_labs = c("fio2_gt_21", "lpm_gt_0", "room_air", "carried_forward", "defaulted_zero")

news_o2_resolution =
  fcount(scores, ca_01, o2_branch, name = "n") |>
  ftransform(resolution_branch = o2_branch_labs[o2_branch]) |>
  fselect(ca_01, resolution_branch, n) |>
  ftransform(site = site_lowercase) |>
  roworder(ca_01, resolution_branch)

## QC: the diagnostic must account for every score row exactly once, and every
## o2_branch value must fall in 1:5 so that the label lookup resolves.

if (sum(news_o2_resolution$n) != nrow(scores)) {
  stop("news_o2_resolution branch counts do not sum to nrow(scores).", call. = FALSE)
}
if (anyNA(news_o2_resolution$resolution_branch)) {
  stop("resolution_branch contains NA; o2_branch took a value outside 1:5.", call. = FALSE)
}

if (!dir.exists(here(BOX_DIR, "diagnostics"))) {
  dir.create(here(BOX_DIR, "diagnostics"), recursive = TRUE)
}

fwrite(
  news_o2_resolution,
  here(BOX_DIR, "diagnostics", paste0("news_o2_resolution-", site_lowercase, ".csv"))
)

message(
  sprintf("✅ NEWS O2 resolution diagnostic written | %d branch × cancer cells.",
          nrow(news_o2_resolution))
)

## QC: the diagnostic carries ten cells (five branches × two cancer strata)
## unless a branch is empty at this site; name any empty cells rather than
## failing, since an empty branch is legitimate.

o2_expected_cells = expand.grid(
  ca_01             = c(0L, 1L),
  resolution_branch = o2_branch_labs,
  stringsAsFactors  = FALSE
)

o2_expected_keys = paste0(o2_expected_cells$ca_01, "|", o2_expected_cells$resolution_branch)
o2_observed_keys = paste0(news_o2_resolution$ca_01, "|", news_o2_resolution$resolution_branch)
o2_missing_keys  = setdiff(o2_expected_keys, o2_observed_keys)

if (length(o2_missing_keys) == 0L) {
  message("✅ NEWS O2 diagnostic complete | 10 of 10 branch × cancer cells populated.")
} else {
  message(sprintf(
    "  note: %d of 10 branch × cancer cells are empty at this site (ca_01|branch): %s",
    length(o2_missing_keys), paste(o2_missing_keys, collapse = ", ")
  ))
}

# sparse-oxygen guard ----------------------------------------------------------
# 00_setup.R fails a site whose respiratory_support lacks lpm_set or fio2_set. A
# site where those columns exist but are entirely null passes that check, sends
# every score row to the defaulted_zero branch, and runs NEWS without its
# supplemental-oxygen item without any downstream complaint. The defaulted
# fraction is therefore printed on every run, warned on above 0.50, and treated
# as fatal above 0.90. Reference site: 17.6% of non-cancer rows and 11.9% of
# cancer rows.

n_rows_all = nrow(scores)
n_rows_ca  = sum(scores$ca_01 == 1L)
n_rows_no  = sum(scores$ca_01 == 0L)

n_defaulted_all = sum(scores$o2_branch == 5L)
n_defaulted_ca  = sum(scores$o2_branch == 5L & scores$ca_01 == 1L)
n_defaulted_no  = sum(scores$o2_branch == 5L & scores$ca_01 == 0L)

frac_defaulted_all = if (n_rows_all > 0L) n_defaulted_all / n_rows_all else NA_real_
frac_defaulted_ca  = if (n_rows_ca  > 0L) n_defaulted_ca  / n_rows_ca  else NA_real_
frac_defaulted_no  = if (n_rows_no  > 0L) n_defaulted_no  / n_rows_no  else NA_real_

message(sprintf(
  "  NEWS O2 defaulted to room air | overall %.1f%% | non-cancer %.1f%% | cancer %.1f%%",
  100 * frac_defaulted_all, 100 * frac_defaulted_no, 100 * frac_defaulted_ca
))

o2_sparse_msg = paste(
  "",
  "================================================================================",
  sprintf("  SPARSE SUPPLEMENTAL-OXYGEN DATA: %.1f%% of score rows have no oxygen",
          100 * frac_defaulted_all),
  "  measurement within 6 hours and were scored as room air. NEWS is largely",
  "  running without its supplemental-oxygen item at this site.",
  "",
  "  Check lpm_set, fio2_set, and device_category in respiratory_support. The",
  "  columns pass the 00_setup.R presence check but may be empty or unmapped.",
  "",
  "  Contact the coordinating center before uploading these artifacts.",
  "================================================================================",
  "",
  sep = "\n"
)

if (isTRUE(frac_defaulted_all > 0.50)) {
  warning(o2_sparse_msg, call. = FALSE, immediate. = TRUE)
}

if (isTRUE(frac_defaulted_all > 0.90) && !allow_sparse_o2) {
  stop(
    paste0(
      o2_sparse_msg,
      "  The oxygen fields are effectively empty rather than merely sparse, so the\n",
      "  run has been stopped. Set allow_sparse_o2: true in\n",
      "  config/config_clif_oncrisk.yaml only after the coordinating center confirms\n",
      "  that this site's oxygen data are genuinely absent.\n"
    ),
    call. = FALSE
  )
}

## o2_branch is a per-row helper for the diagnostic only; drop it so
## scores_full.parquet gains only news_o2 and the outcome columns.
scores = fselect(scores, -o2_branch)

rm(cancer, ed, news_o2_resolution, o2_branch_labs,
   o2_expected_cells, o2_expected_keys, o2_observed_keys, o2_missing_keys,
   n_rows_all, n_rows_ca, n_rows_no,
   n_defaulted_all, n_defaulted_ca, n_defaulted_no,
   frac_defaulted_all, frac_defaulted_ca, frac_defaulted_no, o2_sparse_msg)

# reconcile cohort against the score table -------------------------------------
# A small number of encounters clear every 01_cohort.R exclusion and still
# contribute no score rows. The "< 6h data" filter asks whether a qualifying
# vital exists at least 6h after first ward arrival; it does not require that
# vital to precede the outcome. An encounter whose composite event lands between
# hour 6 and its next qualifying vital therefore survives the cohort filters and
# then loses every row to the `time < end_dttm` truncation above.
#
# Rather than predict which encounters those are, drop them by observation: any
# cohort encounter with zero score rows is removed from cohort.parquet and
# recorded as a sixth step in the inclusion flow. This also catches any other
# route to an empty score set. Because all five end_*_dttm resolve to the
# earliest of ICU / ward death / hospice, an encounter with zero rows under the
# composite has zero rows under every outcome, so one pass covers all of them.

scored_jids = funique(scores$joined_hosp_id)
orphan_jids = setdiff(cohort$joined_hosp_id, scored_jids)

orphan_ca = fsubset(cohort, joined_hosp_id %in% orphan_jids & ca_01 == 1L)
orphan_no = fsubset(cohort, joined_hosp_id %in% orphan_jids & ca_01 == 0L)

message(sprintf(
  "  cohort/score reconciliation: %d encounter(s) with no calculable score (%d cancer, %d non-cancer).",
  length(orphan_jids), nrow(orphan_ca), nrow(orphan_no)
))

cohort = fsubset(cohort, !joined_hosp_id %in% orphan_jids)

write_parquet(cohort, here("proj_tables", "cohort.parquet"))

## append the reconciliation step to the inclusion flow ------------------------
# Rows 2-5 of the flow count ED-admit encounters only, so row 6 does too, even
# though the drop above applies to the full cohort (the broader cohort feeds the
# se_no_ed_req variant).

flow_path = here(BOX_DIR, paste0("figure_s01_flow_", site_lowercase, ".csv"))
flow_df   = fread(flow_path)

ed_reconciled = fsubset(cohort, ed_admit_01 == 1L)
n_remain_ca   = fnunique(fsubset(ed_reconciled, ca_01 == 1L)$joined_hosp_id)
n_remain_no   = fnunique(fsubset(ed_reconciled, ca_01 == 0L)$joined_hosp_id)

flow_row = tidytable(
  step           = "After excluding encounters with no calculable score before the outcome",
  n_remaining_ca = n_remain_ca,
  n_excluded_ca  = flow_df$n_remaining_ca[nrow(flow_df)] - n_remain_ca,
  n_remaining_no = n_remain_no,
  n_excluded_no  = flow_df$n_remaining_no[nrow(flow_df)] - n_remain_no
)

flow_df = rbind(flow_df, flow_row)

fwrite(flow_df, flow_path)

## QC --------------------------------------------------------------------------
# The cohort and the score table must now name exactly the same encounters, in
# both directions, and agree on the ED-admit denominator that every downstream
# artifact family reports.

stopifnot(length(setdiff(cohort$joined_hosp_id, scored_jids)) == 0L)
stopifnot(length(setdiff(scored_jids, cohort$joined_hosp_id)) == 0L)

n_ed_cohort = fnunique(ed_reconciled$joined_hosp_id)
n_ed_scores = fnunique(fsubset(scores, ed_admit_01 == 1L)$joined_hosp_id)

stopifnot(n_ed_cohort == n_ed_scores)

message(sprintf(
  "✅ cohort and score table reconciled | %s encounters, %s ED-admit.",
  format(nrow(cohort), big.mark = ","), format(n_ed_cohort, big.mark = ",")
))

## set the run-log counts reported by run_all ----------------------------------
# The composite outcome rate is set further down, in the descriptive block, where
# its numerator is already computed from the same reconciled sample.

run_log$n_cohort       = nrow(cohort)
run_log$n_cancer       = sum(cohort$ca_01 == 1L)
run_log$n_ed_admit     = n_ed_cohort
run_log$n_no_score_row = length(orphan_jids)

rm(scored_jids, orphan_jids, orphan_ca, orphan_no, flow_path, flow_df,
   flow_row, ed_reconciled, n_remain_ca, n_remain_no, n_ed_cohort, n_ed_scores)

# encounter-level descriptive artifacts ----------------------------------------
# Relocated from 01_cohort.R so that every denominator below is the reconciled
# ED-admit cohort -- the encounters surviving both the 01 exclusions and the
# no-calculable-score drop above -- and therefore matches the denominator every
# analysis artifact in 03a/03b reports. Built in 01, these five artifacts ran on
# the pre-reconciliation set and disagreed with the rest of the pipeline.
#
# Second change, at the coordinating center's direction: the encounter-level
# outcome rows of table_02_cat are now the competing-risk indicators from
# outcome_times.parquet (composite, nohospice, wardicu, warddeath, hospicedc),
# under the same outcome keys used in every analysis filename. The round-one
# rows wicu / d_noicu / h_noicu are removed. Those three conditioned on icu_01
# (ever in an ICU, by any route) whereas the analysis conditions on a ward-to-ICU
# transfer, so a ward -> procedural -> ICU path separated them and Table 2
# implied a composite roughly two percent below the one the analysis reports.
# icu_01, dead_01, and hospice_01 stay: they are unambiguous descriptive facts
# about the encounter and carry no outcome definition.

t2_cohort = fsubset(cohort, ed_admit_01 == 1L)
n_t2      = fnobs(t2_cohort$joined_hosp_id)

## attach the competing-risk outcome indicators --------------------------------

outcome_flags =
  read_parquet(here("proj_tables", "outcome_times.parquet")) |>
  as_tidytable() |>
  fselect(
    joined_hosp_id,
    composite_01 = o_composite_01,
    nohospice_01 = o_nohospice_01,
    wardicu_01   = o_wardicu_01,
    warddeath_01 = o_warddeath_01,
    hospicedc_01 = o_hospicedc_01
  )

t2_flag_cols = c(
  "composite_01",
  "nohospice_01",
  "wardicu_01",
  "warddeath_01",
  "hospicedc_01"
)

t2_cohort =
  join(t2_cohort, outcome_flags, how = "left", multiple = FALSE) |>
  mutate(across(
    .cols = all_of(t2_flag_cols),
    .fns  = ~if_else(is.na(.x), 0L, .x)
  ))

## QC: the outcome rows must decompose exactly ---------------------------------
# Same identity 03a asserts across artifacts. An exact timestamp tie between two
# competing events would break it, so failing here rather than in 03a locates the
# problem at the site's data instead of at the artifact matrix.

n_composite = fsum(t2_cohort$composite_01)
n_nohospice = fsum(t2_cohort$nohospice_01)
n_wardicu   = fsum(t2_cohort$wardicu_01)
n_warddeath = fsum(t2_cohort$warddeath_01)
n_hospicedc = fsum(t2_cohort$hospicedc_01)

if (n_composite != n_wardicu + n_warddeath + n_hospicedc) {
  stop(
    sprintf("Table 2 composite (%d) != wardicu + warddeath + hospicedc (%d).",
            n_composite, n_wardicu + n_warddeath + n_hospicedc),
    call. = FALSE
  )
}

if (n_nohospice != n_wardicu + n_warddeath) {
  stop(
    sprintf("Table 2 nohospice (%d) != wardicu + warddeath (%d).",
            n_nohospice, n_wardicu + n_warddeath),
    call. = FALSE
  )
}

message(
  sprintf("✅ Table 2 outcome decomposition | composite %s = wardicu %s + warddeath %s + hospicedc %s",
          format(n_composite, big.mark = ","), format(n_wardicu, big.mark = ","),
          format(n_warddeath, big.mark = ","), format(n_hospicedc, big.mark = ","))
)

## composite outcome rate for the run report -----------------------------------
# Set here rather than in run_all.R after stage 01. There it was computed over
# the pre-reconciliation ED-admit set while the three counts printed beside it in
# run_report were refreshed post-reconciliation, so a site with a larger
# no-calculable-score drop would have exported a mismatched pair. The numerator
# and denominator below are the ones table_02_cat and table_02_cont export.

run_log$outcome_rate = n_composite / n_t2

# missingness characterization -------------------------------------------------

miss_summary =
  fsummarize(
    t2_cohort,
    n_total     = fnobs(joined_hosp_id),
    n_miss_age  = fsum(miss_age),
    n_miss_sex  = fsum(miss_sex),
    n_miss_race = fsum(miss_race),
    n_miss_eth  = fsum(miss_eth),
    n_miss_vw   = fsum(miss_vw)
  )

n_denom = miss_summary$n_total

miss_summary =
  fselect(miss_summary, -n_total) |>
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") |>
  ftransform(
    variable    = str_remove(variable, "^n_miss_"),
    n_total     = n_denom,
    pct_missing = round(100 * n_missing / n_denom, 2),
    site        = site_lowercase
  )

fwrite(miss_summary, here(BOX_DIR, paste0("missing_demog_", site_lowercase, ".csv")))

# t02. characteristics/outcomes by cancer status (0 = none, 1 = cancer) --------

## prepare table component = continuous variables ------------------------------

t2_cont =
  fgroup_by(t2_cohort, ca_01) |>
  fsummarize(
    n         = fnobs(joined_hosp_id),
    age_sum   = fsum(age),
    age_sumsq = fsum(age^2),
    age_p025  = fquantile(age, 0.025),
    age_p975  = fquantile(age, 0.975),
    vw_sum    = fsum(vw),
    vw_sumsq  = fsum(vw^2),
    los_sum   = fsum(los_hosp_d),
    los_sumsq = fsum(los_hosp_d^2),
    los_p025  = fquantile(los_hosp_d, 0.025),
    los_p975  = fquantile(los_hosp_d, 0.975)
  ) |>
  ftransform(var_type = "continuous") |>
  ftransform(site     = paste0(site_lowercase))

## prepare table component = categorized continuous variables ------------------

### age ------------------------------------------------------------------------
# Bound into t2_cat below. Round one computed these bands and then never bound
# them, so table_02_cat carried vw and los bands but no age bands.

age_breaks = c( 18,      40,      50,      60,      70,      80,  Inf)
age_labs   = c("18_39", "40_49", "50_59", "60_69", "70_79", "80_plus")

ages_cat =
  fselect(t2_cohort, ca_01, age) |>
  fmutate(a = cut(age, breaks = age_breaks, labels = age_labs, right = F)) |>
  fgroup_by(ca_01, a) |>
  fnobs() |>
  select(ca_01, age_cat = a, n = age) |>
  ftransform(age_cat = paste0("age_", age_cat)) |>
  ftransform(var = "age", category = str_remove(age_cat, "age_")) |>
  fselect(ca_01, var, category, n)

### elixhauser -----------------------------------------------------------------

elix_breaks = c(-Inf,  0, 4,  9, 14,  Inf)
elix_labs   = c("<= 0", "1-4", "5-9", "10-14", ">= 15")

elix_cat =
  fselect(t2_cohort, ca_01, vw) |>
  fmutate(a = cut(vw, breaks = elix_breaks, labels = elix_labs, right = F)) |>
  fgroup_by(ca_01, a) |>
  fnobs() |>
  select(ca_01, elix_cat = a, n = vw) |>
  ftransform(elix_cat = paste0("vw_", elix_cat)) |>
  ftransform(var = "vw", category = str_remove(elix_cat, "vw_")) |>
  fselect(ca_01, var, category, n)

### los (days) -----------------------------------------------------------------

l_breaks = c( 0,        2,         4,         7,         14, Inf)
l_labs   = c("0-47h", "48h_96h", "96h_1wk", "1wk_2wk", "2wk_plus")

los_cat =
  fselect(t2_cohort, ca_01, los_hosp_d) |>
  fmutate(los_cat = cut(los_hosp_d, breaks = l_breaks, labels = l_labs, right = F)) |>
  fgroup_by(ca_01, los_cat) |>
  fnobs() |>
  select(ca_01, los_cat, n = los_hosp_d) |>
  ftransform(var = "los", category = los_cat) |>
  fselect(ca_01, var, category, n)

## prepare table component = categorical variables -----------------------------
# Held out deliberately: rank_enc and mets_01 (carried in cohort.parquet, and
# mets_01 is three-valued so a pivot would drop its missing level silently); the
# miss_* indicators (exported above as missing_demog); and wicu_01 / d_noicu_01 /
# h_noicu_01, superseded by the competing-risk indicators attached above.

t2_cat =
  select(t2_cohort,
         -ends_with("id"),
         -ends_with("dttm"),
         -age,
         -vw,
         -los_hosp_d,
         -rank_enc,
         -mets_01,
         -starts_with("miss_"),
         -wicu_01,
         -d_noicu_01,
         -h_noicu_01) |>
  pivot_longer(-ca_01, names_to = "var", values_to = "val") |>
  fsubset(!is.na(val)) |>
  fmutate(n = val) |>
  fgroup_by(ca_01, var, val) |>
  fnobs() |>
  ftransform(var      = str_remove(var, "_01")) |>
  ftransform(category = tolower(str_replace_all(as.character(val), "-", "_"))) |>
  fselect(ca_01, var, category, n) |>
  rowbind(ages_cat) |>
  rowbind(elix_cat) |>
  rowbind(los_cat)

## collapse rare race categories -----------------------------------------------
# Site-level cells for the rarest race categories can be small enough to be
# identifying. American Indian or Alaska Native and Native Hawaiian or Other
# Pacific Islander are folded into the existing other/unknown level at every
# site (three rows become one), so the pooled table carries one consistent
# category set. Alternate spellings of the residual level are folded too.

race_other = c(
  "american indian or alaska native",
  "native hawaiian or other pacific islander",
  "other/unknown",
  "other_unknown",
  "other",
  "unknown",
  "two or more races"
)

t2_cat =
  ftransform(
    t2_cat,
    category = if_else(
      var == "race_category" & category %in% race_other,
      "other/unknown",
      category
    )
  ) |>
  fgroup_by(ca_01, var, category) |>
  fsummarize(n = fsum(n)) |>
  fungroup()

## quality control -------------------------------------------------------------

### check for sample size mismatch in continuous table -------------------------

if (sum(t2_cont$n) != n_t2) {
  stop(
    sprintf("Table 2 sample size mismatch: t2_cont sums to %d, reconciled ED-admit cohort is %d.",
            sum(t2_cont$n), n_t2),
    call. = FALSE
  )
}

### the composite row must equal the analysis denominator's event count --------

t2_composite_n = fsum(fsubset(t2_cat, var == "composite" & category == "1")$n)

if (t2_composite_n != n_composite) {
  stop(
    sprintf("Table 2 composite row (%d) disagrees with the encounter-level flag (%d).",
            t2_composite_n, n_composite),
    call. = FALSE
  )
}

message(sprintf("✅ Table 2 QC passed | %s ED-admit encounters.", format(n_t2, big.mark = ",")))

## export table 2 --------------------------------------------------------------
# Small-cell suppression: every categorical count below 10 is exported as the
# string "<10". The pooled loader (coerce_suppressed_counts) parses this marker
# and imputes a value in [0, 9]; rows are never dropped, so pooled denominators
# remain reconstructible. Applied at export only, AFTER the QC above, which
# needs exact counts.

n_masked = sum(t2_cat$n < 10)

t2_cat_export = ftransform(t2_cat, n = if_else(n < 10, "<10", as.character(n)))

message(sprintf("  table_02_cat: %d cell(s) below 10 masked as \"<10\".", n_masked))

fwrite(t2_cat_export, here(BOX_DIR, paste0("table_02_cat_",  site_lowercase, ".csv")))
fwrite(t2_cont,       here(BOX_DIR, paste0("table_02_cont_", site_lowercase, ".csv")))

# admission diagnosis (reason for admission) -----------------------------------
# Round-two export. Three-character ICD-10-CM stems are exported unmapped; the
# coordinating center maps them to CCSR centrally, which avoids distributing a
# crosswalk to eight sites. Denominators are recoverable from table_02_cont,
# whose n by ca_01 is the same reconciled ed_admit_01 == 1 sample used here.

## does this site populate diagnosis_primary? ----------------------------------
# 00_setup.R does not validate diagnosis_primary, so a site can reach this point
# without it. A descriptive add-on should not halt the primary analysis, so an
# absent or never-positive field yields header-only artifacts and a message.

has_dx_primary = "diagnosis_primary" %in% names(data_list$hospital_diagnosis)

if (!has_dx_primary) {

  message("Note: hospital_diagnosis has no diagnosis_primary column; admission-diagnosis artifacts will be empty.")

  adm_dx_raw = tidytable(
    hospitalization_id = hid_jid_crosswalk$hospitalization_id[0],
    diagnosis_code     = character(0)
  )

} else {

  adm_dx_raw =
    dplyr::filter(data_list$hospital_diagnosis, hospitalization_id %in% cohort_hids) |>
    dplyr::filter(toupper(diagnosis_code_format) == "ICD10CM") |>
    dplyr::select(hospitalization_id, diagnosis_code, diagnosis_primary) |>
    dplyr::collect()

  ### diagnosis_primary arrives as integer, logical, or character across sites
  adm_dx_raw =
    ftransform(adm_dx_raw, dx_primary_int = suppressWarnings(as.integer(diagnosis_primary))) |>
    fsubset(!is.na(dx_primary_int) & dx_primary_int == 1L) |>
    fselect(hospitalization_id, diagnosis_code) |>
    distinct()

  if (nrow(adm_dx_raw) == 0) {
    message("Note: no rows with diagnosis_primary == 1; admission-diagnosis artifacts will be empty.")
  }
}

## first primary code per linked encounter -------------------------------------
# A joined_hosp_id can span several source hospitalizations, each with its own
# primary code. The script's established rule for collapsing constituents is to
# order by time and take the first, so the primary code of the earliest
# constituent admission wins, with ascending code as the deterministic tie-break.

hosp_order =
  dplyr::filter(data_list$hospitalization, hospitalization_id %in% cohort_hids) |>
  dplyr::select(hospitalization_id, admission_dttm) |>
  dplyr::collect() |>
  distinct()

adm_dx_first =
  join(adm_dx_raw, hid_jid_crosswalk, how = "inner", multiple = TRUE) |>
  join(hosp_order, how = "inner", multiple = FALSE) |>
  ftransform(code_clean = toupper(str_remove_all(diagnosis_code, "[^A-Za-z0-9]"))) |>
  fsubset(str_detect(code_clean, "^[A-Z]")) |>
  roworder(admission_dttm, code_clean) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(code_clean = ffirst(code_clean))

## attach cancer status and derive chapter letter and code stem ----------------

adm_dx_enc =
  fselect(t2_cohort, joined_hosp_id, ca_01) |>
  join(adm_dx_first, how = "inner", multiple = FALSE) |>
  ftransform(chapter   = substr(code_clean, 1, 1)) |>
  ftransform(code_stem = substr(code_clean, 1, 3))

## tally at both levels, suppressing cells of five or fewer --------------------

adm_dx_chapter =
  fselect(adm_dx_enc, ca_01, chapter, joined_hosp_id) |>
  fgroup_by(ca_01, chapter) |>
  fnobs() |>
  fselect(ca_01, chapter, n = joined_hosp_id) |>
  fsubset(n > 5) |>
  roworder(ca_01, -n, chapter) |>
  ftransform(site = site_lowercase)

adm_dx_stem =
  fselect(adm_dx_enc, ca_01, code_stem, joined_hosp_id) |>
  fgroup_by(ca_01, code_stem) |>
  fnobs() |>
  fselect(ca_01, code_stem, n = joined_hosp_id) |>
  fsubset(n > 5) |>
  roworder(ca_01, -n, code_stem) |>
  ftransform(site = site_lowercase)

fwrite(
  adm_dx_chapter,
  here(BOX_DIR, paste0("admission_dx_chapter-ca-", site_lowercase, ".csv"))
)

fwrite(
  adm_dx_stem,
  here(BOX_DIR, paste0("admission_dx_stem-ca-", site_lowercase, ".csv"))
)

message(
  sprintf("✅ Admission diagnoses | %s of %s encounters coded | %d chapters | %d stems",
          format(nrow(adm_dx_enc), big.mark = ","),
          format(n_t2, big.mark = ","),
          nrow(adm_dx_chapter),
          nrow(adm_dx_stem))
)

rm(t2_cohort, n_t2, outcome_flags, t2_flag_cols,
   n_composite, n_nohospice, n_wardicu, n_warddeath, n_hospicedc,
   miss_summary, n_denom,
   t2_cont, age_breaks, age_labs, ages_cat, elix_breaks, elix_labs, elix_cat,
   l_breaks, l_labs, los_cat, t2_cat, t2_cat_export, race_other, n_masked,
   t2_composite_n,
   has_dx_primary, adm_dx_raw, hosp_order, adm_dx_first, adm_dx_enc,
   adm_dx_chapter, adm_dx_stem); gc()


# save scores ------------------------------------------------------------------

.n_score_rows = nrow(scores)

write_parquet(scores, here("proj_tables", "scores_full.parquet"))
write_parquet(data.frame(site_lowercase = site_lowercase, stringsAsFactors = FALSE),
              here("proj_tables", "site_lowercase.parquet"))

rm(list = setdiff(
  ls(),
  c("BOX_DIR", "start_time", "run_log", "site_lowercase", "here", ".n_score_rows")
))

gc()

# go to 03

################################################################################
