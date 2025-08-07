
# ====================== Params ================================================

cats             = c("heart_rate","respiratory_rate","sbp","dbp","temp_c","spo2")
union_round_min  = 15L      # set NULL to skip rounding
vital_window_min = 240L
resp_window_min  = 60L
lab_window_wbc   = 720L
lab_window_bands = 720L
lab_window_paco2 = 240L
tz_use           = "UTC"
fio2_units       = "auto"   # "auto" | "fraction" | "percent"
out_parquet      = "scores_ews.parquet"

# ====================== Vitals: pull, coerce, de-outlier ======================

vs = 
  dplyr::filter(data_list$vitals, vital_category %in% cats) |>
  dplyr::select(hospitalization_id, recorded_dttm, vital_category, vital_value) |>
  dplyr::collect()

setDT(vs)
vs[, `:=`(
  recorded_dttm = as.POSIXct(recorded_dttm, tz = tz_use),
  vital_value   = suppressWarnings(as.numeric(vital_value))
)]

# Drop implausible values BEFORE LOCF
vs[vital_category == "heart_rate",       vital_value := fifelse(vital_value < 20 | vital_value > 250, NA_real_, vital_value)]
vs[vital_category == "respiratory_rate", vital_value := fifelse(vital_value < 5  | vital_value > 60,  NA_real_, vital_value)]
vs[vital_category == "sbp",              vital_value := fifelse(vital_value < 50 | vital_value > 300, NA_real_, vital_value)]
vs[vital_category == "dbp",              vital_value := fifelse(vital_value < 20 | vital_value > 200, NA_real_, vital_value)]
vs[vital_category == "temp_c",           vital_value := fifelse(vital_value < 25 | vital_value > 45,  NA_real_, vital_value)]
vs[vital_category == "spo2",             vital_value := fifelse(vital_value < 50 | vital_value > 100, NA_real_, vital_value)]

# ====================== Build union (optionally rounded) grid =================

grid = unique(vs[, .(hospitalization_id, recorded_dttm)])
if (!is.null(union_round_min)) {
  step = union_round_min * 60
  grid[, recorded_dttm := as.POSIXct((as.numeric(recorded_dttm) %/% step) * step,
                                     origin = "1970-01-01", tz = tz_use)]
  setkey(grid, hospitalization_id, recorded_dttm)
  grid = unique(grid)
}
setkey(grid, hospitalization_id, recorded_dttm)

# ====================== Roll each vital once ==================================

cat_map = c(heart_rate="hr", respiratory_rate="rr", sbp="sbp",
             dbp="dbp", temp_c="temp_c", spo2="spo2")

for (c_in in names(cat_map)) {
  c_out = cat_map[[c_in]]
  x = vs[vital_category == c_in, .(hospitalization_id, recorded_dttm, vital_value)]
  setkey(x, hospitalization_id, recorded_dttm)
  j = x[grid, on = .(hospitalization_id, recorded_dttm), roll = TRUE]
  age_min = as.numeric(difftime(grid$recorded_dttm, j$recorded_dttm, units = "mins"))
  grid[, (c_out) := ifelse(is.na(age_min) | age_min > vital_window_min, NA_real_, j$vital_value)]
}
rm(vs, x, j); gc()

DT = grid[]; rm(grid); gc()

# (Optional second-pass plausibility guard after LOCF)
DT[, `:=`(
  hr     = ifelse(hr     < 20 | hr     > 250, NA_real_, hr),
  rr     = ifelse(rr     < 5  | rr     > 60,  NA_real_, rr),
  sbp    = ifelse(sbp    < 50 | sbp    > 300, NA_real_, sbp),
  dbp    = ifelse(dbp    < 20 | dbp    > 200, NA_real_, dbp),
  temp_c = ifelse(temp_c < 25 | temp_c > 45,  NA_real_, temp_c),
  spo2   = ifelse(spo2   < 50 | spo2   > 100, NA_real_, spo2)
)]

# ====================== Component points (vital-only) =========================

DT[, `:=`(
  # qSOFA
  qsofa_rr_i  = as.integer(!is.na(rr)  & rr  >= 22),
  qsofa_sbp_i = as.integer(!is.na(sbp) & sbp <= 100),
  
  # MEWS
  mews_rr_i   = fcase(is.na(rr), NA_integer_,
                      rr <= 8, 3L, rr <= 14, 0L, rr <= 20, 1L, rr <= 29, 2L, default = 3L),
  mews_hr_i   = fcase(is.na(hr), NA_integer_,
                      hr <= 40, 2L, hr <= 50, 1L, hr <= 100, 0L, hr <= 110, 1L,
                      hr <= 129, 2L, default = 3L),
  mews_sbp_i  = fcase(is.na(sbp), NA_integer_,
                      sbp <= 70, 3L, sbp <= 80, 2L, sbp <= 100, 1L, sbp <= 199, 0L, default = 2L),
  mews_temp_i = fcase(is.na(temp_c), NA_integer_,
                      temp_c <= 35, 2L, temp_c < 38.5, 0L, default = 1L),
  
  # NEWS2 (Scale 1 SpO2)
  news_rr_i   = fcase(is.na(rr), NA_integer_,
                      rr <= 8, 3L, rr <= 11, 1L, rr <= 20, 0L, rr <= 24, 2L, default = 3L),
  news_spo2_i = fcase(is.na(spo2), NA_integer_,
                      spo2 <= 91, 3L, spo2 <= 93, 2L, spo2 <= 95, 1L, spo2 >= 96, 0L, default = NA_integer_),
  news_temp_i = fcase(is.na(temp_c), NA_integer_,
                      temp_c <= 35, 3L, temp_c <= 36, 1L, temp_c <= 38, 0L, temp_c <= 39, 1L, default = 2L),
  news_sbp_i  = fcase(is.na(sbp), NA_integer_,
                      sbp <= 90, 3L, sbp <= 100, 2L, sbp <= 110, 1L, sbp <= 219, 0L, default = 3L),
  news_hr_i   = fcase(is.na(hr), NA_integer_,
                      hr <= 40, 3L, hr <= 50, 1L, hr <= 90, 0L, hr <= 110, 1L, hr <= 130, 2L, default = 3L),
  
  # SIRS vital components
  sirs_temp_i = as.integer(!is.na(temp_c) & (temp_c < 36 | temp_c > 38)),
  sirs_hr_i   = as.integer(!is.na(hr)     & hr > 90),
  sirs_rr_i   = as.integer(!is.na(rr)     & rr >= 20),
  
  # Placeholders until you attach mentation
  qsofa_ment_i = NA_integer_,
  mews_avpu_i  = NA_integer_,
  news_conf_i  = NA_integer_
)]

# ====================== Respiratory support (FiO2/flow) =======================

rs = data_list$respiratory_support |>
  select(hospitalization_id, recorded_dttm, lpm_set, fio2_set) |>
  collect()
setDT(rs)

if (nrow(rs)) {
  rs[, `:=`(
    recorded_dttm = as.POSIXct(recorded_dttm, tz = tz_use),
    lpm_set  = suppressWarnings(as.numeric(lpm_set)),
    fio2_set = suppressWarnings(as.numeric(fio2_set))
  )]
  rs[, lpm_set := fifelse(lpm_set < 0 | lpm_set > 200, NA_real_, lpm_set)]
  setkey(rs, hospitalization_id, recorded_dttm)
  
  rj = rs[DT, on = .(hospitalization_id, recorded_dttm), roll = TRUE]
  age_min = as.numeric(difftime(DT$recorded_dttm, rj$recorded_dttm, units = "mins"))
  lpm  = ifelse(is.na(age_min) | age_min > resp_window_min, NA_real_, rj$lpm_set)
  fio2 = ifelse(is.na(age_min) | age_min > resp_window_min, NA_real_, rj$fio2_set)
  
  # Normalize FiO2 → fraction
  fio2_frac = if (fio2_units == "fraction") fio2 else if (fio2_units == "percent") fio2/100 else {
    med = suppressWarnings(median(fio2, na.rm = TRUE))
    if (is.finite(med) && med <= 1) fio2 else fio2/100
  }
  fio2_frac = ifelse(is.na(fio2_frac), NA_real_, pmin(pmax(fio2_frac, 0.21), 1.00))
  on_o2     = as.integer((!is.na(lpm) & lpm > 0) | (!is.na(fio2_frac) & fio2_frac > 0.21 + 1e-9))
  news_o2_i = fifelse(is.na(on_o2), NA_integer_, ifelse(on_o2 == 1L, 2L, 0L))
  sf_ratio  = ifelse(!is.na(DT$spo2) & !is.na(fio2_frac), DT$spo2 / fio2_frac, NA_real_)
  
  DT[, `:=`(lpm_set = lpm, fio2_raw = fio2, fio2_frac = fio2_frac,
            on_o2 = on_o2, news_o2_i = news_o2_i, sf_ratio = sf_ratio)]
}
rm(rs, rj, age_min, lpm, fio2, fio2_frac, on_o2, news_o2_i, sf_ratio); gc()

# ====================== MEWS-SF (bins: >315=0, 236–315=2, ≤235=3) =============

mews_sf_pts = rep(NA_integer_, nrow(DT))
ok = !is.na(DT$sf_ratio)
if (any(ok)) {
  bins = cut(DT$sf_ratio[ok], breaks = c(-Inf, 235, 315, Inf),
              right = TRUE, include.lowest = TRUE, labels = FALSE)
  mews_sf_pts[ok] = c(3L, 2L, 0L)[bins]
}
DT[, mews_sf_i := mews_sf_pts]
rm(mews_sf_pts, ok, bins); gc()

# ====================== Labs: WBC, Bands, PaCO2 ===============================

labs = data_list$labs |>
  select(hospitalization_id, lab_result_dttm, lab_category, lab_value_numeric) |>
  collect()
setDT(labs)

if (nrow(labs)) {
  labs[, `:=`(
    lab_result_dttm  = as.POSIXct(lab_result_dttm, tz = tz_use),
    lab_value_numeric = suppressWarnings(as.numeric(lab_value_numeric))
  )]
  
  # WBC
  wbc_tbl = labs[lab_category == "wbc"]
  if (nrow(wbc_tbl)) {
    setkey(wbc_tbl, hospitalization_id, lab_result_dttm)
    j = wbc_tbl[DT, on = .(hospitalization_id, lab_result_dttm = recorded_dttm), roll = TRUE]
    age = as.numeric(difftime(DT$recorded_dttm, j$lab_result_dttm, units = "mins"))
    wbc_vals = ifelse(is.na(age) | age > lab_window_wbc, NA_real_, j$lab_value_numeric)
    wbc_vals = ifelse(wbc_vals < 0.1 | wbc_vals > 200, NA_real_, wbc_vals)
    DT[, `:=`(wbc = wbc_vals, sirs_wbc_i = as.integer(!is.na(wbc_vals) & (wbc_vals < 4 | wbc_vals > 12)))]
  } else DT[, `:=`(wbc = NA_real_, sirs_wbc_i = NA_integer_)]
  
  # Bands
  bands_tbl = labs[lab_category %in% c("bands","band_neutrophils","immature_neutrophils")]
  if (nrow(bands_tbl)) {
    setkey(bands_tbl, hospitalization_id, lab_result_dttm)
    j = bands_tbl[DT, on = .(hospitalization_id, lab_result_dttm = recorded_dttm), roll = TRUE]
    age = as.numeric(difftime(DT$recorded_dttm, j$lab_result_dttm, units = "mins"))
    bands_vals = ifelse(is.na(age) | age > lab_window_bands, NA_real_, j$lab_value_numeric)
    bands_vals = ifelse(bands_vals < 0 | bands_vals > 100, NA_real_, bands_vals)
    DT[, `:=`(bands = bands_vals, sirs_bands_i = as.integer(!is.na(bands_vals) & bands_vals > 10))]
  } else DT[, `:=`(bands = NA_real_, sirs_bands_i = NA_integer_)]
  
  # PaCO2
  pco_tbl = labs[lab_category == "paco2"]
  if (nrow(pco_tbl)) {
    setkey(pco_tbl, hospitalization_id, lab_result_dttm)
    j = pco_tbl[DT, on = .(hospitalization_id, lab_result_dttm = recorded_dttm), roll = TRUE]
    age = as.numeric(difftime(DT$recorded_dttm, j$lab_result_dttm, units = "mins"))
    pco_vals = ifelse(is.na(age) | age > lab_window_paco2, NA_real_, j$lab_value_numeric)
    pco_vals = ifelse(pco_vals < 10 | pco_vals > 120, NA_real_, pco_vals)
    DT[, `:=`(paco2 = pco_vals, sirs_paco2_i = as.integer(!is.na(pco_vals) & pco_vals < 32))]
  } else DT[, `:=`(paco2 = NA_real_, sirs_paco2_i = NA_integer_)]
}
rm(labs, wbc_tbl, bands_tbl, pco_tbl, j, age, wbc_vals, bands_vals, pco_vals); gc()

# ====================== Totals helper & final scores ==========================

rsum_or_na = function(DT, cols, out) {
  a = as.matrix(DT[, ..cols]); all_na = rowSums(!is.na(a)) == 0
  DT[, (out) := as.integer(ifelse(all_na, NA, rowSums(a, na.rm = TRUE)))]
}

rsum_or_na(DT, c("qsofa_rr_i","qsofa_sbp_i","qsofa_ment_i"), "qsofa")
rsum_or_na(DT, c("mews_rr_i","mews_hr_i","mews_sbp_i","mews_temp_i","mews_avpu_i"), "mews")
rsum_or_na(DT, c("news_rr_i","news_spo2_i","news_temp_i","news_sbp_i","news_hr_i","news_conf_i","news_o2_i"), "news")
rsum_or_na(DT, c("sirs_temp_i","sirs_hr_i","sirs_rr_i","sirs_wbc_i","sirs_bands_i","sirs_paco2_i"), "sirs")
rsum_or_na(DT, c("mews_rr_i","mews_hr_i","mews_sbp_i","mews_temp_i","mews_avpu_i","mews_sf_i"), "mews_sf_total")

# ====================== Write Parquet =========================================

# Keep integer columns as integers (Arrow will store them efficiently)
write_parquet(DT, out_parquet)
cat("Wrote:", out_parquet, "rows:", nrow(DT), "encounters:", length(unique(DT$hospitalization_id)), "\n")
DT[]
