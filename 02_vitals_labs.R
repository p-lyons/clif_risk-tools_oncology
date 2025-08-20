# labs + vitals EWS scoring (chunked, memory-efficient)
# Requires 01_* complete (expects cohort + hid_jid_crosswalk in memory or on disk)

# ------------------------------------------------------------------------------
# parameters
# ------------------------------------------------------------------------------

union_round_min   = 15L       # grid rounding (minutes); set NULL to skip rounding
vital_window_min  = 240L      # max age (min) to carry forward vitals
resp_window_min   = 60L       # max age (min) to carry forward respiratory support
lab_window_wbc    = 720L      # max age (min) for WBC
lab_window_bands  = 720L      # max age (min) for bands
lab_window_paco2  = 240L      # max age (min) for PaCO2
tz_use            = if (exists("site_time_zone")) site_time_zone else "UTC"
fio2_units        = "auto"    # "auto" | "fraction" | "percent"
out_parquet       = "scores_ews.parquet"

# chunking
chunk_size        = 50L
final_dir         = here("proj_tables", tools::file_path_sans_ext(out_parquet)) # dataset dir
temp_dir          = here("proj_temp", "ews_chunks") # temporary storage per chunk

# ------------------------------------------------------------------------------
# inputs from 00, 01 (reuse in-memory; otherwise load)
# ------------------------------------------------------------------------------

if (!exists("req_vitals")) {
  req_vitals = c("heart_rate","respiratory_rate","sbp","dbp","temp_c","spo2")
}

if (!exists("hid_jid_crosswalk")) {
  cw_path = here("proj_tables", "hid_jid_crosswalk.parquet")
  if (file.exists(cw_path)) hid_jid_crosswalk = arrow::read_parquet(cw_path)
}
if (!exists("hid_jid_crosswalk"))
  stop("hid_jid_crosswalk is required (hospitalization_id -> joined_hosp_id). Run 01_* or load crosswalk.", call. = FALSE)

if (!exists("cohort")) {
  coh_path = here("proj_tables", "cohort.parquet")
  if (file.exists(coh_path)) cohort = arrow::read_parquet(coh_path)
}
if (!exists("cohort"))
  stop("cohort is required. Run 01_*.", call. = FALSE)

setDT(hid_jid_crosswalk)
hid_jid_crosswalk = unique(hid_jid_crosswalk[, .(hospitalization_id, joined_hosp_id, admission_dttm, discharge_dttm)])
hid_jid_crosswalk[, `:=`(
  admission_dttm = as.POSIXct(admission_dttm, tz = tz_use),
  discharge_dttm = as.POSIXct(discharge_dttm, tz = tz_use)
)]

# ------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------

safe_rowsum_cols <- function(dt, score_cols) {
  cols_exist = intersect(score_cols, names(dt))
  if (!length(cols_exist)) return(rep(NA_integer_, nrow(dt)))
  m = as.matrix(dt[, ..cols_exist])
  has_any = rowSums(!is.na(m)) > 0
  out = integer(nrow(dt))
  out[has_any] = as.integer(rowSums(m[has_any, , drop = FALSE], na.rm = TRUE))
  out[!has_any] = NA_integer_
  out
}

norm_fio2_fraction <- function(x, mode = "auto") {
  if (is.null(x)) return(x)
  if (mode == "fraction") return(pmin(pmax(x, 0.21), 1.00))
  if (mode == "percent")  return(pmin(pmax(x/100, 0.21), 1.00))
  # auto
  med = suppressWarnings(median(x, na.rm = TRUE))
  frac = if (is.finite(med) && med <= 1) x else x/100
  pmin(pmax(frac, 0.21), 1.00)
}

# ------------------------------------------------------------------------------
# per-patient processor (works on pre-filtered chunk tables)
# ------------------------------------------------------------------------------

process_patient_ews <- function(jid, V, L, R, W) {
  # W: ward intervals for this jid  (joined_hosp_id, in_dttm, out_dttm)
  if (nrow(W) == 0L) return(NULL)
  
  # build time grid from all signals + ward boundaries
  tp = list()
  k  = 1L
  if (nrow(V)) { tp[[k]] = V$time; k = k + 1L }
  if (nrow(L)) { tp[[k]] = L$time; k = k + 1L }
  if (nrow(R)) { tp[[k]] = R$time; k = k + 1L }
  tp[[k]] = c(W$in_dttm, W$out_dttm)
  
  times = unique(do.call(c, tp))
  if (!length(times)) return(NULL)
  
  if (!is.null(union_round_min)) {
    step = union_round_min * 60
    times = unique(as.POSIXct((as.numeric(times) %/% step) * step, origin = "1970-01-01", tz = tz_use))
  }
  
  G = data.table(joined_hosp_id = jid, recorded_dttm = sort(times))
  setkey(G, recorded_dttm)
  
  # ------------------ vitals (LOCF + expiry) ----------------------------------
  if (nrow(V)) {
    # roll one vital at a time to minimize RAM
    vit_types = intersect(unique(V$vital_category), c("heart_rate","respiratory_rate","sbp","dbp","temp_c","spo2"))
    for (vt in vit_types) {
      X = V[vital_category == vt, .(time, value = vital_value)]
      if (!nrow(X)) { G[, (vt) := NA_real_ ]; next }
      setkey(X, time)
      J = X[G, on = .(time = recorded_dttm), roll = TRUE]
      age = as.numeric(difftime(G$recorded_dttm, J$time, units = "mins"))
      age[age < 0] = NA_real_
      G[, (vt) := fifelse(is.na(age) | age > vital_window_min, NA_real_, J$value)]
    }
  } else {
    for (nm in c("heart_rate","respiratory_rate","sbp","dbp","temp_c","spo2")) G[, (nm) := NA_real_]
  }
  
  # rename to compact names used in scoring
  setnames(G,
           old = c("heart_rate","respiratory_rate","temp_c"),
           new = c("hr","rr","temp_c"),
           skip_absent = TRUE)
  
  # ------------------ labs (ward-only, LOCF + expiry) -------------------------
  if (nrow(L)) {
    # filter to ward intervals via foverlaps
    L2 = copy(L)[, `:=`(start = time, end = time)]
    W2 = W[, .(start = in_dttm, end = out_dttm)]
    setkey(L2, start, end); setkey(W2, start, end)
    L2 = foverlaps(L2, W2, nomatch = 0L)[, .(time = start, lab_category, lab_value_numeric)]
    
    if (nrow(L2)) {
      # roll per lab with distinct windows
      for (lab in c("wbc","bands","paco2")) {
        X = L2[lab_category == lab, .(time, value = lab_value_numeric)]
        if (!nrow(X)) { G[, (lab) := NA_real_ ]; next }
        setkey(X, time)
        J = X[G, on = .(time = recorded_dttm), roll = TRUE]
        age = as.numeric(difftime(G$recorded_dttm, J$time, units = "mins"))
        age[age < 0] = NA_real_
        
        win = switch(lab,
                     wbc   = lab_window_wbc,
                     bands = lab_window_bands,
                     paco2 = lab_window_paco2)
        G[, (lab) := fifelse(is.na(age) | age > win, NA_real_, J$value)]
      }
    } else {
      for (lab in c("wbc","bands","paco2")) G[, (lab) := NA_real_]
    }
  } else {
    for (lab in c("wbc","bands","paco2")) G[, (lab) := NA_real_]
  }
  
  # ------------------ respiratory support (ward-only) -------------------------
  if (nrow(R)) {
    R2 = copy(R)[, `:=`(start = time, end = time)]
    W2 = W[, .(start = in_dttm, end = out_dttm)]
    setkey(R2, start, end); setkey(W2, start, end)
    R2 = foverlaps(R2, W2, nomatch = 0L)[, .(time = start, lpm_set, fio2_set)]
    
    if (nrow(R2)) {
      setkey(R2, time)
      J = R2[G, on = .(time = recorded_dttm), roll = TRUE]
      age = as.numeric(difftime(G$recorded_dttm, J$time, units = "mins"))
      age[age < 0] = NA_real_
      lpm  = fifelse(is.na(age) | age > resp_window_min, NA_real_, J$lpm_set)
      fio2 = fifelse(is.na(age) | age > resp_window_min, NA_real_, J$fio2_set)
      
      fio2_frac = norm_fio2_fraction(fio2, mode = fio2_units)
      on_o2     = as.integer((!is.na(lpm) & lpm > 0) | (!is.na(fio2_frac) & fio2_frac > 0.21 + 1e-9))
      news_o2_i = fifelse(is.na(on_o2), NA_integer_, fifelse(on_o2 == 1L, 2L, 0L))
      sf_ratio  = fifelse(!is.na(G$spo2) & !is.na(fio2_frac), G$spo2 / fio2_frac, NA_real_)
      
      G[, `:=`(lpm_set = lpm, fio2_raw = fio2, fio2_frac = fio2_frac,
               on_o2 = on_o2, news_o2_i = news_o2_i, sf_ratio = sf_ratio)]
    } else {
      G[, `:=`(lpm_set = NA_real_, fio2_raw = NA_real_, fio2_frac = NA_real_,
               on_o2 = NA_integer_, news_o2_i = NA_integer_, sf_ratio = NA_real_)]
    }
  } else {
    G[, `:=`(lpm_set = NA_real_, fio2_raw = NA_real_, fio2_frac = NA_real_,
             on_o2 = NA_integer_, news_o2_i = NA_integer_, sf_ratio = NA_real_)]
  }
  
  # ------------------ component scores (vital-only + O2) ----------------------
  
  # qSOFA
  G[, `:=`(
    qsofa_rr_i   = fcase(!is.na(rr)  & rr  >= 22, 1L, !is.na(rr),  0L, default = NA_integer_),
    qsofa_sbp_i  = fcase(!is.na(sbp) & sbp <= 100, 1L, !is.na(sbp), 0L, default = NA_integer_),
    qsofa_ment_i = NA_integer_
  )]
  
  # SIRS
  G[, `:=`(
    sirs_temp_i = fcase(!is.na(temp_c) & (temp_c > 38 | temp_c < 36), 1L,
                        !is.na(temp_c), 0L, default = NA_integer_),
    sirs_hr_i   = fcase(!is.na(hr)     & hr > 90, 1L, !is.na(hr), 0L, default = NA_integer_),
    sirs_rr_i   = fcase(!is.na(rr)     & rr > 20, 1L, !is.na(rr), 0L, default = NA_integer_)
  )]
  
  # MEWS
  G[, `:=`(
    mews_rr_i = fcase(is.na(rr), NA_integer_,
                      rr >= 30, 3L, rr >= 21, 2L, rr >= 15, 1L, rr >= 9, 0L, rr <= 8, 2L),
    mews_hr_i = fcase(is.na(hr), NA_integer_,
                      hr >= 130, 3L, hr >= 111, 2L, hr >= 101, 1L, hr >= 51, 0L, hr >= 41, 1L, hr <= 40, 2L),
    mews_sbp_i = fcase(is.na(sbp), NA_integer_,
                       sbp >= 200, 2L, sbp >= 101, 0L, sbp >= 81, 1L, sbp >= 71, 2L, sbp <= 70, 3L),
    mews_temp_i = fcase(is.na(temp_c), NA_integer_,
                        temp_c >= 38.5, 2L, temp_c > 35.0, 0L, temp_c <= 35.0, 2L),
    mews_avpu_i = NA_integer_
  )]
  
  # NEWS (Scale 1 for SpO2; O2 added above)
  G[, `:=`(
    news_rr_i = fcase(is.na(rr), NA_integer_,
                      rr >= 25, 3L, rr >= 21, 2L, rr >= 12, 0L, rr >= 9, 1L, rr <= 8, 3L),
    news_spo2_i = fcase(is.na(spo2), NA_integer_,
                        spo2 >= 96, 0L, spo2 >= 94, 1L, spo2 >= 92, 2L, spo2 <= 91, 3L),
    news_hr_i = fcase(is.na(hr), NA_integer_,
                      hr >= 131, 3L, hr >= 111, 2L, hr >= 91, 1L, hr >= 51, 0L, hr >= 41, 1L, hr <= 40, 3L),
    news_sbp_i = fcase(is.na(sbp), NA_integer_,
                       sbp >= 220, 3L, sbp >= 111, 0L, sbp >= 101, 1L, sbp >= 91, 2L, sbp <= 90, 3L),
    news_temp_i = fcase(is.na(temp_c), NA_integer_,
                        temp_c >= 39.1, 2L, temp_c >= 38.1, 1L, temp_c >= 36.1, 0L, temp_c >= 35.1, 1L, temp_c <= 35.0, 3L),
    news_conf_i = NA_integer_
  )]
  
  # MEWS-SF (optional)
  if ("sf_ratio" %in% names(G)) {
    G[, mews_sf_i := fcase(is.na(sf_ratio), NA_integer_,
                           sf_ratio > 315, 0L, sf_ratio > 235, 2L, sf_ratio <= 235, 3L)]
  } else G[, mews_sf_i := NA_integer_]
  
  # ------------------ SIRS lab components -------------------------------------
  
  if ("wbc" %in% names(G))   G[, sirs_wbc_i   := fcase(!is.na(wbc)   & (wbc < 4 | wbc > 12), 1L, !is.na(wbc), 0L, default = NA_integer_)]
  if ("bands" %in% names(G)) G[, sirs_bands_i := fcase(!is.na(bands) & bands > 10, 1L, !is.na(bands), 0L, default = NA_integer_)]
  if ("paco2" %in% names(G)) G[, sirs_paco2_i := fcase(!is.na(paco2) & paco2 < 32, 1L, !is.na(paco2), 0L, default = NA_integer_)]
  
  # ------------------ totals ---------------------------------------------------
  
  G[, `:=`(
    qsofa_total = safe_rowsum_cols(.SD, c("qsofa_rr_i","qsofa_sbp_i","qsofa_ment_i")),
    sirs_total  = safe_rowsum_cols(.SD, c("sirs_temp_i","sirs_hr_i","sirs_rr_i","sirs_wbc_i","sirs_bands_i","sirs_paco2_i")),
    mews_total  = safe_rowsum_cols(.SD, c("mews_rr_i","mews_hr_i","mews_sbp_i","mews_temp_i","mews_avpu_i")),
    mews_sf_tot = safe_rowsum_cols(.SD, c("mews_rr_i","mews_hr_i","mews_sbp_i","mews_temp_i","mews_avpu_i","mews_sf_i")),
    news_total  = safe_rowsum_cols(.SD, c("news_rr_i","news_spo2_i","news_temp_i","news_sbp_i","news_hr_i","news_conf_i","news_o2_i"))
  )]
  
  # ------------------ mask to ward-only ---------------------------------------
  
  G2 = copy(G)[, `:=`(start = recorded_dttm, end = recorded_dttm)]
  W2 = copy(W)[, .(start = in_dttm, end = out_dttm)]
  setkey(G2, start, end); setkey(W2, start, end)
  OUT = foverlaps(G2, W2, nomatch = 0L)[, !"start|end|i.start|i.end"]
  
  OUT[]
}

# ------------------------------------------------------------------------------
# chunked pipeline
# ------------------------------------------------------------------------------

# ensure output dirs
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

all_jids = unique(hid_jid_crosswalk$joined_hosp_id)
n_pat    = length(all_jids)
n_chunks = ceiling(n_pat / chunk_size)

chunk_files = character(0)

for (chunk_i in seq_len(n_chunks)) {
  idx_start = (chunk_i - 1L) * chunk_size + 1L
  idx_end   = min(chunk_i * chunk_size, n_pat)
  chunk_jids = all_jids[idx_start:idx_end]
  
  # map to HIDs (limit all pulls to this chunk)
  chunk_hids = hid_jid_crosswalk[joined_hosp_id %in% chunk_jids, unique(hospitalization_id)]
  
  # ward times for chunk (ADT → ward-only → within admit/discharge)
  WT = data_list$adt |>
    dplyr::filter(tolower(location_category) %in% "ward",
                  hospitalization_id %in% chunk_hids) |>
    dplyr::select(hospitalization_id, in_dttm, out_dttm) |>
    dplyr::collect() |>
    setDT()
  
  if (nrow(WT)) {
    WT[, `:=`(in_dttm  = as.POSIXct(in_dttm,  tz = tz_use),
              out_dttm = as.POSIXct(out_dttm, tz = tz_use))]
    WT = merge(WT, hid_jid_crosswalk[, .(hospitalization_id, joined_hosp_id, admission_dttm, discharge_dttm)],
               by = "hospitalization_id")
    WT = WT[in_dttm >= admission_dttm & out_dttm <= discharge_dttm,
            .(joined_hosp_id, hospitalization_id, in_dttm, out_dttm)]
  } else WT = WT[, .(joined_hosp_id = integer(), hospitalization_id = integer(), in_dttm = as.POSIXct(character()), out_dttm = as.POSIXct(character()))]
  
  # vitals for chunk
  V = data_list$vitals |>
    dplyr::filter(hospitalization_id %in% chunk_hids,
                  vital_category %in% req_vitals) |>
    dplyr::select(hospitalization_id, recorded_dttm, vital_category, vital_value) |>
    dplyr::collect() |>
    setDT()
  
  if (nrow(V)) {
    V[, `:=`(time = as.POSIXct(recorded_dttm, tz = tz_use),
             vital_value = suppressWarnings(as.numeric(vital_value)))]
    V = merge(V, hid_jid_crosswalk[, .(hospitalization_id, joined_hosp_id, admission_dttm, discharge_dttm)],
              by = "hospitalization_id")
    V = V[time >= admission_dttm & time <= discharge_dttm,
          .(joined_hosp_id, time, vital_category, vital_value)]
    
    # plausibility clip (once, pre-LOCF)
    V[vital_category == "heart_rate"      & (vital_value < 20 | vital_value > 250), vital_value := NA_real_]
    V[vital_category == "respiratory_rate" & (vital_value < 5  | vital_value > 60),  vital_value := NA_real_]
    V[vital_category == "sbp"             & (vital_value < 50 | vital_value > 300), vital_value := NA_real_]
    V[vital_category == "dbp"             & (vital_value < 20 | vital_value > 200), vital_value := NA_real_]
    V[vital_category == "temp_c"          & (vital_value < 25 | vital_value > 45),  vital_value := NA_real_]
    V[vital_category == "spo2"            & (vital_value < 50 | vital_value > 100), vital_value := NA_real_]
  }
  
  # labs for chunk
  L = data_list$labs |>
    dplyr::filter(hospitalization_id %in% chunk_hids,
                  lab_category %in% c("wbc","bands","band_neutrophils","immature_neutrophils","paco2")) |>
    dplyr::select(hospitalization_id, lab_result_dttm, lab_category, lab_value_numeric) |>
    dplyr::collect() |>
    setDT()
  
  if (nrow(L)) {
    L[lab_category %in% c("bands","band_neutrophils","immature_neutrophils"), lab_category := "bands"]
    L[, `:=`(time = as.POSIXct(lab_result_dttm, tz = tz_use),
             lab_value_numeric = suppressWarnings(as.numeric(lab_value_numeric)))]
    L = merge(L, hid_jid_crosswalk[, .(hospitalization_id, joined_hosp_id, admission_dttm, discharge_dttm)],
              by = "hospitalization_id")
    L = L[time >= admission_dttm & time <= discharge_dttm,
          .(joined_hosp_id, time, lab_category, lab_value_numeric)]
    
    # plausibility clip
    L[lab_category == "wbc"  & (lab_value_numeric < 0.1 | lab_value_numeric > 200), lab_value_numeric := NA_real_]
    L[lab_category == "bands"& (lab_value_numeric < 0   | lab_value_numeric > 100), lab_value_numeric := NA_real_]
    L[lab_category == "paco2"& (lab_value_numeric < 10  | lab_value_numeric > 120), lab_value_numeric := NA_real_]
  }
  
  # respiratory support for chunk
  R = data_list$respiratory_support |>
    dplyr::filter(hospitalization_id %in% chunk_hids) |>
    dplyr::select(hospitalization_id, recorded_dttm, lpm_set, fio2_set) |>
    dplyr::collect() |>
    setDT()
  
  if (nrow(R)) {
    R[, `:=`(time = as.POSIXct(recorded_dttm, tz = tz_use),
             lpm_set = suppressWarnings(as.numeric(lpm_set)),
             fio2_set = suppressWarnings(as.numeric(fio2_set)))]
    R = merge(R, hid_jid_crosswalk[, .(hospitalization_id, joined_hosp_id, admission_dttm, discharge_dttm)],
              by = "hospitalization_id")
    R = R[time >= admission_dttm & time <= discharge_dttm,
          .(joined_hosp_id, time, lpm_set, fio2_set)]
  }
  
  # process each jid in this chunk
  chunk_out = vector("list", length(chunk_jids))
  for (i in seq_along(chunk_jids)) {
    jid = chunk_jids[i]
    Wj  = WT[joined_hosp_id == jid, .(joined_hosp_id, in_dttm, out_dttm)]
    Vj  = if (nrow(V)) V[joined_hosp_id == jid] else V
    Lj  = if (nrow(L)) L[joined_hosp_id == jid] else L
    Rj  = if (nrow(R)) R[joined_hosp_id == jid] else R
    res = process_patient_ews(jid, Vj, Lj, Rj, Wj)
    if (!is.null(res) && nrow(res)) chunk_out[[i]] = res
  }
  
  # write this chunk (multi-file dataset; one file per chunk)
  ck = rbindlist(chunk_out, use.names = TRUE, fill = TRUE)
  if (nrow(ck)) {
    fp = file.path(temp_dir, sprintf("ews_chunk_%03d.parquet", chunk_i))
    arrow::write_parquet(ck, fp)
    chunk_files = c(chunk_files, fp)
  }
  
  rm(WT, V, L, R, chunk_out, ck); gc()
}

# move chunk files to final dataset dir (leave as multi-file parquet dataset)
if (length(chunk_files)) {
  file.copy(from = chunk_files, to = final_dir, overwrite = TRUE, recursive = FALSE)
}

# optional: write a single-file parquet (not recommended for huge cohorts)
# ds <- arrow::open_dataset(final_dir, format = "parquet")
# arrow::write_parquet(dplyr::collect(ds), here("proj_tables", out_parquet))

# cleanup temp
unlink(temp_dir, recursive = TRUE)
