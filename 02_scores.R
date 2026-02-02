
# create risk scores by looping over vital sign and lab components
# process one vital/lab at a time, assign all relevant scores, then combine

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

# missingness characterization (encounter-level) -------------------------------

message("\n== Characterizing encounter-level missingness for vitals/labs ==")

## scope: ED admits only
scope_encs = fsubset(cohort, ed_admit_01 == 1)$joined_hosp_id

## check which encounters have ≥1 measurement of each type

# vitals - just need joined_hosp_id presence in each list element
has_hr   = funique(vitals_list$heart_rate$joined_hosp_id)
has_rr   = funique(vitals_list$respiratory_rate$joined_hosp_id)
has_temp = funique(vitals_list$temp_c$joined_hosp_id)
has_spo2 = funique(vitals_list$spo2$joined_hosp_id)
has_gcs  = funique(vitals_list$gcs$joined_hosp_id)

# wbc - check labs before the pivot
has_wbc = 
  fsubset(labs, !is.na(sirs_wbc)) |>
  fselect(joined_hosp_id) |>
  funique() |>
  tibble::deframe()

## summarize
n_total = length(scope_encs)

miss_vitals_labs = tidytable(
  variable  = c("heart_rate", "resp_rate", "temp", "spo2", "gcs", "wbc"),
  n_total   = n_total,
  n_missing = c(
    sum(!scope_encs %in% has_hr),
    sum(!scope_encs %in% has_rr),
    sum(!scope_encs %in% has_temp),
    sum(!scope_encs %in% has_spo2),
    sum(!scope_encs %in% has_gcs),
    sum(!scope_encs %in% has_wbc)
  )
) |>
  ftransform(
    pct_missing = round(100 * n_missing / n_total, 2),
    site        = site_lowercase
  )

fwrite(miss_vitals_labs,here("upload_to_box", paste0("missing_vlab_", site_lowercase, ".csv")))

rm(has_hr, has_rr, has_temp, has_spo2, has_gcs, has_wbc, miss_vitals_labs)

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

## run locf: vitals 4h, labs 12h -----------------------------------------------

for (cn in vital_cols) scores = locf(scores, cn, hours = 4)
for (cn in lab_cols)   scores = locf(scores, cn, hours = 12)

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

# scores on the wards only -----------------------------------------------------

scores = 
  tidytable(scores) |>
  ftransform(mews_sf_total = mews_total + sf) |>
  join(ward_times, how = "inner", multiple = T) |>
  fsubset(time >= in_dttm & time <= out_dttm) |>
  select(-hospitalization_id) |>
  fgroup_by(joined_hosp_id, in_dttm, out_dttm, time) |>
  fmax() |>
  fgroup_by(joined_hosp_id) |>
  fmutate(
    in_dttm  = fmin(in_dttm),
    out_dttm = fmax(out_dttm)
  ) |>
  fungroup() |>
  funique()

# add primary outcome ----------------------------------------------------------

outcomes = fread(here("proj_tables/outcome_times.csv"))

scores = 
  join(scores, outcomes, how = "left", multiple = T) |>
  fmutate(o_primary_01 = if_else(!is.na(outcome_dttm),         1L, 0L)) |>
  fmutate(o_nohospc_01 = if_else(!is.na(outcome_nohospc_dttm), 1L, 0L)) |>
  fmutate(end_dttm     = pmin(out_dttm, outcome_dttm, na.rm = T)) |>
  fsubset(time < end_dttm) 

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

# save scores ------------------------------------------------------------------

write_parquet(scores, here("proj_tables", "scores_full.parquet"))

rm(list = ls()); gc() # clear memory to help 03 run

# go to 03

################################################################################
