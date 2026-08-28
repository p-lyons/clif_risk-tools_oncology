# Cohort script for CLIF project validating risk tools in oncology.
# Requires data_list to be loaded/validated from 00_*

# check/create output directories ---------------------------------------------

if (!exists("BOX_DIR")) {
  stop("BOX_DIR not found. Did you run 00_setup first?", call. = FALSE)
}

if (!dir.exists(here(BOX_DIR))) {
  dir.create(here(BOX_DIR), recursive = TRUE)
}

if (!dir.exists(here("proj_tables"))) {
  dir.create(here("proj_tables"), recursive = TRUE)
}

# resources for RAM heavy wrangling --------------------------------------------

## cores & RAM (reuse from 00 if available) ------------------------------------

os_type   = if (exists("os_type"))   os_type   else Sys.info()[["sysname"]]
all_cores = if (exists("all_cores")) all_cores else {
  x = parallel::detectCores(logical = TRUE); if (is.na(x)) 1L else as.integer(x)
}

if (!exists("avail_ram_gb") || !is.finite(avail_ram_gb)) {
  get_ram_gb = function() {
    tryCatch({
      if (Sys.info()[["sysname"]] == "Darwin") {
        bytes = suppressWarnings(as.numeric(system("sysctl -n hw.memsize", intern = TRUE)))
        if (length(bytes) > 0 && !is.na(bytes)) bytes / 1024^3 else NA_real_
      } else if (file.exists("/proc/meminfo")) {
        kb = suppressWarnings(as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE)))
        if (length(kb) > 0 && !is.na(kb)) kb / 1024^2 else NA_real_
      } else if (requireNamespace("ps", quietly = TRUE)) {
        ps::ps_system_memory()[["available"]] / 1024^3
      } else NA_real_
    }, error = function(e) NA_real_)
  }
  avail_ram_gb = get_ram_gb()
}

## choose threads (2 GB per thread; leave 1 core free) -------------------------

reserve_cores  = 1L
gb_per_thread  = 2.0
max_by_cores   = max(1L, all_cores - reserve_cores)
max_by_memory  = if (is.finite(avail_ram_gb)) max(1L, floor(avail_ram_gb / gb_per_thread)) else max_by_cores
n_threads      = as.integer(max(1L, min(max_by_cores, max_by_memory)))
n_math_threads = as.integer(max(1L, min(n_threads, 8L)))

## apply thread settings -------------------------------------------------------

data.table::setDTthreads(threads = n_threads)
collapse::set_collapse(nthreads  = n_threads)
options(arrow.use_threads        = TRUE)
Sys.setenv(ARROW_NUM_THREADS     = n_threads)
options(mc.cores                 = n_threads)

message(
  sprintf("01 resources | OS=%s | Cores=%d | Threads=%d | MathThreads=%d | Avail RAM≈%s GB (rule: 2 GB/core)",
          os_type, all_cores, n_threads, n_math_threads,
          ifelse(is.finite(avail_ram_gb), round(avail_ram_gb, 1), "NA"))
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# reclassify stepdown as ward throughout ADT -----------------------------------
# stepdown beds are functionally wards for deterioration monitoring purposes

data_list$adt = data_list$adt |>
  dplyr::mutate(
    location_category = dplyr::case_when(
      location_category %in% c("stepdown", "Stepdown") ~ "ward",
      TRUE ~ location_category
    )
  )

message("✅ Reclassified stepdown → ward in ADT")

# cohort identification --------------------------------------------------------

## start by linking contiguous hospitalizations --------------------------------

### encounters with age >= 18 and dates within study window --------------------

### exclude hospitalizations carrying more than one patient_id ----------------
# CLIF defines hospitalization_id as a key of one, so a hospitalization_id
# appearing under two patient_id values is an electronic-health-record artifact
# rather than a data-model feature. The known mechanism is an unidentified patient
# admitted under a placeholder chart whose real chart is created alongside it; the
# two merge on the source system's next refresh, but a snapshot taken before that
# refresh carries both.
#
# Left in place the defect does not surface as duplicate rows, which is why no
# existing guard catches it. The linkage below groups within patient_id, so one
# hospitalization under two patients yields two distinct joined_hosp_id values.
# hid_jid_crosswalk then maps that one hospitalization to both, and every
# measurement join fans its vitals, labs, and ADT out to both encounters. The
# group-by collapses that follow all run within joined_hosp_id, so the inflation
# appears as extra encounters rather than as extra rows.
#
# These hospitalizations are dropped rather than resolved to one patient_id. CLIF
# carries no field that distinguishes the placeholder chart from the real one, and
# choosing wrong attaches the wrong demographics and the wrong death_dttm to a
# real encounter, which can convert a ward death into a censored discharge. The
# count is carried into the run report so the coordinating center can see which
# sites are affected and at what scale. A site with none sees no change.

hosp_patient_map =
  dplyr::select(data_list$hospitalization, patient_id, hospitalization_id) |>
  dplyr::collect() |>
  distinct()

hid_multi_patient =
  fcount(hosp_patient_map, hospitalization_id) |>
  fsubset(N > 1) |>
  pull(hospitalization_id)

run_log$n_hid_multi_patient = length(hid_multi_patient)

if (length(hid_multi_patient) > 0) {
  n_multi_rows = nrow(fsubset(hosp_patient_map, hospitalization_id %in% hid_multi_patient))
  message(
    sprintf("⚠️  Excluding %s hospitalization_id(s) carrying more than one patient_id (%s chart rows). IDs withheld from the log; inspect hid_multi_patient interactively.",
            format(length(hid_multi_patient), big.mark = ","),
            format(n_multi_rows,              big.mark = ","))
  )
  rm(n_multi_rows)
} else {
  message("✅ No hospitalization_id carries more than one patient_id.")
}

rm(hosp_patient_map)

hosp_blocks =
  dplyr::filter(data_list$hospitalization, age_at_admission >= 18) |>
  dplyr::filter(admission_dttm >= start_date & admission_dttm <= end_date) |>
  dplyr::filter(admission_dttm < discharge_dttm & !is.na(discharge_dttm)) |>
  dplyr::select(patient_id, hospitalization_id, admission_dttm, discharge_dttm) |>
  dplyr::collect() |>
  distinct() |>
  fsubset(!hospitalization_id %in% hid_multi_patient) |>
  roworder(patient_id, admission_dttm)

### use data.table to find joined hospitalizations with <= 6h gaps -------------

link_hours = 6L
linked     = as.data.table(hosp_blocks)
setorder(linked, patient_id, admission_dttm)

#### calculate gaps between encounters and mark encounters that should be linked
linked[, next_admit := shift(admission_dttm, type = "lead"),           by    = patient_id]
linked[, next_gap   := as.numeric(difftime(next_admit, discharge_dttm, units = "hours"))]
linked[, prev_dc    := shift(discharge_dttm, type = "lag"),            by    = patient_id]
linked[, prev_gap   := as.numeric(difftime(admission_dttm, prev_dc,    units = "hours"))]
linked[, link_flag  := (next_gap < link_hours | prev_gap < link_hours)]
linked[is.na(link_flag), link_flag := FALSE]

#### create unique joined hospitalization ID, new group whenever gap > link_hours from  previous discharge
linked[, new_group := is.na(prev_gap) | prev_gap >= link_hours]
linked[, joined_hosp_id := .GRP, by = .(patient_id, cumsum(new_group))]

#### create hid_jid_crosswalk --------------------------------------------------
hid_jid_crosswalk = select(linked, ends_with("id")) |> as_tidytable()

## hospital ward admissions ----------------------------------------------------

### inpatient stay requires wards ----------------------------------------------

inpatient_hids =
  dplyr::filter(data_list$adt, tolower(location_category) %in% c("ward")) |>
  dplyr::select(hospitalization_id) |>
  dplyr::collect()

inpatient_hids =
  funique(inpatient_hids) |>
  tibble::deframe()

inpatient_jids =
  fsubset(hid_jid_crosswalk, hospitalization_id %in% inpatient_hids) |>
  select(joined_hosp_id) |>
  funique()

inpatient_jids =
  funique(inpatient_jids) |>
  tibble::deframe()

### don't want to include obstetrics/psych -------------------------------------

drop_ob =
  dplyr::filter(data_list$adt, tolower(location_category) %in% c("l&d", "psych", "rehab")) |>
  dplyr::select(hospitalization_id) |>
  dplyr::collect()

drop_ob =
  funique(drop_ob) |>
  tibble::deframe()

drop_ob_jids =
  fsubset(hid_jid_crosswalk, hospitalization_id %in% drop_ob) |>
  fselect(joined_hosp_id) |>
  funique() |>
  tibble::deframe()

linked = fsubset(linked,  joined_hosp_id %in% inpatient_jids)
linked = fsubset(linked, !joined_hosp_id %in% drop_ob_jids)
linked = select(linked, ends_with("id"), ends_with("dttm"))

### an admission requires at least 1 full set of vital signs on the wards ------

ward_times =
  dplyr::filter(data_list$adt, tolower(location_category) %in% c("ward")) |>
  dplyr::filter(hospitalization_id %in% inpatient_hids) |>
  dplyr::select(hospitalization_id, in_dttm, out_dttm) |>
  dplyr::collect()

ward_times =
  join(ward_times, linked, how = "inner", multiple = T) |>
  fsubset(in_dttm >= admission_dttm & out_dttm <= discharge_dttm) |>
  fselect(joined_hosp_id, hospitalization_id, in_dttm, out_dttm) |>
  distinct()

has_vital_signs =
  dplyr::filter(data_list$vitals, hospitalization_id %in% inpatient_hids) |>
  dplyr::filter(vital_category %in% req_vitals) |>
  dplyr::select(hospitalization_id, vital_category, recorded_dttm) |>
  dplyr::collect() |>
  distinct()

has_vital_signs =
  join(has_vital_signs, ward_times, how = "inner", multiple = T) |>
  fsubset(recorded_dttm >= in_dttm & recorded_dttm <= out_dttm) |>
  fselect(joined_hosp_id, vital_category) |>
  distinct()

has_vital_signs =
  fgroup_by(has_vital_signs, joined_hosp_id) |>
  fnobs() |>
  fsubset(vital_category== length(req_vitals)) |>
  pull(joined_hosp_id)

### double-check study years ---------------------------------------------------

### data for cohort development ------------------------------------------------

linked            = fsubset(linked, joined_hosp_id %in% has_vital_signs)
hid_jid_crosswalk = select(linked, ends_with("id"))
cohort_hids       = funique(hid_jid_crosswalk$hospitalization_id)
cohort_pats       = funique(hid_jid_crosswalk$patient_id)

### clean up helpers -----------------------------------------------------------

rm(inpatient_hids, inpatient_jids, drop_ob, drop_ob_jids, has_vital_signs, hosp_blocks, link_hours); gc()

## assemble cohort data frame --------------------------------------------------

#### pull additional data for cohort filtering and final variables
cohort_data =
  dplyr::filter(data_list$hospitalization, hospitalization_id %in% cohort_hids) |>
  dplyr::select(ends_with("id"), age_at_admission, discharge_category) |>
  dplyr::collect() |>
  distinct()

# The distinct() above absorbs exactly repeated chart rows, and the multi-patient
# exclusion above removes the placeholder-chart case, so anything still reaching
# this guard is a hospitalization whose attribute rows genuinely conflict. That is
# a defect worth halting on rather than resolving silently.

hid_dups_source =
  fcount(cohort_data, hospitalization_id) |>
  fsubset(N > 1)

if (nrow(hid_dups_source) > 0) {
  stop(
    sprintf("Source has duplicate hospitalization_id: %s",
            paste(head(hid_dups_source$hospitalization_id, 50), collapse = ", ")),
    call. = FALSE
  )
}

#### create final cohort - 1 row per joined_hosp_id
cohort =
  join(linked, cohort_data, how = "left", multiple = T) |>
  roworder(admission_dttm) |>
  fgroup_by(patient_id, joined_hosp_id) |>
  fsummarize(
    age                = ffirst(age_at_admission),
    admission_dttm     = ffirst(admission_dttm),
    discharge_dttm     = flast(discharge_dttm),
    discharge_category = flast(discharge_category)
  )

#### clean up temporary variables
rm(linked, cohort_data); gc()

## quality control -------------------------------------------------------------

### check for duplicates -------------------------------------------------------

dupes = cohort |> janitor::get_dupes(patient_id, admission_dttm)

if (nrow(dupes) > 0) {
  dup_ids = funique(dupes$joined_hosp_id)
  stop(
    sprintf("Found %d duplicate joined_hosp_id(s): %s", length(dup_ids), paste(dup_ids, collapse = ", ")),
    call. = FALSE
  )
}

message("✅ No duplicate joined_hosp_id found.")

### YODO (you only die once) ---------------------------------------------------

#### death duplicates ----------------------------------------------------------

dup_deaths =
  fsubset(cohort, tolower(discharge_category) == "expired") |>
  roworder(admission_dttm, discharge_dttm) |>
  fgroup_by(patient_id) |>
  fmutate(one = 1L) |>
  fmutate(n_deaths = fsum(one), counter = fcumsum(one)) |>
  fungroup() |>
  fsubset(n_deaths > 1)

if (nrow(dup_deaths) > 0) {

  encdrop =
    fsubset(dup_deaths, counter > 1) |>
    pull(joined_hosp_id)

  cohort  = fsubset(cohort, !joined_hosp_id %in% encdrop)
  n_dupes = length(encdrop)
  n_pats  = length(unique(dup_deaths$patient_id))

  message(
    sprintf("Removed %d duplicate deaths for %d patients. IDs withheld from the log; inspect dup_deaths interactively.",
            n_dupes, n_pats)
  )
}

#### readmissions following death ----------------------------------------------

death_times =
  fsubset(cohort, tolower(discharge_category) == "expired") |>
  roworder(discharge_dttm) |>
  fgroup_by(patient_id) |>
  fsummarise(death_instant = ffirst(discharge_dttm))

post_death_admissions =
  join(cohort, death_times, how = "inner", multiple = T) |>
  fsubset(admission_dttm >= death_instant)

if (nrow(post_death_admissions) > 0) {
  cohort = fsubset(cohort, !joined_hosp_id %in% post_death_admissions$joined_hosp_id)

  message(
    sprintf("Removed %d post-death encs for %d patients (? organ donors). IDs withheld from the log; inspect post_death_admissions interactively.",
            nrow(post_death_admissions),
            length(unique(post_death_admissions$patient_id)))
  )
}

message("✅ Cleaned duplicate deaths and post-death encounters.")

suppressWarnings(rm(dupes, dup_deaths, death_times, post_death_admissions))
suppressWarnings(rm(start_date, end_date, encdrop, file_type, n_dupes, n_pats)); gc()

cohort_pats = funique(cohort$patient_id)
cohort_jids = funique(cohort$joined_hosp_id)
cohort_hids = funique(hid_jid_crosswalk$hospitalization_id)
date_frame  = select(cohort, patient_id, joined_hosp_id, ends_with("dttm"))

## identify cancer patients ----------------------------------------------------

### load cancer dx codes -------------------------------------------------------

additional_drops = c(
  "D72.11", "D72.110", "D72.111", "D72.118", "D72.119",
  "D47.9", "D47.Z9", "D45", "D47.3", "D47.02"
)

ca_codes =
  readxl::read_xlsx(here("config/icd10cm_casefinding_2023.xlsx")) |> #### exclude C44, in remission
  janitor::clean_names() |>
  fsubset(general_category %in% c("Hematopoietic neoplasm", "Malignant neoplasm")) |>
  fsubset(is.na(drop) | drop != "drop") |>
  fsubset(!str_detect(icd_10_cm_code_specific, "^Z85")) |>
  fsubset(!str_detect(icd_10_cm_code_definition, "in remission|personal history")) |>
  fsubset(!icd_10_cm_code_specific %in% additional_drops) |>
  fmutate(liquid_01 = if_else(general_category == "Hematopoietic neoplasm", 1L, 0L)) |>
  fselect(diagnosis_code = icd_10_cm_code_specific, liquid_01)

ca_vect  = funique(toupper(ca_codes$diagnosis_code))

### pull all cancer codes from arrow table -------------------------------------

dx =
  dplyr::filter(data_list$hospital_diagnosis, hospitalization_id %in% cohort_hids) |>
  dplyr::filter(toupper(diagnosis_code) %in% ca_vect) |>
  dplyr::select(hospitalization_id, poa_present, diagnosis_code) |>
  dplyr::collect()

dx = distinct(dx)

### assign hierarchy to diagnoses ----------------------------------------------

diagnosis_priority = tribble(
  ~group,              ~pattern,                      ~rank,
  "Metastatic",        "^C7[7-9]|^C80",                1L,
  "Hematologic",       "^C8[1-9]|^C9[0-6]",            2L,
  "High-risk solid",   "^C22|^C25|^C34",               3L,
  "Other solid",       "^C18|^C50|^C61",               4L,
  "Other",             ".*",                           5L
)

dx =
  rowwise(dx) |>
  mutate(
    diag_group = diagnosis_priority$group[
      which.min(if_else(
        str_detect(diagnosis_code, diagnosis_priority$pattern),
        diagnosis_priority$rank,
        Inf
      ))
    ],
    rank = diagnosis_priority$rank[
      which.min(if_else(
        str_detect(diagnosis_code, diagnosis_priority$pattern),
        diagnosis_priority$rank,
        Inf
      ))
    ]
  ) |>
  ungroup()

dx =
  roworder(dx, rank, diagnosis_code) |>
  fgroup_by(hospitalization_id) |>
  fsummarize(
    diagnosis_code = ffirst(diagnosis_code),
    rank           = ffirst(rank)
  ) |>
  join(ca_codes, how = "left", multiple = F)

### cancer dx associated with encounter only -----------------------------------

dx =
  join(dx, hid_jid_crosswalk, how = "left", multiple = F) |>
  fselect(joined_hosp_id, diagnosis_code, liquid_01, rank) |>
  join(cohort, how = "inner", multiple = T)

dx_enc =
  roworder(dx, rank, diagnosis_code) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(
    ca_icd10_enc  = ffirst(diagnosis_code),
    liquid_01_enc = fmax(liquid_01),    # any heme code → liquid, regardless of hierarchy
    rank_enc      = ffirst(rank)
  )

cohort           = join(cohort, dx_enc,  how = "left", multiple = F)
cohort$ca_01     = if_else(is.na(cohort$ca_icd10_enc),  0L, 1L)
cohort$liquid_01 = if_else(is.na(cohort$liquid_01_enc), 0L, cohort$liquid_01_enc)

### tally primary cancer codes (one per encounter using priority) --------------

cancer_code_tally_primary =
  fgroup_by(dx_enc, ca_icd10_enc) |>
  fnobs() |>
  fselect(ca_icd10_enc, n = joined_hosp_id) |>
  fsubset(n > 5) |>
  roworder(-n) |>
  ftransform(site = site_lowercase)

fwrite(
  cancer_code_tally_primary,
  here(BOX_DIR, paste0("cancer_codes_primary_", site_lowercase, ".csv"))
)

### tally hospitals by type (eTable 1) ----------------------------------------
# Reviewer 2 asked for the composition of contributing hospitals. One row per
# hospital_type (academic / community / LTACH, from the ADT table), counting
# distinct hospital_id values among the hospitalizations still in the cohort at
# this point. Counts only; no hospital identifiers leave the site.

hospital_types =
  dplyr::filter(data_list$adt, hospitalization_id %in% cohort_hids) |>
  dplyr::select(hospital_id, hospital_type) |>
  dplyr::distinct() |>
  dplyr::collect() |>
  fgroup_by(hospital_type) |>
  fsummarize(n_hospitals = fnunique(hospital_id)) |>
  ftransform(site = site_lowercase)

fwrite(
  hospital_types,
  here(BOX_DIR, paste0("hospital_types_", site_lowercase, ".csv"))
)

rm(hospital_types)

### tally for inclusion flow diagram -------------------------------------------

fig_s01_01ca = fsubset(cohort, ca_01 == 1) |> fselect(joined_hosp_id) |> fnunique()
fig_s01_01no = fsubset(cohort, ca_01 == 0) |> fselect(joined_hosp_id) |> fnunique()

rm(diagnosis_priority, dx_enc, dx, additional_drops); gc()

## enforce admissions through ED -----------------------------------------------

### identify encounters with ED visit ------------------------------------------

ed_admits =
  dplyr::filter(data_list$adt, hospitalization_id %in% cohort_hids) |>
  dplyr::filter(tolower(location_category) == "ed") |>
  dplyr::select(hospitalization_id) |>
  dplyr::collect()

ed_admits =
  funique(ed_admits) |>
  tibble::deframe()

first_ed_jids =
  fsubset(hid_jid_crosswalk, hospitalization_id %in% ed_admits) |>
  pull(joined_hosp_id)

### apply exclusion to cohort --------------------------------------------------

cohort$ed_admit_01 = if_else(cohort$joined_hosp_id %in% first_ed_jids, 1L, 0L)
fig_s01_02ca       = fsubset(cohort, ca_01 == 1 & ed_admit_01 == 1) |> fselect(joined_hosp_id) |> fnunique()
fig_s01_02no       = fsubset(cohort, ca_01 == 0 & ed_admit_01 == 1) |> fselect(joined_hosp_id) |> fnunique()

## enforce no ICU before touching the wards ------------------------------------

### identify all ICU pre-ward moments ------------------------------------------

icu_before_wards =
  dplyr::filter(data_list$adt, hospitalization_id %in% cohort_hids) |>
  dplyr::filter(tolower(location_category) %in% c("icu", "ward")) |>
  dplyr::collect() |>
  ftransform(location_category = tolower(location_category))

icu_before_wards =
  join(icu_before_wards, hid_jid_crosswalk, how = "inner", multiple = T) |>
  roworder(in_dttm) |>
  fgroup_by(joined_hosp_id, location_category) |>
  fsummarize(min_time = ffirst(in_dttm)) |>
  pivot_wider(
    names_from   = location_category,
    values_from  = min_time,
    names_prefix = "time_"
  ) |>
  fsubset(!is.na(time_icu)) |>
  fsubset(time_icu < time_ward) |>
  pull(joined_hosp_id)

### apply exclusion to cohort --------------------------------------------------

cohort       = fsubset(cohort, !joined_hosp_id %in% icu_before_wards)
cohort_pats  = funique(cohort$patient_id)
cohort_jids  = funique(cohort$joined_hosp_id)
cohort_hids  = funique(hid_jid_crosswalk$hospitalization_id)
date_frame   = select(cohort, patient_id, joined_hosp_id, ends_with("dttm"))
fig_s01_03ca = fsubset(cohort, ca_01 == 1 & ed_admit_01 == 1) |> fselect(joined_hosp_id) |> fnunique()
fig_s01_03no = fsubset(cohort, ca_01 == 0 & ed_admit_01 == 1) |> fselect(joined_hosp_id) |> fnunique()

rm(icu_before_wards, ed_admits, first_ed_jids); gc()

## enforce at least 6h ward data available -------------------------------------

### find the last vital sign measurement time in each encounter ----------------

min_ward_hours = 6L

vmax =
  dplyr::filter(data_list$vitals, hospitalization_id %in% cohort_hids) |>
  dplyr::filter(vital_category %in% req_vitals) |>
  dplyr::select(hospitalization_id, recorded_dttm) |>
  dplyr::collect()

vmax = distinct(vmax)

vmax =
  join(vmax, hid_jid_crosswalk, how = "inner", multiple = T) |>
  roworder(recorded_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(vtime = flast(recorded_dttm))

### is the last vital time within 6h of first ward time? -----------------------

first_ward =
  roworder(ward_times, in_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(first_ward_dttm = ffirst(in_dttm))

vmax =
  join(vmax, first_ward, how = "inner", multiple = F) |>
  fsubset(vtime < first_ward_dttm + lubridate::dhours(min_ward_hours)) |>
  pull(joined_hosp_id)

### apply exclusion to cohort --------------------------------------------------

cohort       = fsubset(cohort, !joined_hosp_id %in% vmax)
cohort       = fsubset(cohort, discharge_dttm >= admission_dttm + lubridate::dhours(min_ward_hours))
cohort_pats  = funique(cohort$patient_id)
cohort_jids  = funique(cohort$joined_hosp_id)
cohort_hids  = funique(hid_jid_crosswalk$hospitalization_id)
date_frame   = select(cohort, patient_id, joined_hosp_id, ends_with("dttm"))
fig_s01_04ca = fsubset(cohort, ca_01 == 1 & ed_admit_01 == 1) |> fselect(joined_hosp_id) |> fnunique()
fig_s01_04no = fsubset(cohort, ca_01 == 0 & ed_admit_01 == 1) |> fselect(joined_hosp_id) |> fnunique()

rm(vmax); gc()

## time starts at the first ward moment ----------------------------------------

cohort = join(cohort, first_ward, how = "left", multiple = F)

## enforce no outcomes before first score --------------------------------------

### early deaths are already covered by no LOS < 6h...

### identify first icu times ---------------------------------------------------

icu =
  dplyr::filter(data_list$adt, hospitalization_id %in% cohort_hids) |>
  dplyr::filter(tolower(location_category) == "icu") |>
  dplyr::select(hospitalization_id, in_dttm) |>
  dplyr::collect()

icu =
  join(icu, hid_jid_crosswalk, how = "inner", multiple = T) |>
  roworder(in_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(itime = ffirst(in_dttm))

### is the first ICU time within 6h of admission_dttm? -------------------------

icu =
  join(icu, cohort, how = "inner", multiple = F) |>
  fsubset(itime < first_ward_dttm + lubridate::dhours(6)) |>
  pull(joined_hosp_id)

### apply exclusion to cohort --------------------------------------------------

cohort       = fsubset(cohort, !joined_hosp_id %in% icu)
cohort       = fsubset(cohort, discharge_dttm >= first_ward_dttm + lubridate::dhours(min_ward_hours))
cohort_pats  = funique(cohort$patient_id)
cohort_jids  = funique(cohort$joined_hosp_id)
cohort_hids  = funique(hid_jid_crosswalk$hospitalization_id)
date_frame   = select(cohort, patient_id, joined_hosp_id, ends_with("dttm"))
fig_s01_05ca = fsubset(cohort, ca_01 == 1 & ed_admit_01 == 1) |> select(joined_hosp_id) |> fnunique()
fig_s01_05no = fsubset(cohort, ca_01 == 0 & ed_admit_01 == 1) |> select(joined_hosp_id) |> fnunique()

rm(icu); gc()

## save data for inclusion flow diagram ----------------------------------------

step_labels = c(
  "Adult inpatient admissions during study period",
  "After excluding patients not admitted through the ED",
  "After excluding patients who were in the ICU before hitting the wards",
  "After excluding encounters with < 6h data",
  "After excluding encounters with outcomes too early"
)

n_remaining_ca = c(fig_s01_01ca, fig_s01_02ca, fig_s01_03ca, fig_s01_04ca, fig_s01_05ca)
n_remaining_no = c(fig_s01_01no, fig_s01_02no, fig_s01_03no, fig_s01_04no, fig_s01_05no)
n_excluded_ca  = c(NA, diff(n_remaining_ca) * -1)
n_excluded_no  = c(NA, diff(n_remaining_no) * -1)

flow_df = tidytable(
  step = step_labels,
  n_remaining_ca,
  n_excluded_ca,
  n_remaining_no,
  n_excluded_no,
)

fwrite(flow_df, here(BOX_DIR, paste0("figure_s01_flow_", site_lowercase, ".csv")))

rm(flow_df, step_labels, n_remaining_ca, n_remaining_no, n_excluded_ca, n_excluded_no)
gc()

# prepare additional cohort details --------------------------------------------

## patient demographics --------------------------------------------------------

### pull demographics from arrow table -----------------------------------------

cohort_demographics =
  dplyr::filter(data_list$patient, patient_id %in% cohort_pats) |>
  dplyr::select(patient_id, death_dttm, ends_with("category")) |>
  dplyr::collect() |>
  distinct()

pt_dups =
  fcount(cohort_demographics, patient_id, name = "n") |>
  fsubset(n > 1)

if (nrow(pt_dups) > 0) {
  stop(
    sprintf("Duplicate patient_id(s): %s", paste(pt_dups$patient_id, collapse = ", ")), call. = F
  )
}

### add demographics to cohort df ----------------------------------------------

cohort =
  join(cohort, cohort_demographics, how = "left", multiple = F) |>
  fmutate(miss_age   = as.integer(is.na(age))) |>
  fmutate(miss_sex   = as.integer(is.na(sex_category)       | tolower(sex_category)       == "unknown")) |>
  fmutate(miss_race  = as.integer(is.na(race_category)      | tolower(race_category)      == "unknown")) |>
  fmutate(miss_eth   = as.integer(is.na(ethnicity_category) | tolower(ethnicity_category) == "unknown")) |>
  fmutate(age        = if_else(age > 90, 90.9, age)) |>
  fmutate(female_01  = if_else(tolower(sex_category) == "female", 1L, 0L)) |>
  fmutate(dead_01    = if_else(tolower(discharge_category) == "expired", 1L, 0L)) |>
  fmutate(hospice_01 = if_else(tolower(discharge_category) == "hospice", 1L, 0L)) |>
  fmutate(los_hosp_d = as.numeric(difftime(discharge_dttm, admission_dttm), "hours")/24) |>
  mutate(across(
    .cols = where(is.character) & !any_of("patient_id"),
    .fns  = ~tolower(.x)
  )) |>
  fselect(-sex_category, -discharge_category)

rm(pt_dups, cohort_demographics); gc()

## admission code status -------------------------------------------------------

code_window_hours = 12L

codes =
  data_list[["code_status"]] |>
  dplyr::select(patient_id, start_dttm, code_status_category) |>
  dplyr::collect() |>
  distinct()

codes =
  join(cohort, codes, how = "left", multiple = T) |>
  fsubset(start_dttm >= admission_dttm - lubridate::ddays(1)) |>
  fsubset(start_dttm <= admission_dttm + lubridate::dhours(code_window_hours)) |>
  roworder(start_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(initial_code_status = flast(code_status_category))

cohort =
  join(cohort, codes, how = "left", multiple = F) |>
  ftransform(
    initial_code_status = case_when(
      tolower(initial_code_status) == "dnr" ~ "other",
      is.na(initial_code_status)            ~ "presume_full",
      TRUE                                  ~ tolower(initial_code_status)
    )
  )

rm(codes, code_window_hours); gc()

# outcomes ---------------------------------------------------------------------

## icu admissions --------------------------------------------------------------

### icu data from arrow table --------------------------------------------------

icu =
  dplyr::filter(data_list$adt, hospitalization_id %in% cohort_hids) |>
  dplyr::select(hospitalization_id, in_dttm, location_category) |>
  dplyr::arrange(in_dttm) |>
  dplyr::collect()

icu = funique(icu)

### link to joined_hosp_id -----------------------------------------------------

icu =
  join(icu, hid_jid_crosswalk, how = "left",  multiple = T) |>
  fselect(joined_hosp_id, in_dttm, location_category) |>
  funique() |>
  join(date_frame, how = "inner", multiple = T) |>
  fsubset(in_dttm >= admission_dttm & in_dttm <= discharge_dttm) |>
  fmutate(location_category = tolower(location_category))

### set aside all icu encounters -----------------------------------------------

icu_encs =
  fsubset(icu, location_category == "icu") |>
  fselect(joined_hosp_id) |>
  funique() |>
  tibble::deframe()

### identify ward-icu transfer moments -----------------------------------------

icu =
  roworder(icu, in_dttm) |>
  # keep procedural locations so Ward→Procedural→ICU is not falsely flagged
  fsubset(tolower(location_category) %in% c('icu', 'ward', 'procedural', 'ed')) |>
  group_by(joined_hosp_id) |>
  mutate(prev_loc = lag(location_category)) |>
  mutate(ward_icu = tolower(location_category) == "icu" & tolower(prev_loc) == "ward") |>
  ungroup() |>
  fsubset(ward_icu)

### first ward-icu transfer ----------------------------------------------------

icu =
  roworder(icu, in_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(event_dttm = ffirst(in_dttm))  |>
  ftransform(event = "icu")

## death -----------------------------------------------------------------------

death =
  fsubset(cohort, dead_01 == 1) |>
  fselect(joined_hosp_id, event_dttm = discharge_dttm) |>
  ftransform(event = "death")

#### don't need to specify death on wards,
#### because if it occurs before ICU it must be

## hospice ---------------------------------------------------------------------

hospice =
  fsubset(cohort, hospice_01 == 1) |>
  fselect(joined_hosp_id, event_dttm = discharge_dttm) |>
  fmutate(event = "hospice")

## outcomes data frame ---------------------------------------------------------

### event-level frame (one row per encounter with >= 1 event) ------------------

outcome_events =
  rowbind(icu, death, hospice) |>
  pivot_wider(
    names_from   = event,
    values_from  = event_dttm,
    names_prefix = "time_"
  )

### guarantee that all three event columns exist -------------------------------
# pivot_wider() drops an event type entirely when a site records none of it,
# which breaks the pmin() calls below. Creating the column as all-missing leaves
# every value unchanged at sites that do record the event.

time_cols_all  = c("time_icu", "time_death", "time_hospice")
time_cols_have = intersect(time_cols_all, names(outcome_events))
time_cols_miss = setdiff(time_cols_all,   names(outcome_events))

if (length(time_cols_have) == 0) {
  stop("No ICU, death, or hospice events found; cannot build the outcome frame.", call. = FALSE)
}

if (length(time_cols_miss) > 0) {

  outcome_events = as.data.table(outcome_events)
  na_dttm_vec    = outcome_events[[time_cols_have[1]]]
  na_dttm_vec[]  = NA

  for (mc in time_cols_miss) {
    data.table::set(outcome_events, j = mc, value = na_dttm_vec)
  }

  outcome_events = as_tidytable(outcome_events)

  message(
    sprintf("Note: no %s event(s) at this site; time column(s) created as missing.",
            paste(str_remove(time_cols_miss, "^time_"), collapse = ", "))
  )
}

### round-one composite fields (definitions unchanged) -------------------------

outcome_events =
  fmutate(outcome_events, end_enc = pmin(time_death, time_hospice, na.rm = T)) |>
  fmutate(outcome_dttm         = pmin(time_icu, end_enc,    na.rm = T)) |>
  fmutate(outcome_nohospc_dttm = pmin(time_icu, time_death, na.rm = T)) |>
  fmutate(
    outcome_cat = case_when(
      outcome_dttm == time_icu     ~ "icu",
      outcome_dttm == time_death   ~ "death",
      outcome_dttm == time_hospice ~ "hospice",
      TRUE                         ~ "problem"
    )
  ) |>
  fselect(
    joined_hosp_id,
    time_icu,
    time_death,
    time_hospice,
    outcome_dttm,
    outcome_nohospc_dttm,
    outcome_cat
  )

### expand to one row per cohort encounter -------------------------------------
# Round-two change, and a deliberate one. Through round one this file held only
# encounters with at least one of ICU transfer, death, or hospice discharge. It
# now holds every encounter in the final cohort, because the file is the source
# of the o_*_01 indicators and 02_scores.R attaches it with a LEFT join. Without
# the full spine, those indicators would arrive missing rather than zero for
# event-free encounters, which would silently remove those encounters from the
# hospice-decomposition denominators. Spine-only rows carry missing times and a
# missing outcome_cat, so the "problem" level of outcome_cat remains a genuine
# error sentinel rather than the routine label for an event-free encounter.

df_outcomes =
  fselect(cohort, joined_hosp_id) |>
  funique() |>
  join(outcome_events, how = "left", multiple = FALSE)

### per-outcome event and censoring times --------------------------------------
# The universe of competing events is {icu, death, hospice}. For outcome k with
# event set E_k and complement C_k, let t_e be the earliest time in E_k and t_c
# the earliest time in C_k. Then
#
#   event_k_dttm  = t_e when t_e is present and strictly precedes t_c
#                   (an absent t_c never blocks the event; this covers the
#                    composite, whose complement is structurally empty)
#   censor_k_dttm = t_c when t_c is present and strictly precedes t_e
#                   (an absent t_e never blocks the censoring; this is what
#                    censors, say, a ward death in the ward-to-ICU outcome)
#   o_k_01        = 1 when event_k_dttm is present
#
# Both rules are strict, so an exact tie leaves event and censoring times both
# missing. That is what makes the QC assertion below (never both present) true
# by construction. Ties are counted and reported rather than passed over.

tie_counts = integer(0)

#### composite: events = icu, death, hospice | complement = none ---------------

t_e    = pmin(df_outcomes$time_icu, df_outcomes$time_death, df_outcomes$time_hospice, na.rm = TRUE)
t_c    = copy(t_e)
t_c[]  = NA
keep_e = !is.na(t_e) & (is.na(t_c) | t_e < t_c)
keep_c = !is.na(t_c) & (is.na(t_e) | t_c < t_e)
n_ties = sum(!is.na(t_e) & !is.na(t_c) & t_e == t_c)

t_e[!keep_e] = NA
t_c[!keep_c] = NA

df_outcomes$event_composite_dttm  = t_e
df_outcomes$censor_composite_dttm = t_c
df_outcomes$o_composite_01        = as.integer(!is.na(t_e))
tie_counts["composite"]           = n_ties

#### nohospice: events = icu, death | complement = hospice ---------------------

t_e    = pmin(df_outcomes$time_icu, df_outcomes$time_death, na.rm = TRUE)
t_c    = copy(df_outcomes$time_hospice)
keep_e = !is.na(t_e) & (is.na(t_c) | t_e < t_c)
keep_c = !is.na(t_c) & (is.na(t_e) | t_c < t_e)
n_ties = sum(!is.na(t_e) & !is.na(t_c) & t_e == t_c)

t_e[!keep_e] = NA
t_c[!keep_c] = NA

df_outcomes$event_nohospice_dttm  = t_e
df_outcomes$censor_nohospice_dttm = t_c
df_outcomes$o_nohospice_01        = as.integer(!is.na(t_e))
tie_counts["nohospice"]           = n_ties

#### wardicu: events = icu | complement = death, hospice -----------------------

t_e    = copy(df_outcomes$time_icu)
t_c    = pmin(df_outcomes$time_death, df_outcomes$time_hospice, na.rm = TRUE)
keep_e = !is.na(t_e) & (is.na(t_c) | t_e < t_c)
keep_c = !is.na(t_c) & (is.na(t_e) | t_c < t_e)
n_ties = sum(!is.na(t_e) & !is.na(t_c) & t_e == t_c)

t_e[!keep_e] = NA
t_c[!keep_c] = NA

df_outcomes$event_wardicu_dttm  = t_e
df_outcomes$censor_wardicu_dttm = t_c
df_outcomes$o_wardicu_01        = as.integer(!is.na(t_e))
tie_counts["wardicu"]           = n_ties

#### warddeath: events = death | complement = icu, hospice ---------------------

t_e    = copy(df_outcomes$time_death)
t_c    = pmin(df_outcomes$time_icu, df_outcomes$time_hospice, na.rm = TRUE)
keep_e = !is.na(t_e) & (is.na(t_c) | t_e < t_c)
keep_c = !is.na(t_c) & (is.na(t_e) | t_c < t_e)
n_ties = sum(!is.na(t_e) & !is.na(t_c) & t_e == t_c)

t_e[!keep_e] = NA
t_c[!keep_c] = NA

df_outcomes$event_warddeath_dttm  = t_e
df_outcomes$censor_warddeath_dttm = t_c
df_outcomes$o_warddeath_01        = as.integer(!is.na(t_e))
tie_counts["warddeath"]           = n_ties

#### hospicedc: events = hospice | complement = icu, death ---------------------

t_e    = copy(df_outcomes$time_hospice)
t_c    = pmin(df_outcomes$time_icu, df_outcomes$time_death, na.rm = TRUE)
keep_e = !is.na(t_e) & (is.na(t_c) | t_e < t_c)
keep_c = !is.na(t_c) & (is.na(t_e) | t_c < t_e)
n_ties = sum(!is.na(t_e) & !is.na(t_c) & t_e == t_c)

t_e[!keep_e] = NA
t_c[!keep_c] = NA

df_outcomes$event_hospicedc_dttm  = t_e
df_outcomes$censor_hospicedc_dttm = t_c
df_outcomes$o_hospicedc_01        = as.integer(!is.na(t_e))
tie_counts["hospicedc"]           = n_ties

## quality control on outcome times --------------------------------------------

### the composite must reproduce the round-one outcome_dttm rate exactly -------

qc_composite_ref = as.integer(!is.na(df_outcomes$outcome_dttm))

if (!identical(df_outcomes$o_composite_01, qc_composite_ref)) {
  stop(
    sprintf("o_composite_01 disagrees with outcome_dttm for %d encounter(s).",
            sum(df_outcomes$o_composite_01 != qc_composite_ref)),
    call. = FALSE
  )
}

### event and censoring times must never both be present -----------------------

outcome_keys = c("composite", "nohospice", "wardicu", "warddeath", "hospicedc")

for (k in outcome_keys) {

  n_both = sum(
    !is.na(df_outcomes[[paste0("event_",  k, "_dttm")]]) &
    !is.na(df_outcomes[[paste0("censor_", k, "_dttm")]])
  )

  if (n_both > 0) {
    stop(
      sprintf("Outcome '%s' has %d encounter(s) carrying both an event and a censoring time.",
              k, n_both),
      call. = FALSE
    )
  }
}

message(
  sprintf("✅ Outcome times QC passed | %s encounters | exact ties: %s",
          format(nrow(df_outcomes), big.mark = ","),
          paste(sprintf("%s=%d", names(tie_counts), tie_counts), collapse = ", "))
)

write_parquet(df_outcomes, here("proj_tables", "outcome_times.parquet"))

ward_icu_tx =
  fsubset(df_outcomes, !is.na(outcome_cat) & outcome_cat == "icu") |>
  pull(joined_hosp_id)

rm(death, hospice, icu, outcome_events)
rm(t_e, t_c, keep_e, keep_c, n_ties, tie_counts, qc_composite_ref, outcome_keys, k, n_both)
rm(time_cols_all, time_cols_have, time_cols_miss)
suppressWarnings(rm(na_dttm_vec, mc)); gc()

# other care processes --------------------------------------------------------#

### vasopressors ---------------------------------------------------------------

va_list = c(
  "norepinephrine",
  "vasopressin",
  "phenylephrine",
  "epinephrine",
  "dopamine"
)

meds =
  dplyr::filter(data_list$medication_admin_continuous, tolower(med_category) %in% va_list) |>
  dplyr::select(hospitalization_id, admin_dttm, ends_with("category")) |>
  dplyr::collect()

meds =
  join(meds, hid_jid_crosswalk, how = "inner", multiple = T) |>
  fselect(joined_hosp_id, admin_dttm, med_category) |>
  distinct()

meds =
  join(meds, date_frame, how = "inner", multiple = T) |>
  fsubset(admin_dttm >= admission_dttm) |>
  fsubset(admin_dttm <= discharge_dttm) |>
  fmutate(med_category = tolower(med_category)) |>
  fmutate(event = "vasopressor") |>
  fselect(joined_hosp_id, event_dttm = admin_dttm, event, med_category) |>
  distinct()

rm(va_list); gc()

### respiratory support --------------------------------------------------------

resp_list = c(
  "imv",
  "nippv",
  "high flow nc"
)

resp =
  dplyr::filter(data_list$respiratory_support, hospitalization_id %in% cohort_hids) |>
  dplyr::filter(tolower(device_category) %in% resp_list) |>
  dplyr::select(hospitalization_id, recorded_dttm, device_category) |>
  dplyr::collect()

resp =
  join(resp, hid_jid_crosswalk, how = "inner", multiple = T) |>
  fselect(joined_hosp_id, recorded_dttm, device_category) |>
  distinct()

resp =
  join(resp, date_frame, how = "inner", multiple = T) |>
  fsubset(recorded_dttm >= admission_dttm) |>
  fsubset(recorded_dttm <= discharge_dttm) |>
  fmutate(event = tolower(device_category)) |>
  fselect(joined_hosp_id, event_dttm = recorded_dttm, event) |>
  distinct()

present_events = sort(funique(tolower(resp$event)))
missing_events = setdiff(resp_list, present_events)

if (length(missing_events) > 0) {
  stop(
    sprintf("Resp need at least one %s, but %s not found.",
            if (length(missing_events) == 1) "event" else "events",
            paste(missing_events, collapse = ", ")
    ),
    call. = FALSE
  )
}

rm(resp_list, present_events, missing_events)
gc()

### combine and save -----------------------------------------------------------

rowbind(meds, resp, fill = T) |>
  write_parquet(here("proj_tables", "careprocess.parquet"))

## cohort (1 row per encounter) ------------------------------------------------

va_encs  =
  fsubset(meds, event == "vasopressor") |>
  fselect(joined_hosp_id) |>
  funique() |>
  tibble::deframe()

imv_encs =
  fsubset(resp, event == "imv") |>
  fselect(joined_hosp_id) |>
  funique() |>
  tibble::deframe()

cohort =
  funique(cohort) |>
  fmutate(wicu_01    = if_else(joined_hosp_id %in% ward_icu_tx,  1L, 0L, 0L)) |>
  fmutate(icu_01     = if_else(joined_hosp_id %in% icu_encs,     1L, 0L, 0L)) |>
  fmutate(imv_01     = if_else(joined_hosp_id %in% imv_encs,     1L, 0L, 0L)) |>
  fmutate(va_01      = if_else(joined_hosp_id %in% va_encs,      1L, 0L, 0L)) |>
  fmutate(d_noicu_01 = if_else(dead_01 == 1    & icu_01 != 1, 1L, 0L, 0L)) |>
  fmutate(h_noicu_01 = if_else(hospice_01 == 1 & icu_01 != 1, 1L, 0L, 0L)) |>
  # Metastatic flag, restricted to solid tumors on purpose. liquid_01_enc is
  # assigned by fmax() across every cancer code on the encounter rather than by
  # the diagnosis hierarchy, so a hematologic encounter can also carry a C77-C80
  # code; leaving those encounters missing keeps the flag interpretable as a
  # solid-tumor stage contrast. rank_enc: 1 metastatic, 2 hematologic, 3
  # high-risk solid (C22/C25/C34), 4 other solid (C18/C50/C61), 5 other.
  fmutate(
    mets_01 = case_when(
      rank_enc == 1L & liquid_01 == 0L                 ~ 1L,
      ca_01 == 1L & liquid_01 == 0L & rank_enc != 1L   ~ 0L,
      TRUE                                             ~ NA_integer_
    )
  ) |>
  select(
    patient_id,
    joined_hosp_id,
    admission_dttm,
    discharge_dttm,
    first_ward_dttm,
    age,
    race_category,
    ethnicity_category,
    ends_with("01"),
    rank_enc,
    initial_code_status,
    los_hosp_d,
    starts_with("miss")
  )

## mets_01 is deliberately three-valued (1 / 0 / missing), so it is held out of
## the sweep that fills missing binary flags with zero.

fill_01_cols = setdiff(grep("01$", names(cohort), value = TRUE), "mets_01")

cohort =
  mutate(cohort, across(
    .cols = ends_with("category"),
    .fns  = ~if_else(is.na(.x), "unknown", tolower(.x))
  )) |>
  mutate(across(
    .cols = all_of(fill_01_cols),
    .fns  = ~if_else(is.na(.x), 0L, .x)
  ))

## quality control on the metastatic flag --------------------------------------

qc_mets_noca  = sum(cohort$ca_01 == 0 & !is.na(cohort$mets_01))
qc_mets_solid = sum(cohort$ca_01 == 1 & cohort$liquid_01 == 0 & is.na(cohort$mets_01))

if (qc_mets_noca > 0 || qc_mets_solid > 0) {
  stop(
    sprintf("mets_01 failed QC: %d non-cancer encounter(s) non-missing, %d solid-tumor encounter(s) missing.",
            qc_mets_noca, qc_mets_solid),
    call. = FALSE
  )
}

message("✅ mets_01 QC passed.")

rm(fill_01_cols, qc_mets_noca, qc_mets_solid)

## add elixhauser --------------------------------------------------------------

### vector common codes to avoid for RAM's sake (to repl w/ keep list) ---------

unused_vect = c(
  "Z79.899",
  "E78.5",
  "Z87.891",
  "K21.9",
  "Z20.822",
  "F17.210",
  "Z79.01",
  "N17.9"
)

### only relevant codes from relevant encounters -------------------------------

elix =
  dplyr::filter(data_list$hospital_diagnosis, hospitalization_id %in% cohort_hids) |>
  dplyr::filter(!toupper(diagnosis_code) %in% unused_vect) |>
  dplyr::select(hospitalization_id, poa_present, diagnosis_code) |>
  dplyr::collect() |>
  distinct()

### assign elixhauser diagnosis dummies ----------------------------------------

elix =
  comorbidity::comorbidity(
    elix,
    id      = "hospitalization_id",
    code    = "diagnosis_code",
    map     = "elixhauser_icd10_quan",
    assign0 = T
  )

### link to joined-hosp-id -----------------------------------------------------

elix =
  join(elix, hid_jid_crosswalk, how = "left", multiple = T) |>
  fselect(-patient_id, -hospitalization_id) |>
  fgroup_by(joined_hosp_id) |>
  fmax()

### calculate vw scores --------------------------------------------------------

vw        = comorbidity::score(elix, weights = "vw", assign0 = T)
elix      = cbind(elix, vw = vw) |> fselect(joined_hosp_id, vw)
cohort    = join(cohort, elix, how = "left", multiple = F)

cohort$miss_vw = as.integer(is.na(cohort$vw))
cohort$vw      = if_else(is.na(cohort$vw), 0L, cohort$vw)

rm(elix, vw); gc()

## add hospital_id -------------------------------------------------------------

hospital =
  dplyr::select(data_list$adt, hospitalization_id, in_dttm, hospital_id) |>
  dplyr::filter(hospitalization_id %in% cohort_hids) |>
  dplyr::collect()

hospital =
  join(hospital, hid_jid_crosswalk, how = "left", multiple = T) |>
  roworder(in_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(hospital_id = ffirst(hospital_id))

cohort = join(cohort, hospital, how = "left", multiple = F)

rm(hospital); gc()

## restrict handoff objects to the final cohort --------------------------------
# The exclusions above (ED, ICU-before-ward, <6h ward data, early outcomes) are
# applied to `cohort`, but hid_jid_crosswalk / cohort_hids / ward_times were built
# earlier and never re-pruned, so they still carry excluded encounters. Prune them
# to the final cohort here — before cohort.parquet and hid_jid_crosswalk.parquet
# are written and before the objects are handed to 02_scores.R — so scoring covers
# exactly the analytic set. Without this, 02 would score a superset of encounters,
# scores_full would hold ~20% extra rows, and those rows (absent from the
# outcome_times spine) would fail the o_composite_01 carry-through QC in 02 on a
# missing-vs-zero mismatch. outcome_times above is already built from the final
# `cohort`, so it needs no pruning. Downstream artifacts are unchanged: 03a
# already inner-joins scores to the cohort, dropping these encounters regardless.

final_jids        = funique(cohort$joined_hosp_id)
hid_jid_crosswalk = fsubset(hid_jid_crosswalk, joined_hosp_id %in% final_jids)
ward_times        = fsubset(ward_times,        joined_hosp_id %in% final_jids)
cohort_hids       = funique(hid_jid_crosswalk$hospitalization_id)
cohort_jids       = final_jids
rm(final_jids)

## merge ward intervals into disjoint blocks per encounter ---------------------
# Up to this point ward_times holds one row per ADT ward stay, so an encounter
# has several rows and adjacent stays share a boundary timestamp: a ward-to-ward
# transfer emits out_dttm on the first row equal to in_dttm on the second. The
# ward restriction in 02_scores.R and 02b_carryforward.R is a closed-interval
# filter, `time >= in_dttm & time <= out_dttm`, so a score computed at exactly
# that shared instant satisfies both rows and the encounter carries the same
# (joined_hosp_id, time) twice. The same happens whenever two ward intervals
# genuinely overlap, which this pipeline can produce on its own: it reclassifies
# stepdown to ward across the whole ADT, so a site emitting both a unit-level
# ward row and an overlapping bed-level stepdown row ends up with two ward
# intervals over one span. Downstream the duplicate was removed only incidentally,
# by the terminal funique() in the ward-restriction block, after in_dttm and
# out_dttm had been overwritten with encounter-level bounds that happened to make
# the two rows identical.
#
# Merging overlapping and touching intervals into disjoint blocks removes the
# duplicate at its source. The union of covered time is unchanged by construction,
# so no score row is gained or lost; only the multiplicity changes. Encounters
# with a genuine off-ward gap (ward, then ICU, then ward) keep separate blocks,
# because those intervals neither overlap nor touch. hospitalization_id is dropped
# here: a merged block can span two constituent hospitalizations, and no consumer
# of ward_times uses the column -- both 02_scores.R and 02b_carryforward.R deleted
# it immediately after the join.

ward_blocks =
  fselect(ward_times, joined_hosp_id, in_dttm, out_dttm) |>
  roworder(joined_hosp_id, in_dttm, out_dttm) |>
  as.data.table()

ward_blocks[, prior_max_out := shift(cummax(as.numeric(out_dttm))), by = joined_hosp_id]
ward_blocks[, new_block     := is.na(prior_max_out) | as.numeric(in_dttm) > prior_max_out]
ward_blocks[, block_id      := cumsum(new_block)]

n_ward_rows_raw = nrow(ward_blocks)

ward_times =
  fgroup_by(ward_blocks, joined_hosp_id, block_id) |>
  fsummarize(
    in_dttm  = fmin(in_dttm),
    out_dttm = fmax(out_dttm)
  ) |>
  fselect(-block_id) |>
  roworder(joined_hosp_id, in_dttm)

message(sprintf("Ward intervals merged | %s ADT rows -> %s disjoint blocks (%s collapsed).",
                format(n_ward_rows_raw,  big.mark = ","),
                format(nrow(ward_times), big.mark = ","),
                format(n_ward_rows_raw - nrow(ward_times), big.mark = ",")))

rm(ward_blocks, n_ward_rows_raw)

message(sprintf("Cohort restriction | %s encounters | %s hospitalizations kept.",
                format(fnunique(cohort$joined_hosp_id), big.mark = ","),
                format(length(cohort_hids),             big.mark = ",")))

## sanity check before saving --------------------------------------------------

props =
  tidytable(
    icu     = fmean(cohort$icu_01,     na.rm=TRUE),
    dead    = fmean(cohort$dead_01,    na.rm=TRUE),
    hospice = fmean(cohort$hospice_01, na.rm=TRUE),
    imv     = fmean(cohort$imv_01,     na.rm=TRUE),
    va      = fmean(cohort$va_01,      na.rm=TRUE)
  )

if (
  props$icu     > 0.50 | props$icu     == 0 |
  props$dead    > 0.20 | props$dead    == 0 |
  props$hospice > 0.10 | props$hospice == 0 |
  props$imv     > 0.40 | props$imv     == 0 |
  props$va      > 0.40 | props$va      == 0) {
  stop("Sanity check failed: Outcome distribution out of expected range.")
}

## add admission year, retained for the descriptive tables built in 02 ---------
# The encounter-level descriptive artifacts (missing_demog, table_02_cat,
# table_02_cont, admission_dx_chapter, admission_dx_stem) moved to the end of
# 02_scores.R, where they run on the cohort AFTER the drop of encounters that
# carry no calculable score. Built here, they ran on the pre-reconciliation set
# and reported a denominator no analysis artifact could reproduce. Two
# consequences for this file: `year` is added before the write so the relocated
# Table 2 still has it, and the miss_* columns are no longer dropped, because
# missing_demog now reads them from cohort.parquet in the next stage.

cohort = fmutate(cohort, year = lubridate::year(admission_dttm))

write_parquet(cohort,            here("proj_tables", "cohort.parquet"))
write_parquet(hid_jid_crosswalk, here("proj_tables", "hid_jid_crosswalk.parquet"))


# final cleanup ----------------------------------------------------------------

keep = c(
  "BOX_DIR",
  "data_list",
  "site_lowercase",
  "allow_sparse_o2",
  "cohort",
  "cohort_hids",
  "cohort_jids",
  "cohort_pats",
  "hid_jid_crosswalk",
  "ward_times",
  "req_vitals",
  "req_labs",
  "start_time",
  "run_log"
)

rm(list = setdiff(ls(), keep)); gc()

# go to 02
