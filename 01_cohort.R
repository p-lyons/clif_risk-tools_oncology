

# Cohort script for CLIF project validating risk tools in oncology.
# Requires data_list to be loaded/validated from 00_*

# check/create output directories ---------------------------------------------

if (!exists("project_location")) {
  stop("project_location not found. Did you run 00_setup first?", call. = FALSE)
}

if (!dir.exists(paste0(project_location, "/upload_to_box"))) {
  dir.create(paste0(project_location, "/upload_to_box"), recursive = TRUE)
}

if (!dir.exists(paste0(project_location, "/proj_tables"))) {
  dir.create(paste0(project_location, "/proj_tables"), recursive = TRUE)
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

# cohort identification --------------------------------------------------------

## start by linking contiguous hospitalizations --------------------------------

### encounters with age >= 18 and dates within study window --------------------

hosp_blocks = 
  dplyr::filter(data_list$hospitalization, age_at_admission >= 18) |>
  dplyr::filter(admission_dttm >= start_date & admission_dttm <= end_date) |> 
  dplyr::filter(admission_dttm < discharge_dttm & !is.na(discharge_dttm)) |>
  dplyr::select(patient_id, hospitalization_id, admission_dttm, discharge_dttm) |>
  dplyr::collect() |>
  roworder(patient_id, admission_dttm)

### use data.table to find joined hospitalizations with <= 6h gaps -------------

link_hours = 6L
linked     = as.data.table(hosp_blocks)
setorder(linked, patient_id, admission_dttm)

#### calculate gaps between encounters
linked[, next_admit := shift(admission_dttm, type = "lead"), by = patient_id]
linked[, next_gap   := as.numeric(difftime(next_admit, discharge_dttm, units = "hours"))]
linked[, prev_dc    := shift(discharge_dttm, type = "lag"), by = patient_id]  
linked[, prev_gap   := as.numeric(difftime(admission_dttm, prev_dc, units = "hours"))]

#### mark encounters that should be linked
linked[, link_flag := (next_gap < link_hours | prev_gap < link_hours)]
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
  dplyr::collect() 

has_vital_signs = distinct(has_vital_signs) 

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
  dplyr::collect()

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
  fsubset(cohort, discharge_category == "Expired") |>
  roworder(admission_dttm, discharge_dttm) |>
  fgroup_by(patient_id) |>
  fmutate(one = 1L) |>
  fmutate(n_deaths = fsum(one), counter = fcumsum(one)) |>
  fungroup() |>
  fsubset(n_deaths > 1) 

if (nrow(dup_deaths) > 0) {
  
  encdrop = 
    fsubset(dup_deaths, counter == 1) |>
    pull(joined_hosp_id) 

  cohort  = fsubset(cohort, !joined_hosp_id %in% encdrop)
  n_dupes = sum(dup_deaths$counter > 1)
  n_pats  = length(unique(dup_deaths$patient_id))
  
  message(
    sprintf("Removed %d duplicate deaths for %d patients: %s",
            n_dupes, n_pats, paste(unique(dup_deaths$patient_id), collapse = ", "))
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
    sprintf("Removed %d post-death encs for %d patients (? organ donors): %s",
            nrow(post_death_admissions),
            length(unique(post_death_admissions$patient_id)),
            paste(unique(post_death_admissions$patient_id), collapse = ", "))
  )
}

message("✅ Cleaned duplicate deaths and post-death encounters.")

rm(dupes, dup_deaths, death_times, post_death_admissions, start_date, end_date)
rm(encdrop, file_type, n_dupes, n_pats); gc()

cohort_pats = funique(cohort$patient_id)
cohort_jids = funique(cohort$joined_hosp_id)
cohort_hids = funique(hid_jid_crosswalk$hospitalization_id)
date_frame  = select(cohort, patient_id, joined_hosp_id, ends_with("dttm"))

## identify cancer patients ----------------------------------------------------

### load cancer dx codes -------------------------------------------------------

ca_codes = 
  readxl::read_xlsx(here("config/icd10cm_casefinding_2023.xlsx")) |> #### exclude C44, in remission 
  janitor::clean_names() |>
  fsubset(general_category %in% c("Hematopoietic neoplasm", "Malignant neoplasm")) |>
  fsubset(is.na(drop) | drop != "drop") |>
  fsubset(!str_detect(icd_10_cm_code_specific, "^Z85")) |>
  fsubset(!str_detect(icd_10_cm_code_definition, "in remission|personal history")) |>
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
    liquid_01_enc = ffirst(liquid_01),
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
  here("upload_to_box", paste0("cancer_codes_primary_", site_lowercase, ".csv"))
)

### tally for inclusion flow diagram -------------------------------------------

fig_s01_01ca = fsubset(cohort, ca_01 == 1) |> fselect(joined_hosp_id) |> fnunique()
fig_s01_01no = fsubset(cohort, ca_01 == 0) |> fselect(joined_hosp_id) |> fnunique()

rm(diagnosis_priority, dx_enc, dx); gc()

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

fwrite(flow_df, here("upload_to_box", paste0("figure_s01_flow_", site_lowercase, ".csv")))

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

codes = 
  data_list[["code_status"]] |>
  dplyr::select(patient_id, start_dttm, code_status_category) |>
  dplyr::collect() |>
  distinct()

codes = 
  join(cohort, codes, how = "left", multiple = T) |>
  fsubset(start_dttm >= admission_dttm - lubridate::ddays(1)) |>
  fsubset(start_dttm <= discharge_dttm) |>
  roworder(start_dttm) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(initial_code_status = ffirst(code_status_category)) |>
  fmutate(
    initial_code_status = case_when(
      tolower(initial_code_status) == "dnr" ~ "other", 
      is.na(initial_code_status)            ~ "unspecified",
      TRUE                                  ~ tolower(initial_code_status)
    )
  )

cohort = join(cohort, codes, how = "left", multiple = F)

rm(codes); gc()

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
  # avoid categorizing radiology/dialysis/other as the event before ward
  fsubset(tolower(location_category) %in% c('icu', 'ward')) |>
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

df_outcomes = 
  rowbind(icu, death, hospice) |>
  pivot_wider(
    names_from   = event,
    values_from  = event_dttm,
    names_prefix = "time_"
  ) |>
  fmutate(end_enc              = pmin(time_death, time_hospice, na.rm = T)) |>
  fmutate(outcome_dttm         = pmin(time_icu,   end_enc,      na.rm = T)) |>
  fmutate(outcome_nohospc_dttm = pmin(time_icu,   time_death,   na.rm = T)) |>
  fmutate(
    outcome_cat = case_when(
      outcome_dttm == time_icu     ~ "icu",
      outcome_dttm == time_death   ~ "death",
      outcome_dttm == time_hospice ~ "hospice",
      TRUE                         ~ "problem"
    )
  ) |>
  select(joined_hosp_id, starts_with("outcome"))

fwrite(df_outcomes, here("proj_tables", "outcome_times.csv"))

ward_icu_tx = 
  fsubset(df_outcomes, outcome_cat == "icu") |>
  pull(joined_hosp_id)

rm(death, hospice, icu); gc()

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
  fmutate(wicu_01 = if_else(joined_hosp_id %in% ward_icu_tx, 1L, 0L, 0L)) |>
  fmutate(icu_01  = if_else(joined_hosp_id %in% icu_encs,    1L, 0L, 0L)) |>
  fmutate(imv_01  = if_else(joined_hosp_id %in% imv_encs,    1L, 0L, 0L)) |>
  fmutate(va_01   = if_else(joined_hosp_id %in% va_encs,     1L, 0L, 0L)) |>
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
    initial_code_status,
    los_hosp_d
  ) |>
  mutate(across(
    .cols = ends_with("category"),
    .fns  = ~if_else(is.na(.x), "unknown", tolower(.x))
  )) |> 
  mutate(across(
    .cols = ends_with("01"),
    .fns  = ~if_else(is.na(.x), 0L, .x)
  ))

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
cohort$vw = if_else(is.na(cohort$vw), 0L, cohort$vw)

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

write_parquet(cohort,            here("proj_tables", "cohort.parquet"))
write_parquet(hid_jid_crosswalk, here("proj_tables", "hid_jid_crosswalk.parquet"))

# t02. characteristics/outcomes by cancer status (0 = none, 1 = cancer) --------

## prepare table component = continuous variables ------------------------------

t2_cont = 
  fsubset(cohort, ed_admit_01 == 1) |>
  fgroup_by(ca_01) |>
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

age_breaks = c( 18,      40,      50,      60,      70,      80,  Inf)
age_labs   = c("18_39", "40_49", "50_59", "60_69", "70_79", "80_plus")

ages_cat = 
  fsubset(cohort, ed_admit_01 == 1) |>
  fselect(ca_01, age) |>
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
  fsubset(cohort, ed_admit_01 == 1) |>
  fselect(ca_01, vw) |>
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
  fsubset(cohort, ed_admit_01 == 1) |>
  fselect(ca_01, los_hosp_d) |>
  fmutate(los_cat = cut(los_hosp_d, breaks = l_breaks, labels = l_labs, right = F)) |>
  fgroup_by(ca_01, los_cat) |>
  fnobs() |>
  select(ca_01, los_cat, n = los_hosp_d) |>
  ftransform(var = "los", category = los_cat) |>
  fselect(ca_01, var, category, n) 

## prepare table component = categorical variables -----------------------------

t2_cat = 
  fsubset(cohort, ed_admit_01 == 1) |>
  select(-ends_with("id"), -ends_with("dttm"),  -age, -vw, -los_hosp_d) |>
  pivot_longer(-ca_01, names_to = "var", values_to = "val") |>
  fsubset(!is.na(val)) |>
  fmutate(n = val) |>
  fgroup_by(ca_01, var, val) |>
  fnobs() |>
  ftransform(var      = str_remove(var, "_01")) |>
  ftransform(category = tolower(str_replace_all(as.character(val), "-", "_"))) |>
  fselect(ca_01, var, category, n) 

## quality control -------------------------------------------------------------

### check for sample size mismatch in continuous table -------------------------

if (sum(t2_cont$n) != nrow(cohort[ed_admit_01 == 1])) {
  stop("ERROR: Sample size mismatch! Sum of t2_cont$n != nrow(df)")
}

## export table 2 --------------------------------------------------------------

fwrite(all_cat, here("upload_to_box", paste0("table_02_cat_",  site_lowercase, ".csv")))
fwrite(t2_cont, here("upload_to_box", paste0("table_02_cont_", site_lowercase, ".csv")))

# final cleanup ----------------------------------------------------------------

keep = c(
  "data_list",
  "site_lowercase",
  "cohort",
  "cohort_hids",
  "cohort_jids",
  "cohort_pats",
  "hid_jid_crosswalk",
  "ward_times",
  "req_vitals",
  "req_labs"
)

rm(list = setdiff(ls(), keep)); gc()

# go to 02
