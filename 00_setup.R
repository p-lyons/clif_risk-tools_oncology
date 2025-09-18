
# Setup script for CLIF project validating risk tools in oncology.

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# setup ------------------------------------------------------------------------

## libraries -------------------------------------------------------------------

packages_to_install =
  c(
    "data.table",
    "tidyverse",
    "tidytable",
    "collapse",
    "janitor",
    "readxl",
    "arrow",
    "rvest", 
    "readr", 
    "yaml",
    "here",
    "fst",
    "ps"
  )

packages_to_load = 
  c(
    "data.table",
    "tidytable",
    "collapse",
    "stringr",
    "arrow",
    "here",
    "ps"
  )

fn_install_if_mi = function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    try(install.packages(p, dependencies = TRUE), silent = TRUE)
  }
}

fn_load_quiet = function(p) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

invisible(lapply(packages_to_install, fn_install_if_mi))
invisible(lapply(packages_to_load,   fn_load_quiet))
options(dplyr.summarise.inform = FALSE)
rm(packages_to_install, packages_to_load, fn_install_if_mi, fn_load_quiet); gc()

## environment -----------------------------------------------------------------

### threads and ram ------------------------------------------------------------

os_type   = Sys.info()[["sysname"]]
all_cores = parallel::detectCores(logical = TRUE)
all_cores = if (is.na(all_cores)) 1L else as.integer(all_cores)

get_ram_gb = function() {
  tryCatch({
    if (os_type == "Darwin") {
      # macOS: sysctl reports total RAM in bytes
      bytes = suppressWarnings(
        as.numeric(system("sysctl -n hw.memsize", intern = TRUE))
      )
      if (length(bytes) > 0 && !is.na(bytes)) bytes / 1024^3 else NA_real_
    } else {
      # Windows & Linux: ps is usually reliable
      val = ps::ps_system_memory()[["available"]] / 1024^3
      if (is.finite(val)) val else {
        # Linux fallback: read /proc/meminfo directly
        if (file.exists("/proc/meminfo")) {
          kb = suppressWarnings(
            as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE))
          )
          if (length(kb) > 0 && !is.na(kb)) kb / 1024^2 else NA_real_
        } else {
          NA_real_
        }
      }
    }
  }, error = function(e) NA_real_)
}

avail_ram_gb = get_ram_gb()

## choose threads conservatively (00 is light; keep headroom) ------------------

reserve_cores   = 1L
gb_per_thread   = 0.50
max_by_cores    = max(1L, all_cores - reserve_cores)
max_by_memory   = if (is.finite(avail_ram_gb)) max(1L, floor(avail_ram_gb / gb_per_thread)) else max_by_cores
n_threads       = as.integer(max(1L, min(max_by_cores, max_by_memory, 8L)))
n_math_threads  = as.integer(max(1L, min(n_threads, 8L)))

## apply thread settings -------------------------------------------------------

data.table::setDTthreads(threads = n_threads)
collapse::set_collapse(nthreads  = n_threads)
options(arrow.use_threads        = TRUE)
Sys.setenv(ARROW_NUM_THREADS     = n_threads)
options(mc.cores                 = n_threads)

# concise summary
message(
  sprintf("Env OK | OS=%s | Cores=%d | Threads=%d | MathThreads=%d | Avail RAM≈%s GB",
          os_type, all_cores, n_threads, n_math_threads,
          ifelse(is.finite(avail_ram_gb), round(avail_ram_gb, 1), "NA"))
)

### site details ---------------------------------------------------------------

config        = yaml::read_yaml(here("config", "config_clif_oncrisk.yaml"))
site_details  = fread(here("config", "clif_sites.csv"))
allowed_sites = site_details$site_name
allowed_files = c("parquet", "csv", "fst")

### user enters site details ---------------------------------------------------

site_lowercase   = tolower(config$site_lowercase)
file_type        = tolower(config$file_type)  
tables_location  = config$clif_data_location # here("../_clif_data/v_2.1") 
project_location = config$project_location

### site_name must be a valid clif site in lowercase ---------------------------

if (!(site_lowercase %in% allowed_sites)) {
  stop(
    paste0(
      "Invalid '", site_lowercase,
      "'. Expected one of: ", 
      paste(allowed_sites, collapse = ", "), call. = F
    )
  )
}

### file_type must be one of c(parquet, csv) -----------------------------------

if (!(file_type %in% allowed_files)) {
  stop(
    paste0(
      "Invalid '", file_type,
      "'. Expected one of: ", 
      paste(allowed_files, collapse = ", "), call. = FALSE
    )
  )
}

### pull time zone from site details -------------------------------------------

# site_time_zone = 
#   fsubset(site_details, site_name == site_lowercase) |>
#   select(tz) |>
#   tibble::deframe()

### file locations -------------------------------------------------------------

if (!dir.exists(paste0(project_location, "/proj_tables"))) {
  dir.create(paste0(project_location, "/proj_tables"))
}
if (!dir.exists(paste0(project_location, "/proj_output"))) {
  dir.create(paste0(project_location, "/proj_output"))
}

### dates ----------------------------------------------------------------------

start_date = as.POSIXct("2016-01-01", tz = "UTC") # could be site_time_zone if we care...
end_date   = as.POSIXct("2024-12-31", tz = "UTC") # could be site_time_zone if we care...
today      = format(Sys.Date(), "%y%m%d")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# data -------------------------------------------------------------------------

## clif tables needed for project ----------------------------------------------

required_tables = 
  c(
    "patient", 
    "hospital_diagnosis",
    "hospitalization", 
    "adt", 
    "vitals", 
    "labs", 
    "medication_admin_continuous", 
    "respiratory_support",
    "code_status"
  )

## check required tables against available tables ------------------------------

clif_table_filenames = 
  list.files(
    path       = tables_location, 
    pattern    = paste0("^clif_.*\\.", file_type, "$"), 
    full.names = TRUE
  )

clif_table_basenames = 
  basename(clif_table_filenames) |>
  str_remove(paste0("\\.", file_type, "$")) |> # remove extension
  str_remove("^clif_") |>                      # remove leading 'clif_'
  str_remove("(_\\d{4}(_\\d{4})*)$")           # remove any date ranges in file name

table_file_map     = setNames(clif_table_filenames, clif_table_basenames)
missing_tables     = setdiff(required_tables, clif_table_basenames)
required_filenames = table_file_map[required_tables]

if (length(missing_tables) > 0) {
  stop(paste("Error: Missing required tables:", paste(missing_tables, collapse = ", ")))
} else {
  message("All required tables are present.")
}

rm(site_details, allowed_files, allowed_sites, missing_tables); gc()

## load tables -----------------------------------------------------------------

if (file_type == "parquet") {
  data_list = lapply(required_filenames, open_dataset)
} else if (file_type == "csv") {
  data_list = lapply(required_filenames, \(f) read_csv_arrow(f))
} else if (file_type == "fst") {
  data_list = lapply(required_filenames, \(f) {
    tmp = read.fst(f, as.data.table = TRUE)
    arrow_table(tmp)
  })
} else {
  stop("Unsupported file format")
}

names(data_list) = names(required_filenames)

## validate table contents -----------------------------------------------------

### function to validate a table -----------------------------------------------
## validate table contents -----------------------------------------------------

### function to validate a table -----------------------------------------------

validate_table = function(tbl, table_name, req_vars = NULL, req_values = list()) {
  
  if (is.null(req_vars)) req_vars = character(0)
  
  problems     = character()
  missing_vars = setdiff(req_vars, names(tbl))
  
  if (length(missing_vars)) {
    problems = c(problems, sprintf("Missing required vars: %s", paste(missing_vars, collapse = ", ")))
  }
  
  for (var in names(req_values)) {
    if (!var %in% names(tbl)) {
      problems = c(problems, sprintf("Missing '%s' needed for value checks.", var))
      next
    }
    
    ### special handling for FileSystemDataset
    if (inherits(tbl, "FileSystemDataset")) {
      present_vals = character(0)
      tryCatch({
        value_counts =
          dplyr::select(tbl, !!rlang::sym(var)) |>
          dplyr::group_by(!!rlang::sym(var)) |>
          dplyr::summarize(n = dplyr::n()) |>
          dplyr::arrange(dplyr::desc(n)) |>
          dplyr::collect()
        if (nrow(value_counts) > 0) {
          present_vals = na.omit(as.character(value_counts[[var]]))
        }
      }, error = function(e) { })
      
      if (length(present_vals) == 0) {
        tryCatch({
          sample_data =
            dplyr::select(tbl, !!rlang::sym(var)) |>
            utils::head(1000) |>
            dplyr::collect()
          if (nrow(sample_data) > 0) {
            present_vals = na.omit(unique(as.character(sample_data[[var]])))
          }
        }, error = function(e) { })
      }
    } else if (!is.data.frame(tbl) && inherits(tbl, "Table")) {
      present_vals = character(0)
      tryCatch({
        distinct_vals = tbl |>
          dplyr::select(!!rlang::sym(var)) |>
          dplyr::distinct() |>
          dplyr::collect()
        if (nrow(distinct_vals) > 0) {
          present_vals = na.omit(as.character(distinct_vals[[var]]))
        }
      }, error = function(e) { })
    } else {
      present_vals = na.omit(unique(as.character(tbl[[var]])))
    }
    
    present_vals  = tolower(trimws(as.character(present_vals)))
    expected_vals = unique(tolower(trimws(req_values[[var]])))
    missing_vals  = setdiff(expected_vals, present_vals)
    
    if (length(missing_vals)) {
      problems = c(
        problems,
        sprintf("Variable '%s' is missing expected values: %s", var, paste(missing_vals, collapse = ", "))
      )
    }
  }
  
  if (length(problems)) {
    return(sprintf("Table '%s':\n- %s", table_name, paste(problems, collapse = "\n- ")))
  }
  
  invisible(NULL)
}

### function to validate all tables --------------------------------------------

validate_all_tables = function(data_list, validation_specs) {
  all_problems = character()
  
  for (spec in validation_specs) {
    tbl_name = spec$table_name
    if (!tbl_name %in% names(data_list)) {
      all_problems = c(all_problems, sprintf("Table '%s' is missing entirely.", tbl_name))
      next
    }
    
    tbl = data_list[[tbl_name]]
    
    tbl_problems = 
      validate_table(
        tbl        = tbl,
        table_name = tbl_name,
        req_vars   = spec$req_vars,
        req_values = spec$req_values
      )
    
    if (!is.null(tbl_problems)) all_problems = c(all_problems, tbl_problems)
  }
  
  if (length(all_problems)) {
    stop("Validation errors found:\n", paste(all_problems, collapse = "\n\n"), call. = FALSE)
  }
  
  message("✅ Validation passed: all required tables are present and contain needed values.")
}

### list the details of each table's required elements -------------------------

#### prep vitals and labs separately, as they're used more than once -----------

req_vitals = 
  c(
    "heart_rate", 
    "respiratory_rate", 
    "sbp", 
    "spo2", 
    "temp_c"
  )

req_labs = 
  c(
    "bilirubin_total",
    "bun",
    "creatinine",
    "hemoglobin",
    "lactate",
    "pco2_arterial",
    "po2_arterial",
    "ph_arterial",
    "platelet_count",
    "so2_arterial",
    "wbc"
  )

#### make lists of other tables' validation requirements -----------------------

patient_list = 
  list(
    table_name = "patient",
    req_vars   = c("patient_id", "race_category", "ethnicity_category", "sex_category"),
    req_values = list(
      sex_category       = c("Female", "Male"),
      race_category      = c("White", "Black or African American", "Asian"),
      ethnicity_category = c("Hispanic", "Non-Hispanic")
    )
  )

hosp_list = 
  list(
    table_name = "hospitalization",
    req_vars   = c(
      "patient_id", 
      "hospitalization_id", 
      "age_at_admission", 
      "admission_dttm", 
      "discharge_dttm", 
      "discharge_category"
    ),
    req_values = list(discharge_category = c("Hospice", "Expired"))
  )

adt_list = 
  list(
    table_name = "adt",
    req_vars   = c("hospitalization_id", "hospital_id", "location_category", "in_dttm", "out_dttm"),
    req_values = list(location_category = c("ed", "icu", "ward"))
  )

dx_list = list(
  table_name = "hospital_diagnosis",
  req_vars   = c("hospitalization_id", "diagnosis_code", "diagnosis_code_format"),
  req_values = list(
    diagnosis_code_format = c("ICD10CM")
  )
)

med_list = 
  list(
    table_name = "medication_admin_continuous",
    req_vars   = c("med_group", "med_category"),
    req_values = list(
      med_group    = c("vasoactives"),
      med_category = c("norepinephrine", "vasopressin", "epinephrine")
    )
  )

resp_list = 
  list(
    table_name = "respiratory_support",
    req_vars   = c("device_category"),
    req_values = list(device_category = c("IMV"))
  )

vitals_list = 
  list(
    table_name = "vitals",
    req_vars   = c("vital_category", "vital_value", "recorded_dttm"),
    req_values = list(vital_category = req_vitals)
  )

labs_list = 
  list(
    table_name = "labs",
    req_vars   = c("lab_category", "lab_value", "lab_result_dttm"),
    req_values = list(lab_category = req_labs)
  )

validation_specs = list(
  patient_list, 
  hosp_list, 
  adt_list, 
  dx_list,
  med_list, 
  resp_list, 
  vitals_list, 
  labs_list
)

## run validation functions ----------------------------------------------------

validate_all_tables(data_list, validation_specs)

rm(labs_list, vitals_list, resp_list, med_list, adt_list, dx_list, hosp_list, patient_list, validation_specs)
rm(clif_table_basenames, clif_table_filenames, required_filenames, required_tables, table_file_map)
gc()

# go to script 01
