

# Setup script for CLIF project validating risk tools in oncology.

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# setup ------------------------------------------------------------------------

## libraries -------------------------------------------------------------------

# uvr owns package management for this pipeline. Dependencies and the R version
# constraint are declared in uvr.toml at the project root; uvr.lock is the source
# of truth; the project library lives in .uvr/library/ and becomes active only
# when the pipeline is launched through uvr. This script installs nothing. It
# confirms that the working set is present in the active library and, if anything
# is missing, stops with the command that repairs it.
#
# Provision once per site, then launch, from the project root:
#
#   uvr sync --frozen         install exactly what uvr.lock names; fail if stale
#   uvr run code/run_all.R    run with the project library and pinned R active
#
# From the RStudio console the equivalents are uvr::sync(frozen = TRUE) and
# uvr::run("code/run_all.R").
#
# Every namespace the pipeline touches is verified here, before stage 01, so a
# site missing one package fails in the first seconds rather than an hour into
# the run. pkgs_attached is attached with library(). pkgs_required is verified
# present and reached with :: at its call sites. glmmTMB (03a) and geepack (03c)
# attach themselves where they are used; both are verified now regardless.

pkgs_attached =
  c(
    "data.table",
    "tidytable",
    "collapse",
    "stringr",
    "arrow",
    "here",
    "pROC",
    "ps"
  )

pkgs_required =
  c(
    "dplyr",
    "tibble",
    "tidyr",
    "rlang",
    "lubridate",
    "yaml",
    "readxl",
    "janitor",
    "comorbidity",
    "glmmTMB",
    "geepack",
    "fst"
  )

pkgs_checked = c(pkgs_attached, pkgs_required)

# Are we running inside the uvr project library? Sourcing this file from a plain
# R session would attach whatever the system library happens to hold, which is
# the exact cross-site version drift uvr exists to prevent, so stop before
# attaching anything. Backslashes are normalized first so the test behaves
# identically on Windows.

lib_paths_norm = gsub("\\\\", "/", .libPaths())
in_uvr_lib     = any(grepl(".uvr/library", lib_paths_norm, fixed = TRUE))

if (!in_uvr_lib) {
  stop(
    paste0(
      "The active R library is not the uvr project library (.uvr/library).\n",
      "Launch the pipeline through uvr so the locked package set and the pinned ",
      "R version are active:\n",
      "  terminal: uvr run code/run_all.R\n",
      "  RStudio:  uvr::run(\"code/run_all.R\")"
    ),
    call. = FALSE
  )
}

# Does the running R match the pinned version? The library guard above is not
# sufficient on its own. uvr 0.4.x links the project library into a plain
# `Rscript` session, so the packages are correct even when the pipeline is
# launched outside `uvr run` — but the R VERSION is whatever that shell resolves.
# Pre-built binaries are compiled against a specific R minor series, so a site
# whose system R is 4.5.x could link 4.6.1 binaries and run them under 4.5.x
# without either guard above noticing. That is the same binary-compatibility
# failure the pin exists to prevent, one layer down.
#
# .r-version is the authority. It is absent only in a checkout that predates the
# pin, which reads as no constraint rather than as a failure, so an older clone
# does not stop at startup.

r_pin_path = here::here(".r-version")
r_pinned   = if (file.exists(r_pin_path)) trimws(readLines(r_pin_path, warn = FALSE)[1L]) else NA_character_
r_running  = paste(R.version$major, R.version$minor, sep = ".")

if (!is.na(r_pinned) && nzchar(r_pinned) && !identical(r_pinned, r_running)) {
  stop(
    paste0(
      "R version mismatch. This pipeline is pinned to R ", r_pinned,
      " in .r-version, but the running session is R ", r_running, ".\n",
      "Package binaries are built against a specific R minor version, so results ",
      "from a mismatched session are not comparable across sites.\n",
      "Install the pinned version and launch through uvr:\n",
      "  uvr r install ", r_pinned, "\n",
      "  uvr run code/run_all.R"
    ),
    call. = FALSE
  )
}

pkgs_missing =
  pkgs_checked[
    !vapply(pkgs_checked, \(p) nzchar(system.file(package = p)), logical(1))
  ]

if (length(pkgs_missing) > 0) {
  stop(
    paste0(
      "Not found in the uvr project library: ",
      paste(pkgs_missing, collapse = ", "), ".\n",
      "The library is out of sync with uvr.lock. Run `uvr sync --frozen` and ",
      "launch again with `uvr run code/run_all.R`."
    ),
    call. = FALSE
  )
}

invisible(
  lapply(pkgs_attached, \(p) suppressPackageStartupMessages(library(p, character.only = TRUE)))
)
options(dplyr.summarise.inform = FALSE)

# Resolved versions for every checked namespace, not only the attached ones. A
# version-driven difference between two sites is otherwise invisible after the
# fact. This table is written to BOX_DIR further down, once the artifact
# directories exist, because the stage keep-lists delete it at the first stage
# boundary.

session_packages =
  data.table::data.table(
    item  = pkgs_checked,
    value = vapply(pkgs_checked, \(p) as.character(utils::packageVersion(p)), character(1)),
    type  = "package"
  )

session_environment =
  data.table::data.table(
    item = c(
      "r_version",
      "r_pinned",
      "platform",
      "os",
      "locale",
      "timezone",
      "library_path"
    ),
    value = c(
      r_running,
      r_pinned,
      R.version$platform,
      paste(Sys.info()[["sysname"]], Sys.info()[["release"]]),
      Sys.getlocale("LC_ALL"),
      Sys.timezone(),
      .libPaths()[1L]
    ),
    type = "environment"
  )

message(
  sprintf("Namespaces verified: %d (%d attached, %d reached with ::) | R %s%s",
          length(pkgs_checked), length(pkgs_attached), length(pkgs_required),
          r_running,
          if (is.na(r_pinned)) " (no .r-version pin)" else " (matches .r-version)")
)

rm(pkgs_attached, pkgs_required, pkgs_checked, pkgs_missing, lib_paths_norm,
   in_uvr_lib, r_pin_path, r_pinned, r_running); gc()

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
      mem_info = ps::ps_system_memory()
      key      = intersect(c("available", "avail"), names(mem_info))
      val      = if (length(key) > 0) mem_info[[key]] / 1024^3 else NA_real_
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
allowed_sites = tolower(trimws(site_details$site_name))
allowed_files = c("parquet", "csv", "fst")

### user enters site details ---------------------------------------------------

# Both keys are trimmed as well as lowercased. A trailing space in the YAML
# value would otherwise fail the validators below with a message that reads
# identical to a correct entry, which is a hard error for a site to diagnose.

site_lowercase   = tolower(trimws(config$site_lowercase))
file_type        = tolower(trimws(config$file_type))
tables_location  = config$clif_data_location # here("../_clif_data/v_2.1")
project_location = config$project_location
allow_sparse_o2  = config$allow_sparse_o2

### artifact destination -------------------------------------------------------
# Single constant governing where every site-level artifact is written, by any
# script in the pipeline. Round one wrote to upload_to_box/; round two writes the
# complete set to upload_to_box_v2/ so the submitted artifacts stay frozen in
# place. Changing rounds means changing this one line.

BOX_DIR = "upload_to_box_v2"

### site_name must be a valid clif site in lowercase ---------------------------

if (!(site_lowercase %in% allowed_sites)) {
  stop(
    paste0(
      "Invalid '", site_lowercase,
      "'. Expected one of: ",
      paste(allowed_sites, collapse = ", ")
    ),
    call. = FALSE
  )
}

### file_type must be one of c(parquet, csv, fst) ------------------------------

if (!(file_type %in% allowed_files)) {
  stop(
    paste0(
      "Invalid '", file_type,
      "'. Expected one of: ",
      paste(allowed_files, collapse = ", ")
    ),
    call. = FALSE
  )
}

### allow_sparse_o2 must be a single TRUE or FALSE -----------------------------
# This key governs the sparse-oxygen guard in 02_scores.R. That guard stops the
# run when more than 90% of NEWS2 score rows have no supplemental-oxygen
# measurement within 6h, which means lpm_set and fio2_set are present as columns
# in respiratory_support but are effectively empty. Setting
#
#   allow_sparse_o2: true
#
# in config/config_clif_oncrisk.yaml downgrades that stop to a warning, and the
# run then finishes with NEWS2 scored without its supplemental-oxygen item. A
# site should set it only after the coordinating center has confirmed that its
# oxygen data are genuinely absent. The key is optional: a config file predating
# it reads as FALSE, so no site fails at startup for holding an older copy.
# Quoted values are rejected rather than coerced, so "false" fails loudly
# instead of evaluating as TRUE.

if (is.null(allow_sparse_o2)) {
  allow_sparse_o2 = FALSE
}

allow_sparse_o2_ok =
  is.logical(allow_sparse_o2) &&
  length(allow_sparse_o2) == 1L &&
  !is.na(allow_sparse_o2)

if (!allow_sparse_o2_ok) {
  stop(
    paste0(
      "Invalid 'allow_sparse_o2' in config_clif_oncrisk.yaml. Expected true or ",
      "false, unquoted. Delete the key to accept the default of false."
    ),
    call. = FALSE
  )
}

rm(allow_sparse_o2_ok)

message(sprintf("Sparse-oxygen override | allow_sparse_o2 = %s", allow_sparse_o2))

### pull time zone from site details -------------------------------------------

# site_time_zone =
#   fsubset(site_details, site_name == site_lowercase) |>
#   select(tz) |>
#   tibble::deframe()

### file locations -------------------------------------------------------------

# Resolved with here() rather than project_location so that directory creation
# and every subsequent write agree on one root.

if (!dir.exists(here("proj_tables"))) {
  dir.create(here("proj_tables"), recursive = TRUE)
}
if (!dir.exists(here(BOX_DIR))) {
  dir.create(here(BOX_DIR), recursive = TRUE)
}

# Round-two artifacts are organized into six analysis subdirectories. Creating
# them here means no script needs a defensive dir.create at write time, and a
# site can see the expected layout before the run produces anything.
box_subdirs = c("main", "threshold", "sensitivity", "horizon", "diagnostics", "meta")
for (sd in box_subdirs) {
  if (!dir.exists(here(BOX_DIR, sd))) {
    dir.create(here(BOX_DIR, sd), recursive = TRUE)
  }
}
rm(box_subdirs, sd)

### environment manifest -------------------------------------------------------
# Written here, immediately after the artifact directories exist, because the
# stage keep-lists delete session_packages and session_environment at the first
# stage boundary. One long-format table: environment facts first, then one row
# per verified package. When one site's results diverge from the other seven,
# this file answers the first question asked, which is whether that site ran a
# different environment.

env_manifest = data.table::rbindlist(list(session_environment, session_packages))
env_manifest[, site := site_lowercase]

fwrite(
  env_manifest,
  here(BOX_DIR, paste0("env_manifest_", site_lowercase, ".csv"))
)

message(
  sprintf("Environment manifest written | R %s | %d packages recorded",
          env_manifest$value[env_manifest$item == "r_version"],
          sum(env_manifest$type == "package"))
)

rm(env_manifest, session_environment, session_packages); gc()

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
    "code_status",
    "patient_assessments"
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
  stop(
    paste("Error: Missing required tables:", paste(missing_tables, collapse = ", ")),
    call. = FALSE
  )
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
    tmp = fst::read.fst(f, as.data.table = TRUE)
    arrow_table(tmp)
  })
} else {
  stop("Unsupported file format")
}

names(data_list) = names(required_filenames)

## validate table contents -----------------------------------------------------

### function to validate a table -----------------------------------------------
## validate table contents -----------------------------------------------------

### case-insensitive column resolver -------------------------------------------

.resolve <- function(tbl, v) {
  nm = names(tbl)
  nm[match(tolower(v), tolower(nm))]
}

### function to validate a table -----------------------------------------------

validate_table = function(tbl, table_name, req_vars = NULL, req_values = list()) {

  if (is.null(req_vars)) req_vars = character(0)

  problems     = character()
  actual_names = tolower(names(tbl))
  req_lower    = tolower(req_vars)
  missing_vars = req_vars[!req_lower %in% actual_names]

  if (length(missing_vars)) {
    problems = c(problems, sprintf("Missing required vars: %s", paste(missing_vars, collapse = ", ")))
  }

  for (var in names(req_values)) {
    resolved_var = .resolve(tbl, var)
    if (is.na(resolved_var)) {
      problems = c(problems, sprintf("Missing '%s' needed for value checks.", var))
      next
    }

    ### special handling for FileSystemDataset
    if (inherits(tbl, "FileSystemDataset")) {
      present_vals = character(0)
      tryCatch({
        value_counts =
          dplyr::select(tbl, !!rlang::sym(resolved_var)) |>
          dplyr::group_by(!!rlang::sym(resolved_var)) |>
          dplyr::summarize(n = dplyr::n()) |>
          dplyr::arrange(dplyr::desc(n)) |>
          dplyr::collect()
        if (nrow(value_counts) > 0) {
          present_vals = na.omit(as.character(value_counts[[resolved_var]]))
        }
      }, error = function(e) { })

      if (length(present_vals) == 0) {
        tryCatch({
          sample_data =
            dplyr::select(tbl, !!rlang::sym(resolved_var)) |>
            utils::head(1000) |>
            dplyr::collect()
          if (nrow(sample_data) > 0) {
            present_vals = na.omit(unique(as.character(sample_data[[resolved_var]])))
          }
        }, error = function(e) { })
      }
    } else if (!is.data.frame(tbl) && inherits(tbl, "Table")) {
      present_vals = character(0)
      tryCatch({
        distinct_vals = tbl |>
          dplyr::select(!!rlang::sym(resolved_var)) |>
          dplyr::distinct() |>
          dplyr::collect()
        if (nrow(distinct_vals) > 0) {
          present_vals = na.omit(as.character(distinct_vals[[resolved_var]]))
        }
      }, error = function(e) { })
    } else {
      present_vals = na.omit(unique(as.character(tbl[[resolved_var]])))
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
  req_vars   = c("hospitalization_id", "diagnosis_code", "diagnosis_code_format", "diagnosis_primary"),
  req_values = list(diagnosis_code_format = c("ICD10CM"))
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
    req_vars   = c("device_category", "lpm_set", "fio2_set"),
    req_values = list(device_category = c("IMV", "Room Air"))
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

pa_list =
  list(
    table_name = "patient_assessments",
    req_vars   = c("assessment_category", "numerical_value", "recorded_dttm"),
    req_values = list(assessment_category = "gcs_total")
  )

validation_specs = list(
  patient_list,
  hosp_list,
  adt_list,
  dx_list,
  med_list,
  resp_list,
  vitals_list,
  labs_list,
  pa_list
)

## run validation functions ----------------------------------------------------

validate_all_tables(data_list, validation_specs)

rm(labs_list, vitals_list, resp_list, med_list, adt_list, dx_list, hosp_list, patient_list, validation_specs)
rm(clif_table_basenames, clif_table_filenames, required_filenames, required_tables, table_file_map)
gc()

# go to script 01
