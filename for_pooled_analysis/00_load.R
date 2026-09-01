# 00_load.R
# Centralized data loading and helper functions for pooled analyses.
#
# Round-two revision. Site artifacts follow the filename grammar in
# code/03a_artifacts.R:
#   {artifact}[-{strata}][-{outcome}][-h{horizon}][-{variant}]-{site}.csv
# The token vocabularies are disjoint and the order is fixed, so this file
# parses every filename with a mirror of .parse_filename() and selects files
# by parsed tokens, never by regex prefix. This fixes the round-one bug in
# which a stem like "maxscores-ca" matched all five outcome files and pooled
# composite with warddeath.
#
# Object contract (decision 1a):
#   - Legacy object names (maxscores_ca_raw, sesp_raw, auroc_enc_raw, ...) hold
#     the COMPOSITE outcome only, so 01-06 run unchanged. Each carries an
#     outcome_key column ("composite" throughout).
#   - Parallel *_all_raw objects hold every outcome shipped for that family,
#     keyed by outcome_key, for the new per-outcome tables and analyses.
#   - Round-two files renamed the event-flag column o_primary_01 -> o_out.
#     Legacy objects rename it BACK to o_primary_01 so 03_threshold.R and
#     06_figures.R run unchanged; *_all_raw objects keep the shipped o_out.
#     Flip this rename when those scripts are next revised.
#   - Carry-forward variants (cf2/cf6/cf12) load into maxscores_ca_raw and
#     counts_h24_raw as analysis values alongside the four se_ variants
#     (decision 2a). Their filenames carry no outcome token; 02b_carryforward.R
#     defines them on the composite outcome, so outcome_key is set accordingly.
#   - Site-level AUROC rows flagged estimable == FALSE are dropped from the
#     legacy pooling objects (per the site-code contract: the coordinating
#     center excludes only the estimate) and retained in *_all_raw.

# setup ------------------------------------------------------------------------

library(data.table)
library(tidytable)
library(collapse)
library(stringr)
library(here)

# configuration ----------------------------------------------------------------

ALLOWED_SITES = c(
  "emory",
  "hopkins",
  "michigan",
  "ohsu",
  "rush",
  "ucmc",
  "umn",
  "upenn"
)

ANALYSIS_SUBDIRS = c("main", "threshold", "horizon", "sensitivity", "diagnostics")

STRATA_VOCAB  = c("ca", "liquid", "mets", "components")
OUTCOME_VOCAB = c("composite", "nohospice", "wardicu", "warddeath", "hospicedc")
VARIANT_VOCAB = c(
  "se_fullcode_only",
  "se_no_ed_req",
  "se_win0_96h",
  "se_one_enc_per_pt",
  "cf2",
  "cf6",
  "cf12",
  "boot",
  "bootenc"
)
CF_VARIANTS = c("cf2", "cf6", "cf12")

SCORE_LABS = c(
  sirs    = "SIRS",
  qsofa   = "qSOFA",
  mews    = "MEWS",
  news    = "NEWS",
  mews_sf = "MEWS-SF"
)

OUTCOME_LABS = c(
  composite = "Composite (ward-to-ICU, ward death, or hospice discharge)",
  nohospice = "Composite excluding hospice discharge",
  wardicu   = "Ward-to-ICU transfer",
  warddeath = "Ward death",
  hospicedc = "Hospice discharge"
)

ANALYSIS_LABS = c(
  main           = "Main",
  fullcode_only  = "Full-code only",
  no_ed_req      = "No ED requirement",
  win0_96h       = "0-96 hour window",
  one_enc_per_pt = "One encounter per patient",
  cf2            = "2-hour vitals carry-forward",
  cf6            = "6-hour vitals carry-forward",
  cf12           = "12-hour vitals carry-forward"
)

today = format(Sys.Date(), "%y%m%d")

# helper functions -------------------------------------------------------------

#' Format numbers with commas
format_n = function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

#' Calculate pooled SD from site-level summary statistics
calculate_pooled_sd = function(sd_vec, n_vec) {
  sqrt(sum(sd_vec^2 * (n_vec - 1)) / (sum(n_vec) - length(n_vec)))
}

#' Calculate pooled mean from site-level sums
calculate_pooled_mean = function(sum_vec, n_vec) {
  sum(sum_vec) / sum(n_vec)
}

#' Calculate pooled SD from site-level sum and sum of squares
calculate_sd_from_sums = function(sum_val, sumsq_val, n_val) {
  sqrt((sumsq_val - sum_val^2 / n_val) / (n_val - 1))
}

#' Mirror of .parse_filename() in code/03a_artifacts.R.
#' Any middle token that is not a stratum, outcome, or horizon is the variant.
parse_artifact_filename = function(fn) {

  stem = sub("\\.csv$", "", fn)
  toks = strsplit(stem, "-", fixed = TRUE)[[1]]
  n    = length(toks)

  out = list(
    artifact = toks[1],
    strata   = NA_character_,
    outcome  = NA_character_,
    horizon  = NA_integer_,
    variant  = NA_character_,
    site     = toks[n]
  )

  mid = if (n > 2) toks[2:(n - 1)] else character(0)
  for (tk in mid) {
    if (tk %in% STRATA_VOCAB) {
      out$strata = tk
    } else if (tk %in% OUTCOME_VOCAB) {
      out$outcome = tk
    } else if (grepl("^h[0-9]+$", tk)) {
      out$horizon = as.integer(sub("^h", "", tk))
    } else {
      out$variant = tk
    }
  }
  out
}

#' Map a filename variant token to the analysis label used downstream.
#' NA token (no variant in the filename) is the main analysis; se_ prefixes are
#' stripped to match the round-one labels; cf and boot tokens pass through.
variant_to_analysis = function(variant_token) {
  fifelse(
    is.na(variant_token),
    "main",
    sub("^se_", "", variant_token)
  )
}

#' Coerce columns corrupted by small-cell suppression (e.g., "<5") back to numeric
#' Replaces "<N" with a reproducible random integer in [0, N-1], then coerces to
#' integer. Columns that are already numeric are left untouched.
#'
#' @param dt     data.table (modified in place)
#' @param seed   RNG seed for reproducibility across runs
#' @return       invisible(dt), modified in place
coerce_suppressed_counts = function(dt, seed = 42L) {

  if (nrow(dt) == 0) return(invisible(dt))

  # columns that should be integer counts — covers every site-level artifact

  count_cols = c(
    "n",
    "n_total", "n_outcome", "n_no_outcome", "n_pos", "n_neg",
    "tp", "fp", "tn", "fn",
    "n_at_risk", "n_became_pos", "n_pos_in_bin",
    "n_pos_day0", "n_pos_day1",
    "n_encs", "n_encounters",
    "n_remaining_ca", "n_excluded_ca", "n_remaining_no", "n_excluded_no",
    "n_enc_ca", "n_enc_noca", "n_enc", "n_pat",
    "n_obs", "n_events", "n_ever_pos", "n_with_outcome",
    "n_missing", "n_hospitals", "n_ed_admit", "n_cancer", "n_score_rows",
    "n_hid_multi_pt",
    "site_n"
  )

  # numeric summary columns that could be affected
  numeric_cols = c(
    "age_sum", "age_sumsq", "vw_sum", "vw_sumsq",
    "los_sum", "los_sumsq",
    "h_from_admit_sum", "h_from_admit_sumsq",
    "sum_hours", "sumsq_hours",
    "sum_obs", "sumsq_obs", "sum_ward_hours",
    "sum_rate_per24h", "sumsq_rate_per24h",
    "sd_score", "mean_score",
    "auroc", "auroc_se", "ci_lower", "ci_upper",
    "pct_missing"
  )

  target_cols = intersect(c(count_cols, numeric_cols), names(dt))

  if (length(target_cols) == 0) return(invisible(dt))

  total_replaced = 0L
  set.seed(seed)

  for (col in target_cols) {

    if (!is.character(dt[[col]])) next

    # identify rows with suppression markers like "<5", "<10", etc.
    suppressed_idx = grepl("^<\\d+$", dt[[col]], perl = TRUE)
    n_suppressed   = sum(suppressed_idx)

    if (n_suppressed > 0) {

      # extract the upper bound from each "<N" value
      upper_bounds = as.integer(sub("^<", "", dt[[col]][suppressed_idx]))

      replacements = vapply(upper_bounds, function(ub) {
        sample.int(ub, size = 1L) - 1L   # 0 to (N-1)
      }, integer(1))

      set(dt, which(suppressed_idx), col, as.character(replacements))
      total_replaced = total_replaced + n_suppressed

      message("    -> '", col, "': replaced ", n_suppressed,
              " suppressed cells with random draws")
    }

    # coerce the full column: count cols to integer (via double to handle
    # "1234.0"), others straight to double
    if (col %in% count_cols) {
      set(dt, j = col, value = as.integer(round(as.double(dt[[col]]))))
    } else {
      set(dt, j = col, value = as.double(dt[[col]]))
    }
  }

  if (total_replaced > 0) {
    message("    -> ", total_replaced, " total suppressed cells replaced (seed = ", seed, ")")
  }

  invisible(dt)
}

#' Strip the "_total" suffix from score_name where present
clean_score_names = function(dt) {
  if ("score_name" %in% names(dt) && nrow(dt) > 0) {
    dt[, score_name := str_remove(score_name, "_total")]
  }
  invisible(dt)
}

# file catalog -----------------------------------------------------------------
# One pass over {site}/{subdir}/*.csv builds a catalog of parsed filenames.
# Every subsequent load selects catalog rows by parsed tokens. Files that fail
# the round-two contract (unknown subdirectory such as a stale meta/, unknown
# variant token such as a stale dxgroup file, or a filename site token that
# disagrees with the folder) are excluded with a printed reason and reported as
# validation errors, because they indicate a site that reran without deleting
# upload_to_box_v2/ first.

build_file_catalog = function(main_folder) {

  site_folders = list.dirs(main_folder, recursive = FALSE, full.names = TRUE)
  site_folders = site_folders[basename(site_folders) %in% ALLOWED_SITES]

  catalog_rows  = list()
  rejected_rows = list()

  for (site_folder in site_folders) {

    folder_site = basename(site_folder)
    subdirs     = list.dirs(site_folder, recursive = FALSE, full.names = FALSE)

    # unknown subdirectories (e.g., a stale meta/) are rejected wholesale
    unknown_subdirs = setdiff(subdirs, ANALYSIS_SUBDIRS)
    for (ud in unknown_subdirs) {
      ud_files = list.files(
        file.path(site_folder, ud),
        pattern    = "\\.csv$",
        full.names = FALSE
      )
      for (f in ud_files) {
        rejected_rows[[length(rejected_rows) + 1]] = tidytable(
          site   = folder_site,
          subdir = ud,
          file   = f,
          reason = "unknown analysis subdirectory (stale upload?)"
        )
      }
    }

    for (subdir in intersect(subdirs, ANALYSIS_SUBDIRS)) {

      files = list.files(
        file.path(site_folder, subdir),
        pattern    = "\\.csv$",
        full.names = FALSE
      )

      for (f in files) {

        p = parse_artifact_filename(f)

        reject_reason = fcase(
          !is.na(p$variant) & !(p$variant %in% VARIANT_VOCAB),
            "unknown variant token (stale dxgroup / round-one file?)",
          p$site != folder_site,
            "filename site token disagrees with site folder",
          default = NA_character_
        )

        if (!is.na(reject_reason)) {
          rejected_rows[[length(rejected_rows) + 1]] = tidytable(
            site   = folder_site,
            subdir = subdir,
            file   = f,
            reason = reject_reason
          )
          next
        }

        catalog_rows[[length(catalog_rows) + 1]] = tidytable(
          site     = folder_site,
          subdir   = subdir,
          artifact = p$artifact,
          strata   = p$strata,
          outcome  = p$outcome,
          horizon  = p$horizon,
          variant  = p$variant,
          path     = file.path(site_folder, subdir, f)
        )
      }
    }
  }

  catalog  = if (length(catalog_rows)  > 0) rbindlist(catalog_rows)  else data.table()
  rejected = if (length(rejected_rows) > 0) rbindlist(rejected_rows) else data.table()

  if (nrow(rejected) > 0) {
    message("  WARNING: ", nrow(rejected), " files excluded from the catalog:")
    for (i in seq_len(nrow(rejected))) {
      message("    - ", rejected$site[i], "/", rejected$subdir[i], "/",
              rejected$file[i], " [", rejected$reason[i], "]")
    }
  }

  list(catalog = catalog, rejected = rejected)
}

#' Load and combine all catalog rows matching the given token filters.
#' NA inside a filter vector matches files without that token (base %in%
#' semantics), e.g., variants = NA_character_ selects the main analysis and
#' horizons = NA_integer_ selects encounter-level files. NULL means no filter.
#' Attaches: site, analysis (label), outcome_key. Files that already ship an
#' outcome_key column (leadtime, crossclass) keep it. Carry-forward files ship
#' no outcome token; their outcome_key is set to "composite" per the site code.
load_from_catalog = function(catalog, artifacts, stratas = NULL, outcomes = NULL,
                             horizons = NULL, variants = NULL, label = NULL) {

  sel = catalog[artifact %in% artifacts]
  if (!is.null(stratas))  sel = sel[strata  %in% stratas]
  if (!is.null(outcomes)) sel = sel[outcome %in% outcomes]
  if (!is.null(horizons)) sel = sel[horizon %in% horizons]
  if (!is.null(variants)) sel = sel[variant %in% variants]

  if (is.null(label)) label = paste(artifacts, collapse = "/")

  if (nrow(sel) == 0) {
    warning("No files found in catalog for: ", label)
    return(data.table())
  }

  message("  Found ", nrow(sel), " files matching '", label, "'")

  file_list = lapply(seq_len(nrow(sel)), function(i) {

    file_data = fread(sel$path[i])

    file_data$site         = sel$site[i]
    file_data$analysis     = variant_to_analysis(sel$variant[i])
    file_data$.source_file = basename(sel$path[i])   # for debugging

    if (!("outcome_key" %in% names(file_data))) {
      ok = sel$outcome[i]
      if (is.na(ok) && !is.na(sel$variant[i]) && sel$variant[i] %in% CF_VARIANTS) {
        ok = "composite"
      }
      file_data$outcome_key = ok
    }

    file_data
  })

  combined = rbindlist(file_list, fill = TRUE)

  # coerce any columns corrupted by small-cell suppression (e.g., "<5")
  coerce_suppressed_counts(combined)
  clean_score_names(combined)

  message("  Loaded ", format(nrow(combined), big.mark = ","),
          " rows from ", length(file_list), " files")

  return(combined)
}

#' Composite-only legacy view of a loaded family. Renames the round-two event
#' flag o_out back to o_primary_01 so 03_threshold.R and 06_figures.R run
#' unchanged, and drops AUROC rows flagged non-estimable from pooling.
legacy_view = function(dt) {

  if (nrow(dt) == 0) return(dt)

  out = dt[outcome_key == "composite"]

  if ("o_out" %in% names(out)) {
    setnames(out, "o_out", "o_primary_01")
  }
  if ("estimable" %in% names(out)) {
    out = out[estimable == TRUE | is.na(estimable)]
  }
  out
}

#' Read files from site root folders (not in analysis subfolders)
read_site_root_files = function(main_folder, stem, exclude_pattern = NULL) {

  all_files = character(0)

  # Get site folders
  site_folders = list.dirs(main_folder, recursive = FALSE, full.names = TRUE)
  site_folders = site_folders[basename(site_folders) %in% ALLOWED_SITES]

  for (site_folder in site_folders) {
    pattern_underscore = paste0("^", gsub("-", "_", stem), ".*\\.csv$")
    pattern_hyphen     = paste0("^", gsub("_", "-", stem), ".*\\.csv$")

    files_found = list.files(
      site_folder,
      pattern    = paste0("(", pattern_underscore, "|", pattern_hyphen, ")"),
      full.names = TRUE
    )
    all_files = c(all_files, files_found)
  }

  all_files = unique(all_files)

  if (!is.null(exclude_pattern)) {
    all_files = all_files[!grepl(exclude_pattern, all_files)]
  }

  if (length(all_files) == 0) {
    warning("No files found matching pattern: ", stem)
    return(data.table())
  }

  message("  Found ", length(all_files), " files matching '", stem, "'")

  file_list = lapply(all_files, function(file_path) {
    file_data = fread(file_path)
    path_parts = strsplit(file_path, "/")[[1]]
    site_idx   = which(path_parts %in% ALLOWED_SITES)
    file_data$site = if (length(site_idx) > 0) path_parts[site_idx[1]] else "unknown"
    file_data$analysis = "main"
    file_data$.source_file = basename(file_path)
    file_data
  })

  combined = rbindlist(file_list, fill = TRUE)

  # coerce any columns corrupted by small-cell suppression (e.g., "<5")
  coerce_suppressed_counts(combined)

  message("  Loaded ", format(nrow(combined), big.mark = ","),
          " rows from ", length(file_list), " files")

  return(combined)
}

#' Read table 02 data (categorical and continuous) from site folders
read_table02_data = function(main_folder) {

  site_folders = list.dirs(main_folder, recursive = FALSE, full.names = FALSE)
  site_folders = site_folders[site_folders %in% ALLOWED_SITES]

  cat_list  = list()
  cont_list = list()

  for (site in site_folders) {

    # Try different naming conventions
    cat_patterns = c(
      file.path(main_folder, site, paste0("table_02_cat_",  site, ".csv")),
      file.path(main_folder, site, paste0("table02_cat_",   site, ".csv"))
    )

    cont_patterns = c(
      file.path(main_folder, site, paste0("table_02_cont_", site, ".csv")),
      file.path(main_folder, site, paste0("table02_cont_",  site, ".csv"))
    )

    cat_file  = cat_patterns[file.exists(cat_patterns)][1]
    cont_file = cont_patterns[file.exists(cont_patterns)][1]

    if (!is.na(cat_file)) {
      cat_data          = fread(cat_file)
      cat_data$site     = site
      cat_data$analysis = "main"
      cat_list[[site]]  = cat_data
    }

    if (!is.na(cont_file)) {
      cont_data          = fread(cont_file)
      cont_data$site     = site
      cont_data$analysis = "main"
      cont_list[[site]]  = cont_data
    }
  }

  cat_combined  = rbindlist(cat_list,  fill = TRUE)
  cont_combined = rbindlist(cont_list, fill = TRUE)

  # coerce any columns corrupted by small-cell suppression (e.g., "<5")
  coerce_suppressed_counts(cat_combined)
  coerce_suppressed_counts(cont_combined)

  message("  Loaded table02 data from ", length(cat_list), " sites")

  return(list(cat = cat_combined, cont = cont_combined))
}

# build the catalog ------------------------------------------------------------

message("\n== Building file catalog ==")

catalog_build   = build_file_catalog(here())
FILE_CATALOG    = catalog_build$catalog
REJECTED_FILES  = catalog_build$rejected

if (nrow(FILE_CATALOG) > 0) {
  catalog_summary = FILE_CATALOG[, .N, by = .(site, subdir)]
  message("  Catalog: ", nrow(FILE_CATALOG), " files across ",
          uniqueN(FILE_CATALOG$site), " sites")
  print(dcast(catalog_summary, site ~ subdir, value.var = "N", fill = 0))
} else {
  warning("File catalog is empty; no site subdirectory files were found")
}

# load all data ----------------------------------------------------------------

message("\n== Loading pooled data ==")

## table 02 (characteristics) --------------------------------------------------

table02       = read_table02_data(here())
cat_data_raw  = table02$cat
cont_data_raw = table02$cont

## site root files -------------------------------------------------------------

flow_data_raw      = read_site_root_files(here(), "figure_s01_flow")
run_report_raw     = read_site_root_files(here(), "run_report")
hospital_types_raw = read_site_root_files(here(), "hospital_types")

## admission diagnoses (R09) ---------------------------------------------------
# Round-two exports from 02_scores.R. They sit at the site root rather than in an
# analysis subdirectory, so build_file_catalog() never sees them and they need
# the root reader. Primary hospital diagnosis per encounter, tallied by ICD-10-CM
# chapter letter and by three-character code stem, stratified by cancer status.
# Site tallies of five or fewer are already suppressed at export, so pooled
# percentages are computed over the retained rows.

adm_dx_chapter_raw = read_site_root_files(here(), "admission_dx_chapter")
adm_dx_stem_raw    = read_site_root_files(here(), "admission_dx_stem")

## maxscores (encounter-level) -------------------------------------------------
# maxscores_ca_raw carries analysis = main, the four se_ variants, and the
# three carry-forward variants (decision 2a), composite outcome only.

maxscores_ca_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "maxscores",
  stratas   = "ca",
  label     = "maxscores-ca (all outcomes and variants)"
)
maxscores_ca_raw = legacy_view(maxscores_ca_all_raw)

maxscores_liquid_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "maxscores",
  stratas   = "liquid",
  label     = "maxscores-liquid"
)
maxscores_liquid_raw = legacy_view(maxscores_liquid_all_raw)

maxscores_mets_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "maxscores",
  stratas   = "mets",
  label     = "maxscores-mets"
)

## horizon counts --------------------------------------------------------------
# counts_h24_raw carries main + se + cf variants (cf exists at h24 only).

counts_h24_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "counts",
  stratas   = "ca",
  horizons  = 24L,
  variants  = c(NA_character_, "se_fullcode_only", "se_no_ed_req",
                "se_win0_96h", "se_one_enc_per_pt", CF_VARIANTS),
  label     = "counts-ca h24 (non-bootstrap)"
)
counts_h24_raw = legacy_view(counts_h24_all_raw)

counts_h12_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "counts",
  stratas   = "ca",
  horizons  = 12L,
  variants  = c(NA_character_, "se_fullcode_only", "se_no_ed_req",
                "se_win0_96h", "se_one_enc_per_pt"),
  label     = "counts-ca h12 (non-bootstrap)"
)
counts_h12_raw = legacy_view(counts_h12_all_raw)

counts_liquid_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "counts",
  stratas   = "liquid",
  horizons  = 24L,
  label     = "counts-liquid h24"
)
counts_liquid_raw = legacy_view(counts_liquid_all_raw)

counts_mets_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "counts",
  stratas   = "mets",
  horizons  = 24L,
  label     = "counts-mets h24"
)

## bootstrap counts ------------------------------------------------------------

boot_h24_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "counts",
  stratas   = "ca",
  horizons  = 24L,
  variants  = "boot",
  label     = "counts-ca h24 bootstrap"
)

boot_h12_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "counts",
  stratas   = "ca",
  horizons  = 12L,
  variants  = "boot",
  label     = "counts-ca h12 bootstrap"
)

## Encounter-weighted bootstrap (variant "bootenc"). Same file grammar as
## "boot", different estimand: one observation is drawn per resampled encounter,
## so long stays no longer dominate. "boot" gives the confidence interval for
## the reported observation-weighted fixed-horizon AUROC; "bootenc" is the
## prespecified sensitivity analysis for differential time at risk.

bootenc_h24_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "counts",
  stratas   = "ca",
  horizons  = 24L,
  variants  = "bootenc",
  label     = "counts-ca h24 encounter-weighted bootstrap"
)

bootenc_h12_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "counts",
  stratas   = "ca",
  horizons  = 12L,
  variants  = "bootenc",
  label     = "counts-ca h12 encounter-weighted bootstrap"
)

## threshold analyses ----------------------------------------------------------

ever_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "ever",
  stratas   = "ca",
  label     = "ever-ca"
)
ever_positive_raw = legacy_view(ever_all_raw)

ever_liquid_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "ever",
  stratas   = "liquid",
  label     = "ever-liquid"
)

sesp_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "sesp",
  stratas   = "ca",
  label     = "sesp-ca"
)
sesp_raw = legacy_view(sesp_all_raw)

sesp_liquid_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "sesp",
  stratas   = "liquid",
  label     = "sesp-liquid"
)

cuminc_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "cuminc",
  stratas   = "ca",
  label     = "cuminc-ca"
)
cuminc_raw = legacy_view(cuminc_all_raw)

first_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "first",
  stratas   = "ca",
  label     = "first-ca"
)
first_pos_raw = legacy_view(first_all_raw)

upset_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "upset",
  stratas   = "ca",
  label     = "upset-ca"
)
upset_raw = legacy_view(upset_all_raw)

upset_comp_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "upset",
  stratas   = "components",
  label     = "upset-components"
)
upset_comp_raw = legacy_view(upset_comp_all_raw)

leadtime_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "leadtime",
  stratas   = "ca",
  label     = "leadtime-ca"
)

## Exact per-unit lead-time medians (health system and hospital_id), each with
## its own n_events. These are deliberately NOT pooled into a single median;
## 08_supplement.R reports the range across units after applying an event floor,
## and the cumulative counts in leadtime_raw supply the encounter-weighted
## counterpart independently.

leadtime_median_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "leadtime_median",
  stratas   = "ca",
  label     = "leadtime_median-ca"
)

crossclass_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "crossclass",
  stratas   = "ca",
  label     = "crossclass-ca"
)

## site-level AUROCs -----------------------------------------------------------

message("\n  Loading site-level AUROCs...")

auroc_enc_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "auroc",
  stratas   = "ca",
  horizons  = NA_integer_,
  label     = "auroc-ca (encounter)"
)
auroc_enc_raw = legacy_view(auroc_enc_all_raw)

auroc_h24_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "auroc",
  stratas   = "ca",
  horizons  = 24L,
  label     = "auroc-ca h24"
)
auroc_h24_raw = legacy_view(auroc_h24_all_raw)

auroc_h12_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "auroc",
  stratas   = "ca",
  horizons  = 12L,
  label     = "auroc-ca h12"
)
auroc_h12_raw = legacy_view(auroc_h12_all_raw)

auroc_liquid_enc_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "auroc",
  stratas   = "liquid",
  horizons  = NA_integer_,
  label     = "auroc-liquid (encounter)"
)
auroc_liquid_enc = legacy_view(auroc_liquid_enc_all_raw)

auroc_liquid_h24_all_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "auroc",
  stratas   = "liquid",
  horizons  = 24L,
  label     = "auroc-liquid h24"
)
auroc_liquid_h24 = legacy_view(auroc_liquid_h24_all_raw)

auroc_mets_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "auroc",
  stratas   = "mets",
  horizons  = NA_integer_,
  label     = "auroc-mets (encounter)"
)

auroc_mets_h24_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "auroc",
  stratas   = "mets",
  horizons  = 24L,
  label     = "auroc-mets h24"
)

## subgroup event tallies ------------------------------------------------------

events_liquid_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "events",
  stratas   = "liquid",
  label     = "events-liquid"
)

events_mets_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "events",
  stratas   = "mets",
  label     = "events-mets"
)

## first-observation score distributions ---------------------------------------

firstscore_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "firstscore",
  stratas   = "ca",
  label     = "firstscore-ca"
)

## diagnostics -----------------------------------------------------------------

monitoring_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "monitoring",
  stratas   = "ca",
  label     = "monitoring-ca"
)

monitoring_bins_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "monitoring_bins",
  stratas   = "ca",
  label     = "monitoring_bins-ca"
)

missing_vlab_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "missing_vlab",
  stratas   = "ca",
  label     = "missing_vlab-ca"
)

news_o2_raw = load_from_catalog(
  FILE_CATALOG,
  artifacts = "news_o2_resolution",
  label     = "news_o2_resolution"
)

# ==============================================================================
# DERIVED CONSTANTS: COHORT_N, VARIANT_N, SITE_N
# ==============================================================================
# The round-one diagnostics files (overall, by_cancer, max_scores) no longer
# exist. Sources, per the site code:
#   - COHORT_N / SITE_N: the flow file's final reconciliation row ("no
#     calculable score before the outcome"), which 02_scores.R writes from the
#     reconciled ED-admit cohort and which run_report$n_ed_admit equals by
#     construction.
#   - VARIANT_N: encounter counts summed from maxscores-ca-composite by
#     analysis variant (one score as reference). The carry-forward variants
#     score the identical encounter set (02b_carryforward.R), so they inherit
#     the main-analysis count.

message("\n== Computing cohort and variant Ns ==")

FLOW_FINAL_STEP = "no calculable score"

if (nrow(flow_data_raw) > 0) {

  flow_final = flow_data_raw[step %like% FLOW_FINAL_STEP, .(
    site,
    flow_n_ca    = n_remaining_ca,
    flow_n_no    = n_remaining_no,
    flow_n_total = n_remaining_ca + n_remaining_no
  )]

  COHORT_N_DT = flow_final[, .(
    n_enc_ca   = sum(flow_n_ca),
    n_enc_noca = sum(flow_n_no),
    n_enc      = sum(flow_n_total)
  )]

  # Named vector for format_cohort() function
  COHORT_N = setNames(
    c(COHORT_N_DT$n_enc_noca, COHORT_N_DT$n_enc_ca),
    c("0", "1")
  )

  SITE_N = setNames(flow_final$flow_n_total, flow_final$site)

  message("  Main cohort: ", format_n(COHORT_N["0"]), " non-cancer, ",
          format_n(COHORT_N["1"]), " cancer encounters")
  message("  Site Ns: ", paste(names(SITE_N), "=", format_n(SITE_N), collapse = ", "))

} else {
  COHORT_N = setNames(c(NA_integer_, NA_integer_), c("0", "1"))
  SITE_N   = setNames(rep(NA_integer_, length(ALLOWED_SITES)), ALLOWED_SITES)
  warning("Could not determine COHORT_N / SITE_N: no flow data loaded")
}

if (nrow(maxscores_ca_raw) > 0) {

  VARIANT_N_DT = maxscores_ca_raw[score_name == "sirs", .(
    n_enc      = sum(n),
    n_enc_ca   = sum(n[ca_01 == 1]),
    n_enc_noca = sum(n[ca_01 == 0])
  ), by = analysis]

  setnames(VARIANT_N_DT, "analysis", "variant_clean")

  # carry-forward variants share the main cohort; assert rather than assume
  cf_check = VARIANT_N_DT[variant_clean %in% CF_VARIANTS & 
                            n_enc != VARIANT_N_DT[variant_clean == "main"]$n_enc]
  if (nrow(cf_check) > 0) {
    warning("Carry-forward variant encounter counts differ from main: ",
            paste(cf_check$variant_clean, collapse = ", "))
  }

  VARIANT_N = setNames(VARIANT_N_DT$n_enc, VARIANT_N_DT$variant_clean)

  message("  Variant Ns:")
  for (v in names(VARIANT_N)) {
    message("    ", v, ": ", format_n(VARIANT_N[v]))
  }

} else {
  VARIANT_N = c(main = sum(COHORT_N))
  warning("Could not determine VARIANT_N from maxscores")
}

# summary ----------------------------------------------------------------------

message("\n== Data loading complete ==")
message("  Sites: ", paste(ALLOWED_SITES, collapse = ", "))
message("  Maxscores rows (composite): ", format_n(nrow(maxscores_ca_raw)))
message("  Maxscores rows (all outcomes): ", format_n(nrow(maxscores_ca_all_raw)))
message("  24h counts rows: ", format_n(nrow(counts_h24_raw)))
message("  Site-level AUROCs (encounter, composite): ", format_n(nrow(auroc_enc_raw)))
message("  Site-level AUROCs (encounter, all outcomes): ", format_n(nrow(auroc_enc_all_raw)))
message("  Site-level AUROCs (24h, composite): ", format_n(nrow(auroc_h24_raw)))
message("  Site-level AUROCs (12h, all outcomes): ", format_n(nrow(auroc_h12_all_raw)))
message("  Admission diagnosis rows (chapter / stem): ",
        format_n(nrow(adm_dx_chapter_raw)), " / ", format_n(nrow(adm_dx_stem_raw)))
message("  Table02 cat rows: ", format_n(nrow(cat_data_raw)))

# Validation check -------------------------------------------------------------

message("\n== Validation ==")

validation_errors = character(0)

# Rejected files are a hard error: they indicate a stale upload (handoff item 8)
if (nrow(REJECTED_FILES) > 0) {
  validation_errors = c(
    validation_errors,
    paste0("Stale or malformed files in upload folders (sites must delete ",
           "upload_to_box_v2/ before rerunning): ",
           paste(unique(REJECTED_FILES$site), collapse = ", "))
  )
}

# Flow structure: six rows per site, ending at the reconciliation step
if (nrow(flow_data_raw) > 0) {

  flow_rows = flow_data_raw[, .(
    n_rows    = .N,
    has_final = any(step %like% FLOW_FINAL_STEP)
  ), by = site]

  bad_flow = flow_rows[n_rows != 6L | has_final == FALSE]

  if (nrow(bad_flow) > 0) {
    validation_errors = c(
      validation_errors,
      paste("Flow file does not have six rows ending at the reconciliation",
            "step for sites:", paste(bad_flow$site, collapse = ", "))
    )
  } else {
    message("  Flow files: six rows per site, reconciliation step present")
  }
}

# Flow final row vs run_report n_ed_admit (equal by construction in 02_scores.R)
if (nrow(flow_data_raw) > 0 && nrow(run_report_raw) > 0) {

  flow_vs_report = merge(
    flow_final[, .(site, flow_n_total)],
    run_report_raw[, .(site, n_ed_admit)],
    by  = "site",
    all = TRUE
  )
  flow_vs_report[, match := fifelse(flow_n_total == n_ed_admit, "OK", "MISMATCH")]

  print(flow_vs_report)

  if (any(flow_vs_report$match == "MISMATCH", na.rm = TRUE)) {
    validation_errors = c(
      validation_errors,
      paste("Flow final row disagrees with run_report n_ed_admit for sites:",
            paste(flow_vs_report[match == "MISMATCH", site], collapse = ", "))
    )
  }
}

# Check maxscores by analysis variant
if (nrow(maxscores_ca_raw) > 0) {
  n_by_analysis = maxscores_ca_raw[score_name == "sirs", .(
    total_n = sum(n),
    n_sites = uniqueN(site)
  ), by = analysis]
  message("  Maxscores by analysis (SIRS as reference):")
  print(n_by_analysis)
}

# Cross-check against flow diagram
if (nrow(flow_data_raw) > 0 && nrow(maxscores_ca_raw) > 0) {

  message("\n  Cross-checking maxscores against flow diagram...")

  # Get maxscores N for main analysis (using one score as reference)
  maxscores_main = maxscores_ca_raw[analysis == "main" & score_name == "sirs", .(
    maxscores_n_ca    = sum(n[ca_01 == 1]),
    maxscores_n_no    = sum(n[ca_01 == 0]),
    maxscores_n_total = sum(n)
  ), by = site]

  # Compare (tolerance covers small-cell suppression replacement)
  validation = merge(flow_final, maxscores_main, by = "site", all = TRUE)
  validation[, `:=`(
    diff_ca    = maxscores_n_ca - flow_n_ca,
    diff_no    = maxscores_n_no - flow_n_no,
    diff_total = maxscores_n_total - flow_n_total
  )]
  validation[, match := fifelse(
    abs(diff_ca) <= 10 & abs(diff_no) <= 10,
    "OK",
    "MISMATCH"
  )]

  print(validation[, .(site, flow_n_total, maxscores_n_total, diff_total, match)])

  if (any(validation$match == "MISMATCH", na.rm = TRUE)) {
    validation_errors = c(validation_errors,
                          paste("Flow diagram vs maxscores mismatch for sites:",
                                paste(validation[match == "MISMATCH", site], collapse = ", ")))
  }

  if (any(is.na(validation$match))) {
    validation_errors = c(validation_errors,
                          paste("Missing validation data for sites:",
                                paste(validation[is.na(match), site], collapse = ", ")))
  }
}

# Check for unexpected duplicates
if (nrow(maxscores_ca_raw) > 0) {
  dup_check = maxscores_ca_raw[analysis == "main",
                               .N,
                               by = .(site, score_name, ca_01, max_value, outcome)]
  dups = dup_check[N > 1]
  if (nrow(dups) > 0) {
    message("  WARNING: Duplicate rows detected in main analysis:")
    print(head(dups, 10))
    validation_errors = c(validation_errors,
                          paste("Duplicate rows in main analysis:", nrow(dups), "combinations affected"))
  } else {
    message("  No duplicate rows in main analysis")
  }
}

# Legacy objects must be single-outcome (the round-one pooling bug guard)
legacy_single_outcome = list(
  maxscores_ca_raw = maxscores_ca_raw,
  sesp_raw         = sesp_raw,
  counts_h24_raw   = counts_h24_raw,
  auroc_enc_raw    = auroc_enc_raw
)
for (obj_name in names(legacy_single_outcome)) {
  obj = legacy_single_outcome[[obj_name]]
  if (nrow(obj) > 0 && uniqueN(obj$outcome_key) != 1L) {
    validation_errors = c(
      validation_errors,
      paste0(obj_name, " contains more than one outcome_key: ",
             paste(unique(obj$outcome_key), collapse = ", "))
    )
  }
}
message("  Legacy objects verified single-outcome (composite)")

# Stop if any validation errors
if (length(validation_errors) > 0) {
  stop("\n\nVALIDATION FAILED:\n",
       paste(" - ", validation_errors, collapse = "\n"), "\n\n")
} else {
  message("\n  All validation checks passed")
}
