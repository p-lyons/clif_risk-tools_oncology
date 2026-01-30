# 00_load.R
# Centralized data loading and helper functions for pooled analyses

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
  "ohsu",
  "rush",
  "ucmc",
  "umn",
  "upenn"
)

SCORE_LABS = c(
  sirs    = "SIRS",
  qsofa   = "qSOFA",
  mews    = "MEWS",
  news    = "NEWS",
  mews_sf = "MEWS-SF"
)

ANALYSIS_LABS = c(
  main           = "Main",
  fullcode_only  = "Full-code only",
  no_ed_req      = "No ED requirement",
  win0_96h       = "0-96 hour window",
  one_enc_per_pt = "One encounter per patient"
)

today = format(Sys.Date(), "%y%m%d")

# helper functions -------------------------------------------------------------

#' Parse analysis variant from file path or filename
#' Checks both the folder structure and the filename for sensitivity markers
#' Returns "main" for base analyses, variant name for sensitivity analyses
parse_analysis_from_path = function(file_path) {
  
  # Check filename first (more specific)
  filename = basename(file_path)
  
  variant = fcase(
    str_detect(filename, "se_fullcode_only|se-fullcode-only"),   "fullcode_only",
    str_detect(filename, "se_no_ed_req|se-no-ed-req"),           "no_ed_req",
    str_detect(filename, "se_win0_96h|se-win0-96h"),             "win0_96h",
    str_detect(filename, "se_win0_120h|se-win0-120h"),           "win0_120h",
    str_detect(filename, "se_one_enc_per_pt|se-one-enc-per-pt"), "one_enc_per_pt",
    default = NA_character_
  )
  
  
  # If not found in filename, check folder path
  if (is.na(variant)) {
    variant = fcase(
      str_detect(file_path, "/sensitivity/"), "sensitivity",  # generic sensitivity flag
      str_detect(file_path, "/main/"),        "main",
      str_detect(file_path, "/threshold/"),   "main",
      str_detect(file_path, "/horizon/"),     "main",
      str_detect(file_path, "/meta/"),        "main",
      str_detect(file_path, "/diagnostics/"), "main",
      default = "main"
    )
  }
  
  return(variant)
}

#' Extract clean site name from file path
extract_site_from_path = function(file_path, allowed_sites) {
  path_parts = strsplit(file_path, "/")[[1]]
  site_idx = which(path_parts %in% allowed_sites)
  if (length(site_idx) > 0) {
    return(path_parts[site_idx[1]])
  }
  return("unknown")
}

#' Read files matching a pattern from nested folder structure
#' Structure: {site}/{analysis_type}/{filename}.csv
#' Extracts site and analysis variant from file path
read_grouped_files = function(main_folder, stem, exclude_pattern = NULL,
                              analysis_folders = c("main", "sensitivity", "threshold",
                                                   "horizon", "meta", "diagnostics")) {
  
  all_files = character(0)
  
  # Get site folders (emory, hopkins, etc.)
  site_folders = list.dirs(main_folder, recursive = FALSE, full.names = TRUE)
  site_folders = site_folders[basename(site_folders) %in% ALLOWED_SITES]
  
  # Search inside each site's analysis subfolders
  for (site_folder in site_folders) {
    for (analysis_folder in analysis_folders) {
      search_path = file.path(site_folder, analysis_folder)
      if (dir.exists(search_path)) {
        # Match both underscore and hyphen versions
        pattern_underscore = paste0("^", gsub("-", "_", stem), ".*\\.csv$")
        pattern_hyphen     = paste0("^", gsub("_", "-", stem), ".*\\.csv$")
        
        files_found = list.files(
          search_path,
          pattern    = paste0("(", pattern_underscore, "|", pattern_hyphen, ")"),
          full.names = TRUE
        )
        all_files = c(all_files, files_found)
      }
    }
  }
  
  # Remove duplicates
  all_files = unique(all_files)
  
  # Apply exclusion pattern
  if (!is.null(exclude_pattern)) {
    all_files = all_files[!grepl(exclude_pattern, all_files)]
  }
  
  if (length(all_files) == 0) {
    warning("No files found matching pattern: ", stem)
    return(data.table())
  }
  
  message("  Found ", length(all_files), " files matching '", stem, "'")
  
  # Read and combine, extracting metadata from path
  file_list = lapply(all_files, function(file_path) {
    file_data = fread(file_path)
    
    # Always extract site and analysis from path (more reliable than file contents)
    file_data$site = extract_site_from_path(file_path, ALLOWED_SITES)
    file_data$analysis = parse_analysis_from_path(file_path)
    file_data$.source_file = basename(file_path)  # for debugging
    
    file_data
  })
  
  combined = rbindlist(file_list, fill = TRUE)
  
  message("  Loaded ", format(nrow(combined), big.mark = ","), " rows from ", length(file_list), " files")
  
  # Report what was loaded
  analysis_summary = combined[, .N, by = analysis]
  message("  Analysis variants: ", paste(analysis_summary$analysis, collapse = ", "))
  
  return(combined)
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
    file_data$site = extract_site_from_path(file_path, ALLOWED_SITES)
    file_data$analysis = "main"
    file_data$.source_file = basename(file_path)
    file_data
  })
  
  combined = rbindlist(file_list, fill = TRUE)
  message("  Loaded ", format(nrow(combined), big.mark = ","), " rows from ", length(file_list), " files")
  
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
      cat_data         = fread(cat_file)
      cat_data$site    = site
      cat_data$analysis = "main"
      cat_list[[site]] = cat_data
    }
    
    if (!is.na(cont_file)) {
      cont_data         = fread(cont_file)
      cont_data$site    = site
      cont_data$analysis = "main"
      cont_list[[site]] = cont_data
    }
  }
  
  cat_combined  = rbindlist(cat_list,  fill = TRUE)
  cont_combined = rbindlist(cont_list, fill = TRUE)
  
  message("  Loaded table02 data from ", length(cat_list), " sites")
  
  return(list(cat = cat_combined, cont = cont_combined))
}

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

# load all data ----------------------------------------------------------------

message("\n== Loading pooled data ==")

## table 02 (characteristics) --------------------------------------------------

table02        = read_table02_data(here())
cat_data_raw   = table02$cat
cont_data_raw  = table02$cont

## flow diagram ----------------------------------------------------------------

flow_data_raw = read_site_root_files(here(), "figure_s01_flow")

## maxscores (encounter-level) -------------------------------------------------

maxscores_ca_raw     = read_grouped_files(here(), "maxscores-ca", exclude_pattern = "liquid")
maxscores_liquid_raw = read_grouped_files(here(), "maxscores-liquid")

## horizon counts (24h) --------------------------------------------------------

counts_h24_raw    = read_grouped_files(here(), "counts-ca-h24",    exclude_pattern = "boot")
counts_h12_raw    = read_grouped_files(here(), "counts-ca-h12",    exclude_pattern = "boot")
counts_liquid_raw = read_grouped_files(here(), "counts-liquid-h24", exclude_pattern = "boot")

## bootstrap counts ------------------------------------------------------------

boot_h24_raw = read_grouped_files(here(), "counts-ca-h24-boot")
boot_h12_raw = read_grouped_files(here(), "counts-ca-h12-boot")

## threshold analyses ----------------------------------------------------------

ever_positive_raw = read_grouped_files(here(), "ever-ca")
sesp_raw          = read_grouped_files(here(), "sesp-ca")
cuminc_raw        = read_grouped_files(here(), "cuminc-ca")
first_pos_raw     = read_grouped_files(here(), "first-ca")
upset_raw         = read_grouped_files(here(), "upset-ca", exclude_pattern = "components")
upset_comp_raw    = read_grouped_files(here(), "upset-components")

## meta-analysis inputs --------------------------------------------------------

coef_data_raw = read_grouped_files(here(), "coefficients", exclude_pattern = "proj_output")
score_sds_raw = read_grouped_files(here(), "score_sds",    exclude_pattern = "proj_output")

## site-level AUROCs -----------------------------------------------------------

message("\n  Loading site-level AUROCs...")

# encounter-level AUROCs (from main/ and sensitivity/ folders)
auroc_enc_raw    = read_grouped_files(here(), "auroc-ca", exclude_pattern = "h12|h24|liquid")
auroc_liquid_enc = read_grouped_files(here(), "auroc-liquid", exclude_pattern = "h12|h24")

# horizon AUROCs (from horizon/ and sensitivity/ folders)
auroc_h24_raw    = read_grouped_files(here(), "auroc-ca-h24")
auroc_h12_raw    = read_grouped_files(here(), "auroc-ca-h12")
auroc_liquid_h24 = read_grouped_files(here(), "auroc-liquid-h24")

## diagnostics -----------------------------------------------------------------

diag_overall_raw   = read_grouped_files(here(), "overall")
diag_by_cancer_raw = read_grouped_files(here(), "by_cancer")
diag_max_raw       = read_grouped_files(here(), "max_scores")

# ==============================================================================
# DERIVED CONSTANTS: COHORT_N and VARIANT_N
# ==============================================================================

message("\n== Computing cohort and variant Ns ==")

# COHORT_N: Main analysis encounter counts by cancer status
# This is the canonical source for "n=" labels in figures
if (nrow(diag_overall_raw) > 0) {
  
  # Pool across sites for main analysis
  COHORT_N_DT = diag_overall_raw[variant == "main", .(
    n_enc_ca   = sum(n_enc_ca,   na.rm = TRUE),
    n_enc_noca = sum(n_enc_noca, na.rm = TRUE),
    n_enc      = sum(n_enc,      na.rm = TRUE),
    n_pat      = sum(n_pat,      na.rm = TRUE)
  )]
  
  # Named vector for format_cohort() function
  COHORT_N = setNames(
    c(COHORT_N_DT$n_enc_noca, COHORT_N_DT$n_enc_ca),
    c("0", "1")
  )
  
  message("  Main cohort: ", format_n(COHORT_N["0"]), " non-cancer, ", 
          format_n(COHORT_N["1"]), " cancer encounters")
  
} else {
  # Fallback: calculate from flow diagram
  message("  WARNING: No diagnostics data, falling back to flow diagram")
  
  if (nrow(flow_data_raw) > 0) {
    flow_final = flow_data_raw[step %like% "outcomes too early", .(
      n_ca = sum(n_remaining_ca),
      n_no = sum(n_remaining_no)
    )]
    
    COHORT_N = setNames(
      c(flow_final$n_no, flow_final$n_ca),
      c("0", "1")
    )
  } else {
    COHORT_N = setNames(c(NA_integer_, NA_integer_), c("0", "1"))
    warning("Could not determine COHORT_N from diagnostics or flow diagram")
  }
}

# VARIANT_N: Encounter counts for each sensitivity analysis variant
# Used for SF07 and other sensitivity analysis figures
if (nrow(diag_overall_raw) > 0) {
  
  VARIANT_N_DT = diag_overall_raw[, .(
    n_enc      = sum(n_enc,      na.rm = TRUE),
    n_enc_ca   = sum(n_enc_ca,   na.rm = TRUE),
    n_enc_noca = sum(n_enc_noca, na.rm = TRUE)
  ), by = variant]
  
  # Clean variant names (remove se_ prefix for consistency)
  VARIANT_N_DT[, variant_clean := gsub("^se_", "", variant)]
  
  # Named vector for easy lookup
  VARIANT_N = setNames(VARIANT_N_DT$n_enc, VARIANT_N_DT$variant_clean)
  
  message("  Variant Ns:")
  for (v in names(VARIANT_N)) {
    message("    ", v, ": ", format_n(VARIANT_N[v]))
  }
  
} else {
  VARIANT_N = c(main = sum(COHORT_N))
  warning("Could not determine VARIANT_N from diagnostics")
}

# SITE_N: Encounter counts by site (for SF04, SF06)
if (nrow(diag_overall_raw) > 0) {
  
  SITE_N_DT = diag_overall_raw[variant == "main", .(
    n_enc      = sum(n_enc,      na.rm = TRUE),
    n_enc_ca   = sum(n_enc_ca,   na.rm = TRUE),
    n_enc_noca = sum(n_enc_noca, na.rm = TRUE)
  ), by = site]
  
  SITE_N = setNames(SITE_N_DT$n_enc, SITE_N_DT$site)
  
  message("  Site Ns: ", paste(names(SITE_N), "=", format_n(SITE_N), collapse = ", "))
  
} else {
  SITE_N = setNames(rep(NA_integer_, length(ALLOWED_SITES)), ALLOWED_SITES)
}

# clean up score names ---------------------------------------------------------

message("\n== Processing data ==")

clean_score_names = function(dt) {
  if ("score_name" %in% names(dt) && nrow(dt) > 0) {
    dt[, score_name := str_remove(score_name, "_total")]
  }
  invisible(dt)
}

clean_score_names(maxscores_ca_raw)
clean_score_names(maxscores_liquid_raw)
clean_score_names(counts_h24_raw)
clean_score_names(counts_h12_raw)
clean_score_names(counts_liquid_raw)
clean_score_names(auroc_enc_raw)
clean_score_names(auroc_liquid_enc)
clean_score_names(auroc_h24_raw)
clean_score_names(auroc_h12_raw)
clean_score_names(auroc_liquid_h24)

# summary ----------------------------------------------------------------------

message("\n== Data loading complete ==")
message("  Sites: ", paste(ALLOWED_SITES, collapse = ", "))
message("  Maxscores rows: ", format_n(nrow(maxscores_ca_raw)))
message("  24h counts rows: ", format_n(nrow(counts_h24_raw)))
message("  Site-level AUROCs (encounter): ", format_n(nrow(auroc_enc_raw)))
message("  Site-level AUROCs (24h): ", format_n(nrow(auroc_h24_raw)))
message("  Table02 cat rows: ", format_n(nrow(cat_data_raw)))

# Validation check -------------------------------------------------------------

message("\n== Validation ==")

validation_errors = character(0)

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
  
  # Get final cohort N from flow diagram (last row per site)
  flow_final = flow_data_raw[step %like% "outcomes too early", .(
    flow_n_ca = n_remaining_ca,
    flow_n_no = n_remaining_no,
    flow_n_total = n_remaining_ca + n_remaining_no
  ), by = site]
  
  # Get maxscores N for main analysis (using one score as reference)
  maxscores_main = maxscores_ca_raw[analysis == "main" & score_name == "sirs", .(
    maxscores_n_ca = sum(n[ca_01 == 1]),
    maxscores_n_no = sum(n[ca_01 == 0]),
    maxscores_n_total = sum(n)
  ), by = site]
  
  # Compare
  validation = merge(flow_final, maxscores_main, by = "site", all = TRUE)
  validation[, `:=`(
    diff_ca = maxscores_n_ca - flow_n_ca,
    diff_no = maxscores_n_no - flow_n_no,
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

# Stop if any validation errors
if (length(validation_errors) > 0) {
  stop("\n\nVALIDATION FAILED:\n
", paste(" - ", validation_errors, collapse = "\n"), "\n\n")
} else {
  message("\n  All validation checks passed")
}
