# ==============================================================================
# 00_load.R
# Centralized data loading and helper functions for pooled analyses
# ==============================================================================

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

# helper functions -------------------------------------------------------------

#' Read files matching a pattern from nested folder structure
#' Structure: {site}/{analysis_type}/{filename}.csv
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
  
  # Read and combine
  file_list = lapply(all_files, function(file_path) {
    file_data = fread(file_path)
    
    # Extract site from path if not in data
    if (!"site" %in% names(file_data)) {
      # Site is the folder two levels up from the file
      path_parts = strsplit(file_path, "/")[[1]]
      site_idx   = which(path_parts %in% ALLOWED_SITES)
      if (length(site_idx) > 0) {
        file_data$site = path_parts[site_idx[1]]
      } else {
        file_data$site = "unknown"
      }
    }
    
    file_data
  })
  
  combined = rbindlist(file_list, fill = TRUE)
  
  message("  Loaded ", format(nrow(combined), big.mark = ","), " rows from ", length(file_list), " files")
  
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
    if (!"site" %in% names(file_data)) {
      file_data$site = basename(dirname(file_path))
    }
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
      cat_list[[site]] = cat_data
    }
    
    if (!is.na(cont_file)) {
      cont_data         = fread(cont_file)
      cont_data$site    = site
      cont_list[[site]] = cont_data
    }
  }
  
  cat_combined  = rbindlist(cat_list,  fill = TRUE)
  cont_combined = rbindlist(cont_list, fill = TRUE)
  
  message("  Loaded table02 data from ", length(cat_list), " sites")
  
  return(list(cat = cat_combined, cont = cont_combined))
}

#' Parse analysis variant from site column
#' Returns "main" for base sites, variant name for sensitivity analyses
parse_analysis_variant = function(site_col) {
  
  variant = fcase(
    str_detect(site_col, "se_fullcode_only"),  "fullcode_only",
    str_detect(site_col, "se_no_ed_req"),      "no_ed_req",
    str_detect(site_col, "se_win0_96h"),       "win0_96h",
    str_detect(site_col, "se_win0_120h"),      "win0_120h",
    str_detect(site_col, "se_one_enc_per_pt"), "one_enc_per_pt",
    default = "main"
  )
  
  return(variant)
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

maxscores_ca_raw     = read_grouped_files(here(), "maxscores-ca",     exclude_pattern = "liquid")
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

## diagnostics -----------------------------------------------------------------

diag_overall_raw   = read_grouped_files(here(), "overall")
diag_by_cancer_raw = read_grouped_files(here(), "by_cancer")
diag_max_raw       = read_grouped_files(here(), "max_scores")

# process and add analysis variant ---------------------------------------------

message("\n== Processing data ==")

## add analysis variant to maxscores -------------------------------------------

if (nrow(maxscores_ca_raw) > 0) {
  maxscores_ca_raw[, analysis := parse_analysis_variant(site)]
  maxscores_ca_raw[, score_name := str_remove(score_name, "_total")]
}

## add analysis variant to horizon counts --------------------------------------

if (nrow(counts_h24_raw) > 0) {
  counts_h24_raw[, analysis := parse_analysis_variant(site)]
  counts_h24_raw[, score_name := str_remove(score_name, "_total")]
}

if (nrow(counts_h12_raw) > 0) {
  counts_h12_raw[, analysis := parse_analysis_variant(site)]
  counts_h12_raw[, score_name := str_remove(score_name, "_total")]
}

# summary ----------------------------------------------------------------------

message("\n== Data loading complete ==")
message("  Sites: ", paste(ALLOWED_SITES, collapse = ", "))
message("  Maxscores rows: ", format_n(nrow(maxscores_ca_raw)))
message("  24h counts rows: ", format_n(nrow(counts_h24_raw)))
message("  Table02 cat rows: ", format_n(nrow(cat_data_raw)))
