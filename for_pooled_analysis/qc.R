
# qc_federated_data.R
# Quality control script for CLIF oncology risk project
# Run this after all sites have uploaded their data folders

# setup ------------------------------------------------------------------------

library(data.table)
library(tidytable)
library(collapse)
library(stringr)
library(here)

# configuration ----------------------------------------------------------------

## expected sites (update as sites join) ---------------------------------------

expected_sites = c(
  "ohsu",
  "emory",
  "hopkins",
  "ucmc",
  "upenn",
  "umn",
  "rush"
)

## file specifications: core tables --------------------------------------------

file_specs_core = list(
  
  flow = list(
    pattern      = "figure_s01_flow_{site}.csv",
    required_cols = c("step", "n_remaining_ca", "n_excluded_ca", 
                      "n_remaining_no", "n_excluded_no"),
    location     = "root"
  ),
  
  cat = list(
    pattern      = "table_02_cat_{site}.csv",
    required_cols = c("ca_01", "var", "category", "n"),
    location     = "root"
  ),
  
  cont = list(
    pattern      = "table_02_cont_{site}.csv",
    required_cols = c("ca_01", "n", "age_sum", "age_sumsq", "vw_sum", "vw_sumsq"),
    location     = "root"
  ),
  
  cancer_codes = list(
    pattern      = "cancer_codes_primary_{site}.csv",
    required_cols = c("ca_icd10_enc", "n", "site"),
    location     = "root"
  )
)

## file specifications: analysis outputs ---------------------------------------

file_specs_analysis = list(
  
  # main analysis
  main_maxscores = list(
    pattern      = "maxscores_ca-{site}.csv",
    required_cols = c("score_name", "ca_01", "max_value", "outcome", "n", "site"),
    location     = "main"
  ),
  
  main_maxscores_liquid = list(
    pattern      = "maxscores_liquid-{site}.csv",
    required_cols = c("score_name", "liquid_01", "max_value", "outcome", "n", "site"),
    location     = "main"
  ),
  
  # threshold analyses
  threshold_ever = list(
    pattern      = "ever_ca-{site}.csv",
    required_cols = c("score_name", "ca_01", "ever_positive", "deadhospice_01", "n", "site"),
    location     = "threshold"
  ),
  
  threshold_sesp = list(
    pattern      = "sesp_ca-{site}.csv",
    required_cols = c("score_name", "ca_01", "tp", "fp", "tn", "fn", 
                      "sensitivity", "specificity", "site"),
    location     = "threshold"
  ),
  
  threshold_first = list(
    pattern      = "first_ca_outcome_times-{site}.csv",
    required_cols = c("score", "ca_01", "o_primary_01", "n", 
                      "h_from_admit_sum", "h_from_admit_sumsq"),
    location     = "threshold"
  ),
  
  threshold_upset = list(
    pattern      = "upset_ca-{site}.csv",
    required_cols = c("ca_01", "outcome", "sirs", "qsofa", "mews", "news", "mews_sf", "n"),
    location     = "threshold"
  ),
  
  threshold_upset_components = list(
    pattern      = "upset_components-{site}.csv",
    required_cols = c("ca_01", "o_primary_01", "score", "component", "n"),
    location     = "threshold"
  ),
  
  threshold_cuminc = list(
    pattern      = "cuminc_ca_outcome-{site}.csv",
    required_cols = c("score", "ca_01", "o_primary_01", "time_bin_start", 
                      "n_at_risk", "cum_inc", "site"),
    location     = "threshold"
  ),
  
  # horizon analyses
  horizon_counts_12 = list(
    pattern      = "counts_ca_h12-{site}.csv",
    required_cols = c("score_name", "ca_01", "value", "outcome", "n", "site", "h"),
    location     = "horizon"
  ),
  
  horizon_counts_24 = list(
    pattern      = "counts_ca_h24-{site}.csv",
    required_cols = c("score_name", "ca_01", "value", "outcome", "n", "site", "h"),
    location     = "horizon"
  ),
  
  horizon_counts_liquid_24 = list(
    pattern      = "counts_liquid_h24-{site}.csv",
    required_cols = c("score_name", "liquid_01", "value", "outcome", "n", "site", "h"),
    location     = "horizon"
  ),
  
  # bootstrap counts
  horizon_boot_12 = list(
    pattern      = "counts_ca_h12-boot-{site}.csv",
    required_cols = c("score_name", "ca_01", "value", "outcome", "n", "iter"),
    location     = "horizon"
  ),
  
  horizon_boot_24 = list(
    pattern      = "counts_ca_h24-boot-{site}.csv",
    required_cols = c("score_name", "ca_01", "value", "outcome", "n", "iter"),
    location     = "horizon"
  ),
  
  # diagnostics
  diag_overall = list(
    pattern      = "overall-{site}.csv",
    required_cols = c("variant", "n_rows", "n_enc", "n_pat", "site"),
    location     = "diagnostics"
  ),
  
  diag_by_cancer = list(
    pattern      = "by_cancer-{site}.csv",
    required_cols = c("variant", "ca_01", "n_enc", "evt_rate_24h", "site"),
    location     = "diagnostics"
  ),
  
  diag_max_scores = list(
    pattern      = "max_scores-{site}.csv",
    required_cols = c("variant", "ca_01", "n_enc", "event_rate", "site"),
    location     = "diagnostics"
  ),
  
  # meta-analysis inputs
  meta_coefficients = list(
    pattern      = "coefficients-{site}.csv",
    required_cols = c("score", "beta_int", "se_int", "site_n", "n_events", 
                      "converged", "site"),
    location     = "meta"
  ),
  
  meta_score_sds = list(
    pattern      = "score_sds-{site}.csv",
    required_cols = c("score", "sd_score", "mean_score", "n_encounters", "site"),
    location     = "meta"
  )
)

## sensitivity analysis variants -----------------------------------------------

sensitivity_variants = c(
  "se_no_ed_req", 
  "se_fullcode_only", 
  "se_win0_96h",        
  "se_one_enc_per_pt"
)

## clinical plausibility bounds ------------------------------------------------

bounds = list(
  
  # outcome rates (proportion)
  mortality_min    = 0.01,
  mortality_max    = 0.15,
  icu_min          = 0.03,
  icu_max          = 0.30,
  hospice_min      = 0.005,
  hospice_max      = 0.08,
  
  # cancer prevalence
  ca_prev_min      = 0.05,
  ca_prev_max      = 0.40,
  
  # demographics
  age_mean_min     = 50,
  age_mean_max     = 75,
  female_prop_min  = 0.40,
  female_prop_max  = 0.60,
  
  # van walraven
  vw_mean_min      = -5,
  vw_mean_max      = 20,
  
  # sample size
  n_min_per_site   = 1000,
  n_max_per_site   = 500000,
  
  # exclusion percentages (max reasonable % excluded at any step)
  excl_pct_max     = 50,
  
  # sensitivity/specificity bounds
  sens_min         = 0.10,
  sens_max         = 0.95,
  spec_min         = 0.30,
  spec_max         = 0.99,
  
  # time to first positive (hours)
  time_to_pos_min  = 0,
  time_to_pos_max  = 168,
  
  # meta-analysis coefficient bounds
  beta_int_min     = -2.0,
  beta_int_max     = 2.0,
  se_int_max       = 1.0
)

## score-specific valid ranges -------------------------------------------------

score_ranges = list(
  sirs_total    = c(min = 0L, max = 4L),
  qsofa_total   = c(min = 0L, max = 3L),
  mews_total    = c(min = 0L, max = 14L),
  news_total    = c(min = 0L, max = 20L),
  mews_sf_total = c(min = 0L, max = 17L)
)

# also create version without _total suffix for some outputs
score_ranges_short = list(
  sirs    = c(min = 0L, max = 4L),
  qsofa   = c(min = 0L, max = 3L),
  mews    = c(min = 0L, max = 14L),
  news    = c(min = 0L, max = 20L),
  mews_sf = c(min = 0L, max = 17L)
)

## outlier detection threshold (SDs from pooled mean) --------------------------

outlier_sd_threshold = 2.5

# helper functions -------------------------------------------------------------

## null coalescing operator ----------------------------------------------------

`%||%` = function(a, b) if (is.null(a)) b else a

## discover sites from folder structure ----------------------------------------

discover_sites = function(main_folder) {
  
  folders = list.dirs(main_folder, recursive = FALSE, full.names = FALSE)
  
  # filter to folders that look like site names
  folders = folders[!str_starts(folders, "\\.")]
  folders = folders[!folders %in% c("config", "proj_tables", "upload_to_box", 
                                    "main", "threshold", "horizon", 
                                    "sensitivity", "diagnostics", "meta")]
  folders
}

## build file path based on location -------------------------------------------

build_filepath = function(main_folder, site, file_spec) {
  
  filename = str_replace(file_spec$pattern, "\\{site\\}", site)
  
  if (file_spec$location == "root") {
    file.path(main_folder, site, filename)
  } else {
    file.path(main_folder, site, file_spec$location, filename)
  }
}

## read file for a site --------------------------------------------------------

read_site_file = function(main_folder, site, file_spec) {
  
  filepath = build_filepath(main_folder, site, file_spec)
  
  if (!file.exists(filepath)) {
    return(list(success = FALSE, error = paste("File not found:", filepath), data = NULL))
  }
  
  tryCatch({
    dt = fread(filepath)
    
    # check required columns
    missing_cols = setdiff(file_spec$required_cols, names(dt))
    if (length(missing_cols) > 0) {
      return(list(
        success = FALSE, 
        error   = paste("Missing columns:", paste(missing_cols, collapse = ", ")),
        data    = NULL
      ))
    }
    
    list(success = TRUE, error = NULL, data = dt, filepath = filepath)
    
  }, error = function(e) {
    list(success = FALSE, error = conditionMessage(e), data = NULL)
  })
}

## format numbers for display --------------------------------------------------

fmt_n   = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
fmt_pct = function(x) sprintf("%.1f%%", x * 100)
fmt_dec = function(x, d = 2) sprintf(paste0("%.", d, "f"), x)

## print section header --------------------------------------------------------

print_header = function(num, title, char = "-", width = 60) {
  cat("\n", strrep(char, width), "\n", sep = "")
  cat(num, ". ", toupper(title), "\n", sep = "")
  cat(strrep(char, width), "\n", sep = "")
}

# QC module 1: structural integrity --------------------------------------------

qc_structure = function(main_folder, sites, file_specs_core, file_specs_analysis) {
  
  results = list()
  all_specs = c(file_specs_core, file_specs_analysis)
  
  for (site in sites) {
    site_results = list(site = site, files = list())
    
    for (file_type in names(all_specs)) {
      spec   = all_specs[[file_type]]
      result = read_site_file(main_folder, site, spec)
      
      site_results$files[[file_type]] = list(
        exists   = result$success,
        error    = result$error,
        n_rows   = if (result$success) nrow(result$data) else NA,
        n_cols   = if (result$success) ncol(result$data) else NA,
        location = spec$location
      )
    }
    
    results[[site]] = site_results
  }
  
  # summarize
  summary_dt = rbindlist(lapply(names(results), function(site) {
    rbindlist(lapply(names(results[[site]]$files), function(ft) {
      f = results[[site]]$files[[ft]]
      tidytable(
        site      = site,
        file_type = ft,
        location  = f$location,
        exists    = f$exists,
        error     = f$error %||% NA_character_,
        n_rows    = f$n_rows,
        n_cols    = f$n_cols
      )
    }))
  }))
  
  list(details = results, summary = summary_dt)
}

# QC module 2: sample sizes and cancer prevalence ------------------------------

qc_sample_sizes = function(main_folder, sites, file_specs_core, bounds) {
  
  results = list()
  
  for (site in sites) {
    
    # read continuous table for N
    cont_result = read_site_file(main_folder, site, file_specs_core$cont)
    
    if (!cont_result$success) {
      results[[site]] = tidytable(site = site, error = cont_result$error)
      next
    }
    
    dt = cont_result$data
    
    n_ca    = dt[ca_01 == 1, sum(n)]
    n_noca  = dt[ca_01 == 0, sum(n)]
    n_total = n_ca + n_noca
    ca_prev = n_ca / n_total
    
    results[[site]] = tidytable(
      site         = site,
      n_total      = n_total,
      n_cancer     = n_ca,
      n_nocancer   = n_noca,
      ca_prev      = ca_prev,
      flag_n_low   = n_total < bounds$n_min_per_site,
      flag_n_high  = n_total > bounds$n_max_per_site,
      flag_ca_low  = ca_prev < bounds$ca_prev_min,
      flag_ca_high = ca_prev > bounds$ca_prev_max
    )
  }
  
  rbindlist(results, fill = TRUE)
}

# QC module 3: exclusion flow --------------------------------------------------

qc_flow = function(main_folder, sites, file_specs_core, bounds, outlier_threshold) {
  
  all_flow = list()
  
  for (site in sites) {
    result = read_site_file(main_folder, site, file_specs_core$flow)
    if (result$success) {
      dt      = copy(result$data)
      dt$site = site
      all_flow[[site]] = dt
    }
  }
  
  if (length(all_flow) == 0) {
    return(list(site_level = NULL, outliers = NULL, pooled = NULL))
  }
  
  flow_combined = rbindlist(all_flow, fill = TRUE)
  
  # assign step order based on expected sequence
  step_order = c(
    "Adult inpatient admissions during study period",
    "After excluding patients not admitted through the ED",
    "After excluding patients who were in the ICU before admission to the wards",
    "After excluding patients who were in the ICU before hitting the wards",
    "After excluding encounters with < 6h data",
    "After excluding encounters with outcomes before prediction data available",
    "After excluding encounters with outcomes too early"
  )
  
  # normalize step names and assign numeric order
  flow_combined[, step_clean := str_squish(tolower(step))]
  step_order_clean = str_squish(tolower(step_order))
  
  flow_combined[, step_num := match(step_clean, step_order_clean)]
  
  # for any unmatched steps, try partial matching
  flow_combined[is.na(step_num), step_num := {
    sapply(step_clean[is.na(step_num)], function(s) {
      matches = which(sapply(step_order_clean, function(so) grepl(substr(s, 1, 30), so) | grepl(substr(so, 1, 30), s)))
      if (length(matches) > 0) matches[1] else NA_integer_
    })
  }]
  
  # fallback: assign order by remaining count (descending)
  flow_combined[is.na(step_num), step_num := frank(-n_remaining_ca, ties.method = "dense") + 100L, by = site]
  
  setorder(flow_combined, site, step_num)
  
  # calculate exclusion percentages using PREVIOUS step's remaining
  flow_combined[, `:=`(
    prev_remaining_ca = shift(n_remaining_ca, type = "lag"),
    prev_remaining_no = shift(n_remaining_no, type = "lag")
  ), by = site]
  
  flow_combined[, `:=`(
    pct_excluded_ca = fifelse(
      !is.na(n_excluded_ca) & !is.na(prev_remaining_ca) & prev_remaining_ca > 0,
      (n_excluded_ca / prev_remaining_ca) * 100,
      NA_real_
    ),
    pct_excluded_no = fifelse(
      !is.na(n_excluded_no) & !is.na(prev_remaining_no) & prev_remaining_no > 0,
      (n_excluded_no / prev_remaining_no) * 100,
      NA_real_
    )
  )]
  
  # sanity check: cap at 100% and flag issues
  flow_combined[pct_excluded_ca > 100, pct_excluded_ca := NA_real_]
  flow_combined[pct_excluded_no > 100, pct_excluded_no := NA_real_]
  
  # pooled stats by step
  pooled = flow_combined[!is.na(pct_excluded_ca), .(
    mean_pct_ca = mean(pct_excluded_ca, na.rm = TRUE),
    sd_pct_ca   = sd(pct_excluded_ca,   na.rm = TRUE),
    mean_pct_no = mean(pct_excluded_no, na.rm = TRUE),
    sd_pct_no   = sd(pct_excluded_no,   na.rm = TRUE),
    n_sites     = .N
  ), by = .(step_num, step)]
  
  setorder(pooled, step_num)
  
  # merge and flag outliers
  flow_combined = merge(flow_combined, 
                        pooled[, .(step, mean_pct_ca, sd_pct_ca, mean_pct_no, sd_pct_no)], 
                        by = "step", all.x = TRUE)
  
  flow_combined[, `:=`(
    z_ca = (pct_excluded_ca - mean_pct_ca) / sd_pct_ca,
    z_no = (pct_excluded_no - mean_pct_no) / sd_pct_no
  )]
  
  flow_combined[, `:=`(
    flag_outlier_ca = abs(z_ca) > outlier_threshold & !is.na(z_ca),
    flag_outlier_no = abs(z_no) > outlier_threshold & !is.na(z_no),
    flag_high_excl  = pmax(pct_excluded_ca, pct_excluded_no, na.rm = TRUE) > bounds$excl_pct_max
  )]
  
  outliers = flow_combined[flag_outlier_ca == TRUE | flag_outlier_no == TRUE | flag_high_excl == TRUE]
  
  list(
    site_level = flow_combined,
    outliers   = outliers,
    pooled     = pooled
  )
}

# QC module 4: outcome rates ---------------------------------------------------

qc_outcomes = function(main_folder, sites, file_specs_core, bounds, outlier_threshold) {
  
  all_cat = list()
  
  for (site in sites) {
    result = read_site_file(main_folder, site, file_specs_core$cat)
    if (result$success) {
      dt = copy(result$data)
      if (!"site" %in% names(dt)) dt$site = site
      all_cat[[site]] = dt
    }
  }
  
  if (length(all_cat) == 0) return(NULL)
  
  cat_combined = rbindlist(all_cat, fill = TRUE)
  
  # get total N per site/cancer status
  totals = cat_combined[var == "female", .(N = sum(n)), by = .(site, ca_01)]
  
  # outcome variables
  outcome_vars = c("dead", "hospice", "icu", "wicu", "va", "imv")
  
  outcomes = cat_combined[var %in% outcome_vars & category == "1"]
  outcomes = merge(outcomes, totals, by = c("site", "ca_01"))
  outcomes[, rate := n / N]
  
  # pivot wider for easier comparison
  outcomes_wide = dcast(outcomes, site + ca_01 ~ var, value.var = "rate")
  
  # flag implausible rates
  outcomes_wide[, `:=`(
    flag_mort_low   = !is.na(dead) & dead < bounds$mortality_min,
    flag_mort_high  = !is.na(dead) & dead > bounds$mortality_max,
    flag_icu_low    = !is.na(icu) & icu < bounds$icu_min,
    flag_icu_high   = !is.na(icu) & icu > bounds$icu_max,
    flag_hosp_low   = !is.na(hospice) & hospice < bounds$hospice_min,
    flag_hosp_high  = !is.na(hospice) & hospice > bounds$hospice_max
  )]
  
  # pooled stats for outlier detection
  pooled = outcomes[, .(
    mean_rate = mean(rate, na.rm = TRUE),
    sd_rate   = sd(rate, na.rm = TRUE)
  ), by = .(var, ca_01)]
  
  outcomes = merge(outcomes, pooled, by = c("var", "ca_01"))
  outcomes[, z := (rate - mean_rate) / sd_rate]
  outcomes[, flag_outlier := abs(z) > outlier_threshold & !is.na(z)]
  
  list(
    long     = outcomes,
    wide     = outcomes_wide,
    pooled   = pooled,
    outliers = outcomes[flag_outlier == TRUE]
  )
}

# QC module 5: demographics ----------------------------------------------------

qc_demographics = function(main_folder, sites, file_specs_core, bounds, outlier_threshold) {
  
  all_cont = list()
  all_cat  = list()
  
  for (site in sites) {
    
    cont_result = read_site_file(main_folder, site, file_specs_core$cont)
    if (cont_result$success) {
      dt = copy(cont_result$data)
      if (!"site" %in% names(dt)) dt$site = site
      all_cont[[site]] = dt
    }
    
    cat_result = read_site_file(main_folder, site, file_specs_core$cat)
    if (cat_result$success) {
      dt = copy(cat_result$data)
      if (!"site" %in% names(dt)) dt$site = site
      all_cat[[site]] = dt
    }
  }
  
  results = list()
  
  # continuous: age, vw
  if (length(all_cont) > 0) {
    cont_combined = rbindlist(all_cont, fill = TRUE)
    
    cont_combined[, `:=`(
      age_mean = age_sum / n,
      age_sd   = sqrt((age_sumsq - age_sum^2 / n) / (n - 1)),
      vw_mean  = vw_sum / n,
      vw_sd    = sqrt((vw_sumsq - vw_sum^2 / n) / (n - 1))
    )]
    
    cont_summary = cont_combined[, .(site, ca_01, n, age_mean, age_sd, vw_mean, vw_sd)]
    
    cont_summary[, `:=`(
      flag_age_low  = age_mean < bounds$age_mean_min,
      flag_age_high = age_mean > bounds$age_mean_max,
      flag_vw_low   = vw_mean < bounds$vw_mean_min,
      flag_vw_high  = vw_mean > bounds$vw_mean_max
    )]
    
    results$continuous = cont_summary
  }
  
  # categorical: sex, race, ethnicity
  if (length(all_cat) > 0) {
    cat_combined = rbindlist(all_cat, fill = TRUE)
    
    totals = cat_combined[var == "female", .(N = sum(n)), by = .(site, ca_01)]
    
    # female proportion
    female = cat_combined[var == "female" & category == "1"]
    female = merge(female, totals, by = c("site", "ca_01"))
    female[, prop := n / N]
    
    female[, `:=`(
      flag_female_low  = prop < bounds$female_prop_min,
      flag_female_high = prop > bounds$female_prop_max
    )]
    
    results$female = female[, .(site, ca_01, n_female = n, N, prop_female = prop, 
                                flag_female_low, flag_female_high)]
    
    # race distribution
    race = cat_combined[var == "race_category"]
    race = merge(race, totals, by = c("site", "ca_01"))
    race[, prop := n / N]
    results$race = race[, .(site, ca_01, category, n, N, prop)]
    
    # ethnicity distribution
    ethnicity = cat_combined[var == "ethnicity_category"]
    ethnicity = merge(ethnicity, totals, by = c("site", "ca_01"))
    ethnicity[, prop := n / N]
    results$ethnicity = ethnicity[, .(site, ca_01, category, n, N, prop)]
  }
  
  results
}

# QC module 6: cross-site heterogeneity ----------------------------------------

qc_heterogeneity = function(main_folder, sites, file_specs_core) {
  
  all_cat = list()
  
  for (site in sites) {
    result = read_site_file(main_folder, site, file_specs_core$cat)
    if (result$success) {
      dt = copy(result$data)
      if (!"site" %in% names(dt)) dt$site = site
      all_cat[[site]] = dt
    }
  }
  
  if (length(all_cat) < 2) {
    return(list(message = "Need at least 2 sites for heterogeneity analysis"))
  }
  
  cat_combined = rbindlist(all_cat, fill = TRUE)
  totals = cat_combined[var == "female", .(N = sum(n)), by = .(site, ca_01)]
  
  # key proportions to check
  key_vars = c("dead", "hospice", "icu", "wicu", "va", "imv", "female")
  
  props = cat_combined[var %in% key_vars & category == "1"]
  props = merge(props, totals, by = c("site", "ca_01"))
  props[, prop := n / N]
  
  # coefficient of variation by var and cancer status
  cv_summary = props[, .(
    n_sites   = .N,
    mean_prop = mean(prop, na.rm = TRUE),
    sd_prop   = sd(prop, na.rm = TRUE),
    min_prop  = min(prop, na.rm = TRUE),
    max_prop  = max(prop, na.rm = TRUE),
    cv        = sd(prop, na.rm = TRUE) / mean(prop, na.rm = TRUE)
  ), by = .(var, ca_01)]
  
  # flag high heterogeneity (CV > 0.5)
  cv_summary[, flag_high_cv := cv > 0.5]
  
  cv_summary
}

# QC module 7: score value ranges ----------------------------------------------

qc_score_ranges = function(main_folder, sites, file_specs_analysis, score_ranges_short) {
  
  results = list()
  
  # check maxscores files
  for (site in sites) {
    
    result = read_site_file(main_folder, site, file_specs_analysis$main_maxscores)
    
    if (!result$success) next
    
    dt = result$data
    
    # check each score
    for (score_name in names(score_ranges_short)) {
      
      range_info = score_ranges_short[[score_name]]
      
      # find rows for this score (handle different naming conventions)
      score_rows = dt[grepl(score_name, score_name, ignore.case = TRUE)]
      
      if (nrow(score_rows) == 0) next
      
      min_val = min(score_rows$max_value, na.rm = TRUE)
      max_val = max(score_rows$max_value, na.rm = TRUE)
      
      out_of_range = score_rows[max_value < range_info["min"] | max_value > range_info["max"]]
      
      if (nrow(out_of_range) > 0 || min_val < range_info["min"] || max_val > range_info["max"]) {
        results[[paste(site, score_name, sep = "_")]] = tidytable(
          site           = site,
          score          = score_name,
          observed_min   = min_val,
          observed_max   = max_val,
          expected_min   = range_info["min"],
          expected_max   = range_info["max"],
          n_out_of_range = nrow(out_of_range),
          flag           = TRUE
        )
      }
    }
  }
  
  # check horizon counts files
  for (site in sites) {
    
    for (h in c("12", "24")) {
      
      spec_name = paste0("horizon_counts_", h)
      result    = read_site_file(main_folder, site, file_specs_analysis[[spec_name]])
      
      if (!result$success) next
      
      dt = result$data
      
      for (score_name in names(score_ranges_short)) {
        
        range_info = score_ranges_short[[score_name]]
        score_rows = dt[grepl(score_name, score_name, ignore.case = TRUE)]
        
        if (nrow(score_rows) == 0) next
        
        min_val = min(score_rows$value, na.rm = TRUE)
        max_val = max(score_rows$value, na.rm = TRUE)
        
        if (min_val < range_info["min"] || max_val > range_info["max"]) {
          key = paste(site, score_name, "h", h, sep = "_")
          results[[key]] = tidytable(
            site           = site,
            score          = score_name,
            horizon        = as.integer(h),
            observed_min   = min_val,
            observed_max   = max_val,
            expected_min   = range_info["min"],
            expected_max   = range_info["max"],
            flag           = TRUE
          )
        }
      }
    }
  }
  
  if (length(results) == 0) {
    return(tidytable(message = "All score values within expected ranges"))
  }
  
  rbindlist(results, fill = TRUE)
}

# QC module 8: threshold analysis outputs --------------------------------------

qc_threshold_outputs = function(main_folder, sites, file_specs_analysis, bounds) {
  
  results = list()
  
  # sensitivity/specificity checks
  sesp_results = list()
  
  for (site in sites) {
    
    result = read_site_file(main_folder, site, file_specs_analysis$threshold_sesp)
    
    if (!result$success) next
    
    dt = result$data
    
    # flag implausible values
    dt[, `:=`(
      flag_sens_low  = sensitivity < bounds$sens_min,
      flag_sens_high = sensitivity > bounds$sens_max,
      flag_spec_low  = specificity < bounds$spec_min,
      flag_spec_high = specificity > bounds$spec_max,
      flag_sum_wrong = abs((tp + fn + fp + tn) - n_total) > 1  # allow for rounding
    )]
    
    flagged = dt[flag_sens_low == TRUE | flag_sens_high == TRUE | 
                   flag_spec_low == TRUE | flag_spec_high == TRUE |
                   flag_sum_wrong == TRUE]
    
    if (nrow(flagged) > 0) {
      sesp_results[[site]] = flagged
    }
  }
  
  results$sesp = if (length(sesp_results) > 0) rbindlist(sesp_results, fill = TRUE) else NULL
  
  # cumulative incidence checks
  cuminc_results = list()
  
  for (site in sites) {
    
    result = read_site_file(main_folder, site, file_specs_analysis$threshold_cuminc)
    
    if (!result$success) next
    
    dt = result$data
    
    # cum_inc should be monotonically increasing within groups
    setorder(dt, score, ca_01, o_primary_01, time_bin_start)
    
    dt[, prev_cuminc := shift(cum_inc, type = "lag"), 
       by = .(score, ca_01, o_primary_01)]
    
    dt[, flag_nonmonotonic := !is.na(prev_cuminc) & cum_inc < prev_cuminc]
    
    # cum_inc should be between 0 and 1
    dt[, flag_bounds := cum_inc < 0 | cum_inc > 1]
    
    flagged = dt[flag_nonmonotonic == TRUE | flag_bounds == TRUE]
    
    if (nrow(flagged) > 0) {
      cuminc_results[[site]] = flagged[, .(site, score, ca_01, o_primary_01, 
                                           time_bin_start, cum_inc, 
                                           flag_nonmonotonic, flag_bounds)]
    }
  }
  
  results$cuminc = if (length(cuminc_results) > 0) rbindlist(cuminc_results, fill = TRUE) else NULL
  
  # time to first positive checks
  first_results = list()
  
  for (site in sites) {
    
    result = read_site_file(main_folder, site, file_specs_analysis$threshold_first)
    
    if (!result$success) next
    
    dt = result$data
    
    # calculate mean time
    dt[, mean_time := h_from_admit_sum / n]
    
    # flag implausible times
    dt[, flag_time := mean_time < bounds$time_to_pos_min | mean_time > bounds$time_to_pos_max]
    
    flagged = dt[flag_time == TRUE]
    
    if (nrow(flagged) > 0) {
      first_results[[site]] = flagged
    }
  }
  
  results$first = if (length(first_results) > 0) rbindlist(first_results, fill = TRUE) else NULL
  
  results
}

# QC module 9: meta-analysis inputs --------------------------------------------

qc_meta_inputs = function(main_folder, sites, file_specs_analysis, bounds) {
  
  results = list()
  
  # coefficients
  coef_results = list()
  
  for (site in sites) {
    
    result = read_site_file(main_folder, site, file_specs_analysis$meta_coefficients)
    
    if (!result$success) next
    
    dt = result$data
    
    # flag issues
    dt[, `:=`(
      flag_beta_extreme = beta_int < bounds$beta_int_min | beta_int > bounds$beta_int_max,
      flag_se_large     = se_int > bounds$se_int_max,
      flag_se_zero      = se_int <= 0,
      flag_not_converged = converged == FALSE
    )]
    
    flagged = dt[flag_beta_extreme == TRUE | flag_se_large == TRUE | 
                   flag_se_zero == TRUE | flag_not_converged == TRUE]
    
    if (nrow(flagged) > 0) {
      coef_results[[site]] = flagged
    }
    
    # also collect all for pooling comparison
    dt$site = site
    results$all_coefs = rbindlist(list(results$all_coefs, dt), fill = TRUE)
  }
  
  results$coef_flags = if (length(coef_results) > 0) rbindlist(coef_results, fill = TRUE) else NULL
  
  # check for coefficient consistency across sites
  if (!is.null(results$all_coefs) && nrow(results$all_coefs) > 0) {
    
    coef_summary = results$all_coefs[, .(
      n_sites   = .N,
      mean_beta = mean(beta_int, na.rm = TRUE),
      sd_beta   = sd(beta_int, na.rm = TRUE),
      min_beta  = min(beta_int, na.rm = TRUE),
      max_beta  = max(beta_int, na.rm = TRUE)
    ), by = score]
    
    # flag if betas have different signs across sites
    sign_check = results$all_coefs[, .(
      n_positive = sum(beta_int > 0, na.rm = TRUE),
      n_negative = sum(beta_int < 0, na.rm = TRUE)
    ), by = score]
    
    sign_check[, flag_mixed_signs := n_positive > 0 & n_negative > 0]
    
    results$coef_summary   = coef_summary
    results$coef_sign_check = sign_check
  }
  
  # score SDs
  sd_results = list()
  
  for (site in sites) {
    
    result = read_site_file(main_folder, site, file_specs_analysis$meta_score_sds)
    
    if (!result$success) next
    
    dt       = result$data
    dt$site  = site
    sd_results[[site]] = dt
  }
  
  if (length(sd_results) > 0) {
    all_sds = rbindlist(sd_results, fill = TRUE)
    
    # check for SD consistency
    sd_summary = all_sds[, .(
      n_sites = .N,
      mean_sd = mean(sd_score, na.rm = TRUE),
      cv_sd   = sd(sd_score, na.rm = TRUE) / mean(sd_score, na.rm = TRUE)
    ), by = score]
    
    # flag high variation in SDs
    sd_summary[, flag_high_cv := cv_sd > 0.3]
    
    results$sd_summary = sd_summary
  }
  
  results
}

# QC module 10: sensitivity analysis consistency -------------------------------

qc_sensitivity_analyses = function(main_folder, sites, file_specs_analysis, 
                                   sensitivity_variants) {
  
  results = list()
  
  # for each site, check that main and sensitivity variants have expected relationships
  
  for (site in sites) {
    
    # read main analysis
    main_result = read_site_file(main_folder, site, file_specs_analysis$main_maxscores)
    
    if (!main_result$success) next
    
    main_n = main_result$data[, .(n_main = sum(n)), by = .(score_name, ca_01)]
    
    # check each sensitivity variant
    for (variant in sensitivity_variants) {
      
      # construct the expected filename pattern
      variant_spec = list(
        pattern       = paste0("maxscores_ca-", variant, "-{site}.csv"),
        required_cols = file_specs_analysis$main_maxscores$required_cols,
        location      = "sensitivity"
      )
      
      var_result = read_site_file(main_folder, site, variant_spec)
      
      if (!var_result$success) next
      
      var_n = var_result$data[, .(n_var = sum(n)), by = .(score_name, ca_01)]
      
      # compare
      comparison = merge(main_n, var_n, by = c("score_name", "ca_01"), all = TRUE)
      comparison[, ratio := n_var / n_main]
      comparison[, site := site]
      comparison[, variant := variant]
      
      # flag unexpected patterns
      # - se_no_ed_req should have MORE patients than main
      # - all others should have FEWER or equal
      if (variant == "se_no_ed_req") {
        comparison[, flag := n_var < n_main * 0.95]  # allow 5% tolerance
      } else {
        comparison[, flag := n_var > n_main * 1.05]  # shouldn't be larger
      }
      
      results[[paste(site, variant, sep = "_")]] = comparison
    }
  }
  
  if (length(results) == 0) return(NULL)
  
  all_comparisons = rbindlist(results, fill = TRUE)
  
  list(
    all     = all_comparisons,
    flagged = all_comparisons[flag == TRUE]
  )
}

# QC module 11: internal consistency checks ------------------------------------

qc_internal_consistency = function(main_folder, sites, file_specs_core, file_specs_analysis) {
  
  results = list()
  
  for (site in sites) {
    
    # get N from different sources and compare
    
    # source 1: table_02_cont
    cont_result = read_site_file(main_folder, site, file_specs_core$cont)
    n_cont = if (cont_result$success) cont_result$data[, sum(n)] else NA
    
    # source 2: flow diagram final step
    flow_result = read_site_file(main_folder, site, file_specs_core$flow)
    if (flow_result$success) {
      flow_dt = flow_result$data
      final_step = flow_dt[.N]
      n_flow = final_step$n_remaining_ca + final_step$n_remaining_no
    } else {
      n_flow = NA
    }
    
    # source 3: diagnostics overall (main variant)
    diag_result = read_site_file(main_folder, site, file_specs_analysis$diag_overall)
    if (diag_result$success) {
      n_diag = diag_result$data[variant == "main", n_enc]
    } else {
      n_diag = NA
    }
    
    # compare
    n_values = c(cont = n_cont, flow = n_flow, diag = n_diag)
    n_values = n_values[!is.na(n_values)]
    
    if (length(n_values) >= 2) {
      max_diff_pct = (max(n_values) - min(n_values)) / mean(n_values) * 100
      
      results[[site]] = tidytable(
        site         = site,
        n_cont       = n_cont,
        n_flow       = n_flow,
        n_diag       = n_diag,
        max_diff_pct = max_diff_pct,
        flag         = max_diff_pct > 5  # flag if >5% discrepancy
      )
    }
  }
  
  if (length(results) == 0) return(NULL)
  
  rbindlist(results, fill = TRUE)
}

# main QC runner ---------------------------------------------------------------

run_qc = function(main_folder              = here(), 
                  expected_sites           = NULL,
                  core_specs               = NULL,
                  analysis_specs           = NULL,
                  plausibility_bounds      = NULL,
                  score_ranges             = NULL,
                  sens_variants            = NULL,
                  outlier_threshold        = 2.5,
                  check_analysis_outputs   = TRUE) {
  
  # use global defaults if not provided
  if (is.null(core_specs))           core_specs           = file_specs_core
  if (is.null(analysis_specs))       analysis_specs       = file_specs_analysis
  if (is.null(plausibility_bounds))  plausibility_bounds  = bounds
  if (is.null(score_ranges))         score_ranges         = score_ranges_short
  if (is.null(sens_variants))        sens_variants        = sensitivity_variants
  
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("FEDERATED DATA QUALITY CONTROL REPORT\n")
  cat(strrep("=", 70), "\n", sep = "")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Main folder:", main_folder, "\n\n")
  
  # discover or use expected sites
  discovered_sites = discover_sites(main_folder)
  
  if (is.null(expected_sites) || length(expected_sites) == 0) {
    sites = discovered_sites
    cat("Discovered", length(sites), "site folders:", paste(sites, collapse = ", "), "\n\n")
  } else {
    missing_sites = setdiff(expected_sites, discovered_sites)
    extra_sites   = setdiff(discovered_sites, expected_sites)
    sites         = intersect(expected_sites, discovered_sites)
    
    if (length(missing_sites) > 0) {
      cat("⚠️  MISSING EXPECTED SITES:", paste(missing_sites, collapse = ", "), "\n")
    }
    if (length(extra_sites) > 0) {
      cat("ℹ️  Extra sites found:", paste(extra_sites, collapse = ", "), "\n")
    }
    cat("Processing", length(sites), "sites:", paste(sites, collapse = ", "), "\n\n")
  }
  
  if (length(sites) == 0) {
    stop("No sites to process!")
  }
  
  results  = list()
  n_issues = 0
  
  # -------------------------------------------------------------------------
  print_header(1, "structural integrity")
  
  results$structure = qc_structure(
    main_folder, 
    sites, 
    core_specs, 
    if (check_analysis_outputs) analysis_specs else list()
  )
  
  missing_core = results$structure$summary[location == "root" & exists == FALSE]
  if (nrow(missing_core) > 0) {
    cat("⚠️  Missing CORE files:\n")
    print(missing_core[, .(site, file_type, error)])
    n_issues = n_issues + nrow(missing_core)
  } else {
    cat("✓ All core files present\n")
  }
  
  if (check_analysis_outputs) {
    missing_analysis = results$structure$summary[location != "root" & exists == FALSE]
    if (nrow(missing_analysis) > 0) {
      cat("\n⚠️  Missing ANALYSIS files:\n")
      print(missing_analysis[, .(site, file_type, location, error)][1:min(20, .N)])
      if (nrow(missing_analysis) > 20) cat("... and", nrow(missing_analysis) - 20, "more\n")
      n_issues = n_issues + nrow(missing_analysis)
    } else {
      cat("✓ All analysis output files present\n")
    }
  }
  
  # -------------------------------------------------------------------------
  print_header(2, "sample sizes & cancer prevalence")
  
  results$sample_sizes = qc_sample_sizes(main_folder, sites, core_specs, plausibility_bounds)
  
  print(results$sample_sizes[, .(
    site, 
    n_total  = fmt_n(n_total), 
    n_cancer = fmt_n(n_cancer),
    ca_prev  = fmt_pct(ca_prev)
  )])
  
  flagged_n = results$sample_sizes[flag_n_low == TRUE | flag_n_high == TRUE | 
                                     flag_ca_low == TRUE | flag_ca_high == TRUE]
  if (nrow(flagged_n) > 0) {
    cat("\n⚠️  Flagged sample size/prevalence:\n")
    print(flagged_n[, .(site, n_total, ca_prev, flag_n_low, flag_n_high, 
                        flag_ca_low, flag_ca_high)])
    n_issues = n_issues + nrow(flagged_n)
  }
  
  # -------------------------------------------------------------------------
  print_header(3, "exclusion flow")
  
  results$flow = qc_flow(main_folder, sites, core_specs, plausibility_bounds, outlier_threshold)
  
  if (!is.null(results$flow$pooled)) {
    cat("Pooled exclusion % by step:\n")
    print(results$flow$pooled[, .(
      step    = str_trunc(step, 50),
      mean_ca = fmt_dec(mean_pct_ca, 1),
      sd_ca   = fmt_dec(sd_pct_ca, 1),
      mean_no = fmt_dec(mean_pct_no, 1),
      sd_no   = fmt_dec(sd_pct_no, 1)
    )])
  }
  
  if (!is.null(results$flow$outliers) && nrow(results$flow$outliers) > 0) {
    cat("\n⚠️  Outlier exclusion %:\n")
    print(results$flow$outliers[, .(
      site, 
      step   = str_trunc(step, 30),
      pct_ca = fmt_dec(pct_excluded_ca, 1),
      pct_no = fmt_dec(pct_excluded_no, 1)
    )])
    n_issues = n_issues + nrow(results$flow$outliers)
  } else {
    cat("\n✓ No outlier exclusion percentages\n")
  }
  
  # -------------------------------------------------------------------------
  print_header(4, "outcome rates")
  
  results$outcomes = qc_outcomes(main_folder, sites, core_specs, plausibility_bounds, outlier_threshold)
  
  if (!is.null(results$outcomes$pooled)) {
    cat("Pooled outcome rates:\n")
    print(results$outcomes$pooled[, .(
      var, ca_01,
      mean = fmt_pct(mean_rate),
      sd   = fmt_pct(sd_rate)
    )])
  }
  
  if (!is.null(results$outcomes$outliers) && nrow(results$outcomes$outliers) > 0) {
    cat("\n⚠️  Outlier outcome rates:\n")
    print(results$outcomes$outliers[, .(site, ca_01, var, rate = fmt_pct(rate))])
    n_issues = n_issues + nrow(results$outcomes$outliers)
  } else {
    cat("\n✓ No outlier outcome rates\n")
  }
  
  # -------------------------------------------------------------------------
  print_header(5, "demographics")
  
  results$demographics = qc_demographics(main_folder, sites, core_specs, 
                                         plausibility_bounds, outlier_threshold)
  
  if (!is.null(results$demographics$continuous)) {
    cat("Age and Van Walraven:\n")
    print(results$demographics$continuous[, .(
      site, ca_01,
      n    = fmt_n(n),
      age  = paste0(fmt_dec(age_mean, 1), " (", fmt_dec(age_sd, 1), ")"),
      vw   = paste0(fmt_dec(vw_mean, 1), " (", fmt_dec(vw_sd, 1), ")")
    )])
    
    flagged_demo = results$demographics$continuous[
      flag_age_low == TRUE | flag_age_high == TRUE | 
        flag_vw_low == TRUE | flag_vw_high == TRUE
    ]
    if (nrow(flagged_demo) > 0) {
      cat("\n⚠️  Flagged demographics:\n")
      print(flagged_demo[, .(site, ca_01, age_mean, vw_mean)])
      n_issues = n_issues + nrow(flagged_demo)
    }
  }
  
  # -------------------------------------------------------------------------
  print_header(6, "cross-site heterogeneity")
  
  results$heterogeneity = qc_heterogeneity(main_folder, sites, core_specs)
  
  if (is.data.table(results$heterogeneity)) {
    cat("Coefficient of variation for key proportions:\n")
    print(results$heterogeneity[, .(
      var, ca_01,
      mean  = fmt_pct(mean_prop),
      range = paste0(fmt_pct(min_prop), "-", fmt_pct(max_prop)),
      cv    = fmt_dec(cv, 2),
      flag  = ifelse(flag_high_cv, "⚠️", "✓")
    )])
    
    high_cv  = results$heterogeneity[flag_high_cv == TRUE]
    n_issues = n_issues + nrow(high_cv)
  }
  
  # -------------------------------------------------------------------------
  print_header(7, "score value ranges")
  
  results$score_ranges = qc_score_ranges(main_folder, sites, analysis_specs, score_ranges)
  
  if ("message" %in% names(results$score_ranges)) {
    cat("✓", results$score_ranges$message, "\n")
  } else if (nrow(results$score_ranges) > 0) {
    cat("⚠️  Score values outside expected ranges:\n")
    print(results$score_ranges)
    n_issues = n_issues + nrow(results$score_ranges)
  }
  
  if (check_analysis_outputs) {
    
    # -----------------------------------------------------------------------
    print_header(8, "threshold analysis outputs")
    
    results$threshold = qc_threshold_outputs(main_folder, sites, analysis_specs, plausibility_bounds)
    
    if (!is.null(results$threshold$sesp) && nrow(results$threshold$sesp) > 0) {
      cat("⚠️  Sensitivity/specificity flags:\n")
      print(results$threshold$sesp[, .(site, score_name, ca_01, sensitivity, specificity)])
      n_issues = n_issues + nrow(results$threshold$sesp)
    } else {
      cat("✓ Sensitivity/specificity values plausible\n")
    }
    
    if (!is.null(results$threshold$cuminc) && nrow(results$threshold$cuminc) > 0) {
      cat("\n⚠️  Cumulative incidence flags:\n")
      print(results$threshold$cuminc)
      n_issues = n_issues + nrow(results$threshold$cuminc)
    } else {
      cat("✓ Cumulative incidence curves monotonic\n")
    }
    
    # -----------------------------------------------------------------------
    print_header(9, "meta-analysis inputs")
    
    results$meta = qc_meta_inputs(main_folder, sites, analysis_specs, plausibility_bounds)
    
    if (!is.null(results$meta$coef_flags) && nrow(results$meta$coef_flags) > 0) {
      cat("⚠️  Coefficient flags:\n")
      print(results$meta$coef_flags[, .(site, score, beta_int, se_int, converged)])
      n_issues = n_issues + nrow(results$meta$coef_flags)
    } else {
      cat("✓ All regression coefficients plausible\n")
    }
    
    if (!is.null(results$meta$coef_sign_check)) {
      mixed = results$meta$coef_sign_check[flag_mixed_signs == TRUE]
      if (nrow(mixed) > 0) {
        cat("\n⚠️  Mixed coefficient signs across sites:\n")
        print(mixed)
        n_issues = n_issues + nrow(mixed)
      }
    }
    
    if (!is.null(results$meta$sd_summary)) {
      cat("\nScore SD consistency:\n")
      print(results$meta$sd_summary[, .(
        score, 
        mean_sd = fmt_dec(mean_sd, 2),
        cv      = fmt_dec(cv_sd, 2),
        flag    = ifelse(flag_high_cv, "⚠️", "✓")
      )])
    }
    
    # -----------------------------------------------------------------------
    print_header(10, "sensitivity analysis consistency")
    
    results$sensitivity = qc_sensitivity_analyses(main_folder, sites, 
                                                  analysis_specs, 
                                                  sens_variants)
    
    if (!is.null(results$sensitivity$flagged) && nrow(results$sensitivity$flagged) > 0) {
      cat("⚠️  Unexpected sensitivity analysis sample sizes:\n")
      print(results$sensitivity$flagged[, .(site, variant, score_name, ca_01, 
                                            n_main, n_var, ratio)])
      n_issues = n_issues + nrow(results$sensitivity$flagged)
    } else {
      cat("✓ Sensitivity analysis sample sizes consistent\n")
    }
    
    # -----------------------------------------------------------------------
    print_header(11, "internal consistency")
    
    results$consistency = qc_internal_consistency(main_folder, sites, 
                                                  core_specs, analysis_specs)
    
    if (!is.null(results$consistency)) {
      flagged_cons = results$consistency[flag == TRUE]
      if (nrow(flagged_cons) > 0) {
        cat("⚠️  Sample size discrepancies across outputs:\n")
        print(flagged_cons)
        n_issues = n_issues + nrow(flagged_cons)
      } else {
        cat("✓ Sample sizes consistent across outputs\n")
      }
    }
  }
  
  # -------------------------------------------------------------------------
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("SUMMARY\n")
  cat(strrep("=", 70), "\n", sep = "")
  
  cat("Sites processed:", length(sites), "\n")
  cat("Total flags:", n_issues, "\n\n")
  
  if (n_issues == 0) {
    cat("✓ No major issues detected\n")
  } else {
    cat("⚠️  Review flagged items before pooling\n")
  }
  
  cat("\n")
  
  invisible(results)
}
  

# run the QC -------------------------------------------------------------------

qc_results = run_qc(
   main_folder           = here(),
   expected_sites        = expected_sites,
   check_analysis_outputs = TRUE
)

# export flagged items to CSV for review
 export_flags = function(results, output_dir = here("qc_output")) {
   if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#   # flow outliers
   if (!is.null(results$flow$outliers) && nrow(results$flow$outliers) > 0) {
     fwrite(results$flow$outliers, file.path(output_dir, "flag_flow_outliers.csv"))
   }
#
#   # outcome outliers
   if (!is.null(results$outcomes$outliers) && nrow(results$outcomes$outliers) > 0) {
     fwrite(results$outcomes$outliers, file.path(output_dir, "flag_outcome_outliers.csv"))
   }
#
#   # score range violations
   if (is.data.table(results$score_ranges) && nrow(results$score_ranges) > 0) {
     fwrite(results$score_ranges, file.path(output_dir, "flag_score_ranges.csv"))
   }
#
#   # meta coefficient flags
   if (!is.null(results$meta$coef_flags) && nrow(results$meta$coef_flags) > 0) {
     fwrite(results$meta$coef_flags, file.path(output_dir, "flag_meta_coefficients.csv"))
   }
#
   message("Flags exported to: ", output_dir)
 }
 
 #################
 
 # qc_report_html.R
 # HTML report generator for federated QC results
 # Requires: ggplot2, htmltools, base64enc
 
 # additional packages for reporting ---------------------------------------------
 
 library(ggplot2)
 library(htmltools)
 library(base64enc)
 
 # theme for plots --------------------------------------------------------------
 
 theme_qc = function(base_size = 11) {
   
   theme_bw(base_size = base_size) +
     theme(
       plot.title        = element_text(face = "bold", size = rel(1.1)),
       plot.subtitle     = element_text(color = "gray40"),
       panel.grid.minor  = element_blank(),
       panel.grid.major  = element_line(color = "gray90"),
       axis.title        = element_text(face = "bold"),
       legend.position   = "bottom",
       strip.text        = element_text(face = "bold"),
       strip.background  = element_rect(fill = "gray95", color = "gray70"),
       panel.border      = element_rect(color = "gray70", fill = NA, linewidth = 0.5),
       panel.spacing     = unit(0.5, "lines")
     )
 }
 
 
 # color palettes ---------------------------------------------------------------
 
 pal_cancer   = c("0" = "#4575b4", "1" = "#d73027")
 pal_outcome  = c("dead" = "#d73027", "hospice" = "#fc8d59", "icu" = "#fee090", 
                  "wicu" = "#91bfdb", "va" = "#4575b4", "imv" = "#313695")
 pal_flag     = c("TRUE" = "#d73027", "FALSE" = "#1a9850")
 
 # helper: convert ggplot to base64 PNG -----------------------------------------
 
 plot_to_base64 = function(p, width = 8, height = 5, dpi = 150) {
   
   tmp = tempfile(fileext = ".png")
   
   ggsave(tmp, plot = p, width = width, height = height, dpi = dpi, bg = "white")
   
   b64 = base64enc::base64encode(tmp)
   unlink(tmp)
   
   paste0("data:image/png;base64,", b64)
 }
 
 # helper: create collapsible section -------------------------------------------
 
 collapsible_section = function(title, content, id, open = TRUE) {
   
   tags$details(
     class = "qc-section",
     open = if (open) NA else NULL,
     tags$summary(tags$h2(title)),
     tags$div(class = "section-content", content)
   )
 }
 
 # helper: status badge ---------------------------------------------------------
 
 status_badge = function(n_issues, label = "issues") {
   
   if (n_issues == 0) {
     tags$span(class = "badge badge-ok", paste("✓ No", label))
   } else {
     tags$span(class = "badge badge-warn", paste("⚠️", n_issues, label))
   }
 }
 
 # helper: data table to HTML ---------------------------------------------------
 
 dt_to_html = function(dt, max_rows = 50, caption = NULL) {
   
   if (is.null(dt) || nrow(dt) == 0) {
     return(tags$p(class = "no-data", "No data available"))
   }
   
   dt_display = if (nrow(dt) > max_rows) dt[1:max_rows] else dt
   
   header = tags$tr(lapply(names(dt_display), function(x) tags$th(x)))
   
   rows = lapply(1:nrow(dt_display), function(i) {
     tags$tr(lapply(dt_display[i], function(x) {
       val = as.character(x)
       # highlight flags
       if (grepl("^TRUE$", val)) {
         tags$td(class = "flag-true", val)
       } else if (grepl("^FALSE$", val)) {
         tags$td(class = "flag-false", val)
       } else {
         tags$td(val)
       }
     }))
   })
   
   table_content = tagList(
     if (!is.null(caption)) tags$caption(caption),
     tags$thead(header),
     tags$tbody(rows)
   )
   
   result = tags$div(
     class = "table-wrapper",
     tags$table(class = "qc-table", table_content)
   )
   
   if (nrow(dt) > max_rows) {
     result = tagList(
       result,
       tags$p(class = "truncated", 
              paste("Showing", max_rows, "of", nrow(dt), "rows"))
     )
   }
   
   result
 }
 
 # plot: sample sizes -----------------------------------------------------------
 
 plot_sample_sizes = function(sample_sizes) {
   
   if (is.null(sample_sizes) || nrow(sample_sizes) == 0) return(NULL)
   
   dt = copy(sample_sizes)
   dt_long = melt(dt, 
                  id.vars = "site", 
                  measure.vars = c("n_cancer", "n_nocancer"),
                  variable.name = "group", 
                  value.name = "n")
   
   dt_long[, group := fifelse(group == "n_cancer", "Cancer", "No Cancer")]
   
   p = ggplot(dt_long, aes(x = reorder(site, -n), y = n, fill = group)) +
     geom_col(position = "stack", alpha = 0.9) +
     scale_fill_manual(values = c("Cancer" = "#d73027", "No Cancer" = "#4575b4")) +
     scale_y_continuous(labels = scales::comma) +
     labs(
       title = "Sample Size by Site",
       x = NULL,
       y = "Number of Encounters",
       fill = NULL
     ) +
     theme_qc() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   
   p
 }
 
 # plot: cancer prevalence ------------------------------------------------------
 
 plot_cancer_prevalence = function(sample_sizes, bounds) {
   
   if (is.null(sample_sizes) || nrow(sample_sizes) == 0) return(NULL)
   
   dt = copy(sample_sizes)
   dt[, site := reorder(site, ca_prev)]
   
   p = ggplot(dt, aes(x = ca_prev, y = site)) +
     annotate("rect", 
              xmin = bounds$ca_prev_min, xmax = bounds$ca_prev_max,
              ymin = -Inf, ymax = Inf, 
              fill = "green", alpha = 0.1) +
     geom_segment(aes(x = 0, xend = ca_prev, yend = site), color = "gray70") +
     geom_point(aes(color = ca_prev < bounds$ca_prev_min | ca_prev > bounds$ca_prev_max), 
                size = 4) +
     scale_color_manual(values = c("FALSE" = "#4575b4", "TRUE" = "#d73027"), guide = "none") +
     scale_x_continuous(labels = scales::percent, limits = c(0, NA)) +
     geom_vline(xintercept = c(bounds$ca_prev_min, bounds$ca_prev_max), 
                linetype = "dashed", color = "gray50") +
     labs(
       title = "Cancer Prevalence by Site",
       subtitle = paste0("Expected range: ", scales::percent(bounds$ca_prev_min), 
                         " - ", scales::percent(bounds$ca_prev_max)),
       x = "Cancer Prevalence",
       y = NULL
     ) +
     theme_qc()
   
   p
 }
 
 # plot: exclusion flow ---------------------------------------------------------
 
 plot_exclusion_flow = function(flow_results) {
   
   if (is.null(flow_results$site_level) || nrow(flow_results$site_level) == 0) return(NULL)
   
   dt = copy(flow_results$site_level)
   dt = dt[!is.na(pct_excluded_ca)]
   
   # simplify step names
   dt[, step_short := str_remove(step, "^After excluding ")]
   dt[, step_short := str_remove(step_short, "patients (who were |not )")]
   dt[, step_short := str_trunc(step_short, 35)]
   
   # ensure step order
   if ("step_num" %in% names(dt)) {
     dt[, step_short := factor(step_short, levels = unique(step_short[order(step_num)]))]
   }
   
   dt_long = melt(dt, 
                  id.vars = c("site", "step_short", "step_num"),
                  measure.vars = c("pct_excluded_ca", "pct_excluded_no"),
                  variable.name = "group",
                  value.name = "pct")
   
   dt_long = dt_long[!is.na(pct)]
   dt_long[, group := fifelse(group == "pct_excluded_ca", "Cancer", "No Cancer")]
   
   p = ggplot(dt_long, aes(x = pct, y = step_short, color = site, shape = group)) +
     geom_point(size = 3, alpha = 0.8, position = position_dodge(width = 0.6)) +
     scale_x_continuous(
       labels = function(x) paste0(round(x), "%"),
       limits = c(0, NA),
       expand = expansion(mult = c(0, 0.05))
     ) +
     scale_shape_manual(values = c("Cancer" = 17, "No Cancer" = 16)) +
     labs(
       title = "Exclusion Percentages by Step",
       subtitle = "Percentage of previous step's cohort excluded",
       x = "Percent Excluded",
       y = NULL,
       color = "Site",
       shape = "Group"
     ) +
     theme_qc() +
     theme(legend.position = "right")
   
   p
 }
 
 # plot: outcome rates ----------------------------------------------------------
 
 plot_outcome_rates = function(outcomes_results, bounds) {
   
   if (is.null(outcomes_results$long) || nrow(outcomes_results$long) == 0) return(NULL)
   
   dt = copy(outcomes_results$long)
   dt[, ca_label := fifelse(ca_01 == 1, "Cancer", "No Cancer")]
   
   # add pooled mean
   pooled = outcomes_results$pooled
   
   p = ggplot(dt, aes(x = rate, y = reorder(site, rate), color = var)) +
     geom_point(size = 3, alpha = 0.8) +
     geom_vline(data = pooled, aes(xintercept = mean_rate, color = var), 
                linetype = "dashed", alpha = 0.5) +
     scale_x_continuous(labels = scales::percent) +
     scale_color_manual(values = pal_outcome) +
     facet_grid(var ~ ca_label, scales = "free_x") +
     labs(
       title = "Outcome Rates by Site",
       subtitle = "Dashed lines show pooled mean",
       x = "Rate",
       y = NULL
     ) +
     theme_qc() +
     theme(legend.position = "none",
           strip.text.y = element_text(angle = 0))
   
   p
 }
 
 # plot: outcome rates comparison (forest-plot style) ---------------------------
 
 plot_outcome_forest = function(outcomes_results) {
   
   if (is.null(outcomes_results$long) || nrow(outcomes_results$long) == 0) return(NULL)
   
   dt = copy(outcomes_results$long)
   dt[, ca_label := fifelse(ca_01 == 1, "Cancer", "No Cancer")]
   
   # focus on key outcomes
   key_outcomes = c("dead", "hospice", "icu")
   dt = dt[var %in% key_outcomes]
   
   # nicer labels
   dt[, var_label := fcase(
     var == "dead", "Mortality",
     var == "hospice", "Hospice",
     var == "icu", "ICU Admission",
     default = var
   )]
   
   pooled = outcomes_results$pooled[var %in% key_outcomes]
   pooled[, var_label := fcase(
     var == "dead", "Mortality",
     var == "hospice", "Hospice",
     var == "icu", "ICU Admission",
     default = var
   )]
   
   p = ggplot(dt, aes(x = rate, y = site)) +
     geom_vline(data = pooled, aes(xintercept = mean_rate), 
                linetype = "dashed", color = "gray50", linewidth = 0.5) +
     geom_point(aes(color = ca_label, shape = ca_label), size = 3, alpha = 0.8,
                position = position_dodge(width = 0.6)) +
     scale_x_continuous(labels = scales::percent) +
     scale_color_manual(values = c("Cancer" = "#d73027", "No Cancer" = "#4575b4")) +
     scale_shape_manual(values = c("Cancer" = 17, "No Cancer" = 16)) +
     facet_wrap(~ var_label, scales = "free_x", nrow = 1) +
     labs(
       title = "Key Outcome Rates by Site and Cancer Status",
       subtitle = "Dashed line = pooled mean",
       x = "Rate",
       y = NULL,
       color = NULL,
       shape = NULL
     ) +
     theme_qc() +
     theme(
       legend.position = "bottom",
       panel.spacing = unit(1, "lines")
     )
   
   p
 }
 
 # plot: demographics -----------------------------------------------------------
 
 plot_demographics = function(demographics_results) {
   
   if (is.null(demographics_results$continuous) || 
       nrow(demographics_results$continuous) == 0) return(NULL)
   
   dt = copy(demographics_results$continuous)
   dt[, ca_label := fifelse(ca_01 == 1, "Cancer", "No Cancer")]
   
   # age plot
   p_age = ggplot(dt, aes(x = age_mean, y = reorder(site, age_mean), color = ca_label)) +
     geom_errorbarh(aes(xmin = age_mean - age_sd, xmax = age_mean + age_sd), 
                    height = 0.3, alpha = 0.5,
                    position = position_dodge(width = 0.5)) +
     geom_point(size = 3, position = position_dodge(width = 0.5)) +
     scale_color_manual(values = c("Cancer" = "#d73027", "No Cancer" = "#4575b4")) +
     labs(
       title = "Age Distribution by Site",
       subtitle = "Mean ± SD",
       x = "Age (years)",
       y = NULL,
       color = NULL
     ) +
     theme_qc()
   
   # VW plot
   p_vw = ggplot(dt, aes(x = vw_mean, y = reorder(site, vw_mean), color = ca_label)) +
     geom_errorbarh(aes(xmin = vw_mean - vw_sd, xmax = vw_mean + vw_sd), 
                    height = 0.3, alpha = 0.5,
                    position = position_dodge(width = 0.5)) +
     geom_point(size = 3, position = position_dodge(width = 0.5)) +
     scale_color_manual(values = c("Cancer" = "#d73027", "No Cancer" = "#4575b4")) +
     labs(
       title = "Van Walraven Score by Site",
       subtitle = "Mean ± SD",
       x = "Van Walraven Score",
       y = NULL,
       color = NULL
     ) +
     theme_qc()
   
   list(age = p_age, vw = p_vw)
 }
 
 # plot: heterogeneity heatmap --------------------------------------------------
 
 plot_heterogeneity = function(heterogeneity_results) {
   
   if (!is.data.table(heterogeneity_results) || 
       nrow(heterogeneity_results) == 0) return(NULL)
   
   dt = copy(heterogeneity_results)
   dt[, ca_label := fifelse(ca_01 == 1, "Cancer", "No Cancer")]
   
   p = ggplot(dt, aes(x = var, y = ca_label, fill = cv)) +
     geom_tile(color = "white", linewidth = 0.5) +
     geom_text(aes(label = sprintf("%.2f", cv)), color = "white", fontface = "bold") +
     scale_fill_gradient2(low = "#1a9850", mid = "#ffffbf", high = "#d73027",
                          midpoint = 0.3, limits = c(0, 1),
                          name = "CV") +
     labs(
       title = "Cross-Site Heterogeneity",
       subtitle = "Coefficient of Variation (CV > 0.5 flagged)",
       x = NULL,
       y = NULL
     ) +
     theme_qc() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1),
           panel.grid = element_blank())
   
   p
 }
 
 # plot: meta coefficients ------------------------------------------------------
 
 plot_meta_coefficients = function(meta_results) {
   
   if (is.null(meta_results$all_coefs) || nrow(meta_results$all_coefs) == 0) return(NULL)
   
   dt = copy(meta_results$all_coefs)
   dt[, ci_lower := beta_int - 1.96 * se_int]
   dt[, ci_upper := beta_int + 1.96 * se_int]
   
   # nicer score labels
   dt[, score_label := fcase(
     grepl("sirs", score, ignore.case = TRUE) & !grepl("mews", score, ignore.case = TRUE), "SIRS",
     grepl("qsofa", score, ignore.case = TRUE), "qSOFA",
     grepl("mews_sf|mewssf", score, ignore.case = TRUE), "MEWS-SF",
     grepl("mews", score, ignore.case = TRUE), "MEWS",
     grepl("news", score, ignore.case = TRUE), "NEWS",
     default = score
   )]
   
   # convergence as factor for legend
   dt[, converged_label := fifelse(converged == TRUE, "Converged", "Not converged")]
   
   p = ggplot(dt, aes(x = beta_int, y = site)) +
     geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
     geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = converged_label), 
                    height = 0.3, alpha = 0.7, linewidth = 0.6) +
     geom_point(aes(color = converged_label, shape = converged_label), size = 2.5) +
     scale_color_manual(
       values = c("Converged" = "#4575b4", "Not converged" = "#d73027"),
       name = NULL
     ) +
     scale_shape_manual(
       values = c("Converged" = 16, "Not converged" = 4),
       name = NULL
     ) +
     facet_wrap(~ score_label, scales = "free_x", nrow = 1) +
     labs(
       title = "Score × Cancer Interaction Coefficients",
       subtitle = "Point estimate ± 95% CI | β > 0 means score less predictive in cancer",
       x = "β (interaction)",
       y = NULL
     ) +
     theme_qc() +
     theme(
       legend.position = "bottom",
       panel.spacing = unit(1, "lines")
     )
   
   p
 }
 
 # plot: race/ethnicity distribution --------------------------------------------
 
 plot_race_distribution = function(demographics_results) {
   
   if (is.null(demographics_results$race) || 
       nrow(demographics_results$race) == 0) return(NULL)
   
   dt = copy(demographics_results$race)
   dt[, ca_label := fifelse(ca_01 == 1, "Cancer", "No Cancer")]
   
   # standardize categories
   dt[, category_clean := fcase(
     grepl("white", category, ignore.case = TRUE), "White",
     grepl("black|african", category, ignore.case = TRUE), "Black",
     grepl("asian", category, ignore.case = TRUE), "Asian",
     grepl("hispanic", category, ignore.case = TRUE), "Hispanic",
     grepl("pacific|hawaiian", category, ignore.case = TRUE), "Pacific Islander",
     grepl("native|american indian|alaska", category, ignore.case = TRUE), "Native American",
     grepl("two|multiple|more", category, ignore.case = TRUE), "Multiple",
     default = "Other/Unknown"
   )]
   
   # aggregate by standardized category
   dt_agg = dt[, .(n = sum(n, na.rm = TRUE)), by = .(site, ca_label, category_clean)]
   
   # calculate totals and proportions
   dt_agg[, N := sum(n, na.rm = TRUE), by = .(site, ca_label)]
   dt_agg[, prop := n / N]
   
   # verify proportions sum to 1 (or close)
   check_sums = dt_agg[, .(total_prop = sum(prop)), by = .(site, ca_label)]
   if (any(abs(check_sums$total_prop - 1) > 0.01)) {
     warning("Race proportions don't sum to 100% for some site/group combinations")
   }
   
   # order categories for consistent stacking
   cat_order = c("White", "Black", "Asian", "Hispanic", "Pacific Islander", 
                 "Native American", "Multiple", "Other/Unknown")
   dt_agg[, category_clean := factor(category_clean, levels = cat_order)]
   
   p = ggplot(dt_agg, aes(x = site, y = prop, fill = category_clean)) +
     geom_col(position = "stack", alpha = 0.9, color = "white", linewidth = 0.2) +
     scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0)) +
     scale_fill_brewer(palette = "Set2", name = "Race") +
     facet_wrap(~ ca_label) +
     labs(
       title = "Race Distribution by Site",
       x = NULL,
       y = "Proportion"
     ) +
     theme_qc() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1))
   
   p
 }
 
 # CSS for HTML report ----------------------------------------------------------
 
 qc_css = "
:root {
  --color-ok: #1a9850;
  --color-warn: #d73027;
  --color-info: #4575b4;
  --color-bg: #f8f9fa;
  --color-border: #dee2e6;
}

body {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
  line-height: 1.6;
  max-width: 1400px;
  margin: 0 auto;
  padding: 20px;
  background: #fff;
  color: #333;
}

h1 {
  color: var(--color-info);
  border-bottom: 3px solid var(--color-info);
  padding-bottom: 10px;
}

h2 {
  color: #495057;
  margin: 0;
  font-size: 1.3em;
}

.header-info {
  background: var(--color-bg);
  padding: 15px;
  border-radius: 8px;
  margin-bottom: 20px;
}

.header-info p {
  margin: 5px 0;
}

.summary-box {
  display: flex;
  gap: 20px;
  flex-wrap: wrap;
  margin: 20px 0;
}

.summary-card {
  background: var(--color-bg);
  border-radius: 8px;
  padding: 15px 25px;
  text-align: center;
  min-width: 150px;
}

.summary-card .number {
  font-size: 2.5em;
  font-weight: bold;
  color: var(--color-info);
}

.summary-card.warn .number {
  color: var(--color-warn);
}

.summary-card.ok .number {
  color: var(--color-ok);
}

.summary-card .label {
  color: #6c757d;
  font-size: 0.9em;
}

.qc-section {
  margin: 15px 0;
  border: 1px solid var(--color-border);
  border-radius: 8px;
  overflow: hidden;
}

.qc-section summary {
  background: var(--color-bg);
  padding: 12px 15px;
  cursor: pointer;
  user-select: none;
}

.qc-section summary:hover {
  background: #e9ecef;
}

.qc-section[open] summary {
  border-bottom: 1px solid var(--color-border);
}

.section-content {
  padding: 15px;
}

.badge {
  display: inline-block;
  padding: 4px 10px;
  border-radius: 12px;
  font-size: 0.85em;
  font-weight: 500;
  margin-left: 10px;
}

.badge-ok {
  background: #d4edda;
  color: var(--color-ok);
}

.badge-warn {
  background: #f8d7da;
  color: var(--color-warn);
}

.plot-container {
  margin: 15px 0;
  text-align: center;
}

.plot-container img {
  max-width: 100%;
  height: auto;
  border-radius: 4px;
  box-shadow: 0 2px 4px rgba(0,0,0,0.1);
}

.plot-row {
  display: flex;
  gap: 20px;
  flex-wrap: wrap;
  justify-content: center;
}

.plot-row .plot-container {
  flex: 1;
  min-width: 400px;
}

.table-wrapper {
  overflow-x: auto;
  margin: 15px 0;
}

.qc-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 0.9em;
}

.qc-table th {
  background: var(--color-bg);
  padding: 10px;
  text-align: left;
  border-bottom: 2px solid var(--color-border);
  white-space: nowrap;
}

.qc-table td {
  padding: 8px 10px;
  border-bottom: 1px solid var(--color-border);
}

.qc-table tr:hover {
  background: #f8f9fa;
}

.flag-true {
  background: #f8d7da !important;
  color: var(--color-warn);
  font-weight: bold;
}

.flag-false {
  color: var(--color-ok);
}

.no-data {
  color: #6c757d;
  font-style: italic;
  text-align: center;
  padding: 20px;
}

.truncated {
  color: #6c757d;
  font-size: 0.85em;
  text-align: right;
  font-style: italic;
}

.footnote {
  font-size: 0.85em;
  color: #6c757d;
  margin-top: 30px;
  padding-top: 15px;
  border-top: 1px solid var(--color-border);
}

@media (max-width: 768px) {
  .plot-row .plot-container {
    min-width: 100%;
  }
  
  .summary-box {
    justify-content: center;
  }
}
"

## main report generator --------------------------------------------------------

generate_qc_report = function(qc_results, 
                              output_file         = "qc_report.html",
                              title               = "Federated Data Quality Control Report",
                              plausibility_bounds = NULL) {
  
  # use global default if not provided
  if (is.null(plausibility_bounds)) plausibility_bounds = bounds
  
  cat("Generating HTML report...\n")
  
  # count issues
  n_sites = length(unique(qc_results$sample_sizes$site))
  
  n_missing_files = nrow(qc_results$structure$summary[exists == FALSE])
  n_sample_flags  = nrow(qc_results$sample_sizes[
    flag_n_low == TRUE | flag_n_high == TRUE | flag_ca_low == TRUE | flag_ca_high == TRUE
  ])
  n_flow_outliers    = if (!is.null(qc_results$flow$outliers)) nrow(qc_results$flow$outliers) else 0
  n_outcome_outliers = if (!is.null(qc_results$outcomes$outliers)) nrow(qc_results$outcomes$outliers) else 0
  n_total_issues     = n_missing_files + n_sample_flags + n_flow_outliers + n_outcome_outliers
  
  # header
  header = tags$div(
    class = "header-info",
    tags$p(tags$strong("Generated: "), format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    tags$p(tags$strong("Sites: "), paste(unique(qc_results$sample_sizes$site), collapse = ", "))
  )
  
  # summary cards
  summary_cards = tags$div(
    class = "summary-box",
    tags$div(class = "summary-card",
             tags$div(class = "number", n_sites),
             tags$div(class = "label", "Sites")),
    tags$div(class = paste("summary-card", if (n_missing_files > 0) "warn" else "ok"),
             tags$div(class = "number", n_missing_files),
             tags$div(class = "label", "Missing Files")),
    tags$div(class = paste("summary-card", if (n_total_issues > 0) "warn" else "ok"),
             tags$div(class = "number", n_total_issues),
             tags$div(class = "label", "Total Flags"))
  )
  
  # section 1: sample sizes
  cat("  Creating sample size plots...\n")
  p_sizes = plot_sample_sizes(qc_results$sample_sizes)
  p_prev  = plot_cancer_prevalence(qc_results$sample_sizes, plausibility_bounds)
  
  section_sample = collapsible_section(
    title = "Sample Sizes & Cancer Prevalence",
    id = "sample-sizes",
    content = tagList(
      status_badge(n_sample_flags, "flags"),
      tags$div(
        class = "plot-row",
        tags$div(class = "plot-container", 
                 tags$img(src = plot_to_base64(p_sizes, width = 7, height = 5))),
        tags$div(class = "plot-container", 
                 tags$img(src = plot_to_base64(p_prev, width = 6, height = 5)))
      ),
      tags$h3("Data Table"),
      dt_to_html(qc_results$sample_sizes[, .(
        site, 
        n_total = format(n_total, big.mark = ","),
        n_cancer = format(n_cancer, big.mark = ","),
        ca_prev = sprintf("%.1f%%", ca_prev * 100),
        flag_n_low, flag_n_high, flag_ca_low, flag_ca_high
      )])
    )
  )
  
  # section 2: exclusion flow
  cat("  Creating exclusion flow plot...\n")
  p_flow = plot_exclusion_flow(qc_results$flow)
  
  section_flow = collapsible_section(
    title = "Exclusion Flow",
    id = "exclusion-flow",
    content = tagList(
      status_badge(n_flow_outliers, "outliers"),
      if (!is.null(p_flow)) {
        tags$div(class = "plot-container",
                 tags$img(src = plot_to_base64(p_flow, width = 10, height = 6)))
      },
      tags$h3("Pooled Summary"),
      if (!is.null(qc_results$flow$pooled)) {
        dt_to_html(qc_results$flow$pooled[, .(
          step = stringr::str_trunc(step, 50),
          mean_pct_ca = sprintf("%.1f", mean_pct_ca),
          sd_pct_ca = sprintf("%.1f", sd_pct_ca),
          mean_pct_no = sprintf("%.1f", mean_pct_no),
          sd_pct_no = sprintf("%.1f", sd_pct_no)
        )])
      }
    )
  )
  
  # section 3: outcomes
  cat("  Creating outcome plots...\n")
  p_outcomes = plot_outcome_forest(qc_results$outcomes)
  
  section_outcomes = collapsible_section(
    title = "Outcome Rates",
    id = "outcomes",
    content = tagList(
      status_badge(n_outcome_outliers, "outliers"),
      if (!is.null(p_outcomes)) {
        tags$div(class = "plot-container",
                 tags$img(src = plot_to_base64(p_outcomes, width = 12, height = 5)))
      },
      tags$h3("Pooled Rates"),
      if (!is.null(qc_results$outcomes$pooled)) {
        dt_to_html(qc_results$outcomes$pooled[, .(
          var, ca_01,
          mean_rate = sprintf("%.1f%%", mean_rate * 100),
          sd_rate = sprintf("%.1f%%", sd_rate * 100)
        )])
      }
    )
  )
  
  # section 4: demographics
  cat("  Creating demographic plots...\n")
  p_demo = plot_demographics(qc_results$demographics)
  p_race = plot_race_distribution(qc_results$demographics)
  
  section_demo = collapsible_section(
    title = "Demographics",
    id = "demographics",
    content = tagList(
      if (!is.null(p_demo$age)) {
        tags$div(
          class = "plot-row",
          tags$div(class = "plot-container",
                   tags$img(src = plot_to_base64(p_demo$age, width = 7, height = 5))),
          tags$div(class = "plot-container",
                   tags$img(src = plot_to_base64(p_demo$vw, width = 7, height = 5)))
        )
      },
      if (!is.null(p_race)) {
        tags$div(class = "plot-container",
                 tags$img(src = plot_to_base64(p_race, width = 10, height = 5)))
      },
      tags$h3("Summary Table"),
      if (!is.null(qc_results$demographics$continuous)) {
        dt_to_html(qc_results$demographics$continuous[, .(
          site, ca_01,
          n = format(n, big.mark = ","),
          age_mean = sprintf("%.1f", age_mean),
          age_sd = sprintf("%.1f", age_sd),
          vw_mean = sprintf("%.1f", vw_mean),
          vw_sd = sprintf("%.1f", vw_sd)
        )])
      }
    )
  )
  
  # section 5: heterogeneity
  cat("  Creating heterogeneity plot...\n")
  p_het = plot_heterogeneity(qc_results$heterogeneity)
  
  n_high_cv = if (is.data.table(qc_results$heterogeneity)) {
    nrow(qc_results$heterogeneity[flag_high_cv == TRUE])
  } else 0
  
  section_het = collapsible_section(
    title = "Cross-Site Heterogeneity",
    id = "heterogeneity",
    content = tagList(
      status_badge(n_high_cv, "high CV"),
      if (!is.null(p_het)) {
        tags$div(class = "plot-container",
                 tags$img(src = plot_to_base64(p_het, width = 8, height = 5)))
      },
      tags$h3("Coefficient of Variation"),
      if (is.data.table(qc_results$heterogeneity)) {
        dt_to_html(qc_results$heterogeneity[, .(
          var, ca_01,
          mean_prop = sprintf("%.1f%%", mean_prop * 100),
          min_prop = sprintf("%.1f%%", min_prop * 100),
          max_prop = sprintf("%.1f%%", max_prop * 100),
          cv = sprintf("%.3f", cv),
          flag_high_cv
        )])
      }
    )
  )
  
  # section 6: meta-analysis (if available)
  section_meta = NULL
  if (!is.null(qc_results$meta) && !is.null(qc_results$meta$all_coefs)) {
    cat("  Creating meta-analysis plots...\n")
    p_coef = plot_meta_coefficients(qc_results$meta)
    
    n_coef_flags = if (!is.null(qc_results$meta$coef_flags)) {
      nrow(qc_results$meta$coef_flags)
    } else 0
    
    section_meta = collapsible_section(
      title = "Meta-Analysis Inputs",
      id = "meta",
      content = tagList(
        status_badge(n_coef_flags, "coefficient flags"),
        if (!is.null(p_coef)) {
          tags$div(class = "plot-container",
                   tags$img(src = plot_to_base64(p_coef, width = 12, height = 5)))
        },
        tags$h3("Coefficients by Site"),
        if (!is.null(qc_results$meta$all_coefs)) {
          dt_to_html(qc_results$meta$all_coefs[, .(
            site, score,
            beta_int = sprintf("%.3f", beta_int),
            se_int = sprintf("%.3f", se_int),
            site_n = format(site_n, big.mark = ","),
            converged
          )])
        }
      )
    )
  }
  
  # section 7: score ranges (if issues)
  section_scores = NULL
  if (!is.null(qc_results$score_ranges) && 
      is.data.table(qc_results$score_ranges) &&
      nrow(qc_results$score_ranges) > 0) {
    
    section_scores = collapsible_section(
      title = "Score Range Violations",
      id = "scores",
      content = tagList(
        status_badge(nrow(qc_results$score_ranges), "violations"),
        dt_to_html(qc_results$score_ranges)
      )
    )
  }
  
  # section 8: file inventory
  section_files = collapsible_section(
    title = "File Inventory",
    id = "files",
    open = FALSE,
    content = tagList(
      status_badge(n_missing_files, "missing"),
      dt_to_html(qc_results$structure$summary[, .(
        site, file_type, location, exists, n_rows
      )], max_rows = 100)
    )
  )
  
  # assemble full document
  doc = tags$html(
    tags$head(
      tags$meta(charset = "UTF-8"),
      tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
      tags$title(title),
      tags$style(HTML(qc_css))
    ),
    tags$body(
      tags$h1(title),
      header,
      summary_cards,
      section_sample,
      section_flow,
      section_outcomes,
      section_demo,
      section_het,
      section_meta,
      section_scores,
      section_files,
      tags$div(
        class = "footnote",
        tags$p("Generated by qc_federated_data.R"),
        tags$p(paste("R version:", R.version.string))
      )
    )
  )
  
  # write to file
  cat("  Writing HTML file...\n")
  save_html(doc, file = output_file)
  
  cat("✓ Report saved to:", normalizePath(output_file), "\n")
  
  invisible(output_file)
}

# run -------------------------------------------------------------------------

output = run_qc_with_report(
   main_folder    = here(),
   expected_sites = expected_sites,
   output_file    = here("qc_report.html")
)

# or run separately:
# qc_results = run_qc(main_folder = here())
# generate_qc_report(qc_results, output_file = "qc_report.html")