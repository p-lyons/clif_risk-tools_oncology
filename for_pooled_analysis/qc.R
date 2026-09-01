# qc.R
# Quality control script for the CLIF oncology risk project, round two.
# Run after sites upload their aggregate folders and you copy them into the
# repo root as {site}/{subdir}/file (the layout 00_load.R reads).
#
# Round-two regeneration (decision 4b):
#   - The expected-file catalog is DERIVED from the filename grammar in
#     code/03a_artifacts.R ({artifact}[-{strata}][-{outcome}][-h{horizon}]
#     [-{variant}]-{site}.csv) via a family table, so one row per artifact
#     family covers every outcome/horizon/variant combination and the catalog
#     cannot drift from the site pipeline's naming logic.
#   - Michigan added to expected sites (eight total).
#   - Meta-analysis input checks removed (coefficients / score_sds no longer
#     exist; no longitudinal model in the paper).
#   - New module: MANIFEST reconciliation (listed vs present, row counts, md5),
#     per the coordinating-center role assigned in code/run_all.R.
#   - New check: files present on disk but absent from the expected catalog are
#     flagged as stale (a site that reran without deleting upload_to_box_v2/).
#   - Round-two schemas throughout (o_out replaces o_primary_01; sesp gains
#     ppv/npv and count identities; six-row flow ending at the reconciliation
#     step; outcome variables composite/nohospice/wardicu/warddeath/hospicedc).

# setup ------------------------------------------------------------------------

library(data.table)
library(tidytable)
library(collapse)
library(stringr)
library(here)

# configuration ----------------------------------------------------------------

## expected sites (update as sites join) ---------------------------------------

expected_sites = c(
  "emory",
  "hopkins",
  "michigan",
  "ohsu",
  "rush",
  "ucmc",
  "umn",
  "upenn"
)

## filename grammar (mirrors code/03a_artifacts.R) ------------------------------

ANALYSIS_SUBDIRS = c("main", "threshold", "horizon", "sensitivity", "diagnostics")

STRATA_VOCAB  = c("ca", "liquid", "mets", "components")
OUTCOME_VOCAB = c("composite", "nohospice", "wardicu", "warddeath", "hospicedc")

SE_VARIANTS = c(
  "se_fullcode_only",
  "se_no_ed_req",
  "se_win0_96h",
  "se_one_enc_per_pt"
)
CF_VARIANTS   = c("cf2", "cf6", "cf12")
BOOT_VARIANTS = c("boot", "bootenc")
VARIANT_VOCAB = c(SE_VARIANTS, CF_VARIANTS, BOOT_VARIANTS)

#' Build an artifact filename (mirror of .build_filename in 03a_artifacts.R).
build_artifact_filename = function(artifact, site, strata = NA, outcome = NA,
                                   horizon = NA, variant = NA) {

  ok = function(x) !is.na(x) && nzchar(as.character(x))

  parts = artifact
  if (ok(strata))  parts = paste0(parts, "-", strata)
  if (ok(outcome)) parts = paste0(parts, "-", outcome)
  if (!is.na(horizon)) parts = paste0(parts, "-h", horizon)
  if (ok(variant)) parts = paste0(parts, "-", variant)

  paste0(parts, "-", site, ".csv")
}

#' Parse an artifact filename (mirror of .parse_filename in 03a_artifacts.R).
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

## expected artifact families ---------------------------------------------------
# One row per (subdir, artifact, strata); outcomes / horizons / variants list
# the combinations that family ships. NA means the token is absent from the
# filename. The expansion of this table over the eight sites IS the expected
# file catalog (135 files per site, plus ten root files and the manifest).

.fam = function(subdir, artifact, strata, outcomes = NA_character_,
                horizons = NA_integer_, variants = NA_character_) {
  CJ(
    subdir   = subdir,
    artifact = artifact,
    strata   = strata,
    outcome  = outcomes,
    horizon  = horizons,
    variant  = variants,
    sorted   = FALSE
  )
}

ALL5  = OUTCOME_VOCAB
COMP2 = c("composite", "nohospice")

expected_families = rbindlist(list(

  # main
  .fam("main", "maxscores",  "ca",     ALL5),
  .fam("main", "maxscores",  "liquid", COMP2),
  .fam("main", "maxscores",  "mets",   COMP2),
  .fam("main", "auroc",      "ca",     ALL5),
  .fam("main", "auroc",      "liquid", COMP2),
  .fam("main", "auroc",      "mets",   COMP2),
  .fam("main", "events",     "liquid", ALL5),
  .fam("main", "events",     "mets",   ALL5),
  .fam("main", "firstscore", "ca"),

  # threshold
  .fam("threshold", "ever",       "ca",         ALL5),
  .fam("threshold", "ever",       "liquid",     COMP2),
  .fam("threshold", "sesp",       "ca",         ALL5),
  .fam("threshold", "sesp",       "liquid",     COMP2),
  .fam("threshold", "cuminc",     "ca",         ALL5),
  .fam("threshold", "first",      "ca",         ALL5),
  .fam("threshold", "upset",      "ca",         ALL5),
  .fam("threshold", "upset",      "components", ALL5),
  .fam("threshold", "leadtime",        "ca",    COMP2),
  .fam("threshold", "leadtime_median", "ca",    COMP2),
  .fam("threshold", "crossclass",      "ca",    COMP2),

  # horizon
  .fam("horizon", "auroc",  "ca",     ALL5,        c(12L, 24L)),
  .fam("horizon", "auroc",  "liquid", COMP2,       24L),
  .fam("horizon", "auroc",  "mets",   COMP2,       24L),
  .fam("horizon", "counts", "ca",     ALL5,        c(12L, 24L)),
  .fam("horizon", "counts", "liquid", COMP2,       24L),
  .fam("horizon", "counts", "mets",   COMP2,       24L),
  .fam("horizon", "counts", "ca",     "composite", c(12L, 24L), "boot"),
  .fam("horizon", "counts", "ca",     "composite", c(12L, 24L), "bootenc"),

  # sensitivity
  .fam("sensitivity", "auroc",     "ca", "composite", c(NA_integer_, 12L, 24L), SE_VARIANTS),
  .fam("sensitivity", "counts",    "ca", "composite", c(12L, 24L),              SE_VARIANTS),
  .fam("sensitivity", "counts",    "ca", NA_character_, 24L,                    CF_VARIANTS),
  .fam("sensitivity", "maxscores", "ca", "composite",  NA_integer_,             SE_VARIANTS),
  .fam("sensitivity", "maxscores", "ca", NA_character_, NA_integer_,            CF_VARIANTS),

  # diagnostics
  .fam("diagnostics", "monitoring",         "ca"),
  .fam("diagnostics", "monitoring_bins",    "ca"),
  .fam("diagnostics", "missing_vlab",       "ca"),
  .fam("diagnostics", "news_o2_resolution", NA_character_)
))

#' Expand the family table into the expected file catalog for one site.
expand_expected_files = function(site) {

  out = copy(expected_families)
  out[, site := site]
  out[, filename := mapply(
    build_artifact_filename,
    artifact = artifact,
    site     = site,
    strata   = strata,
    outcome  = outcome,
    horizon  = horizon,
    variant  = variant
  )]
  out[]
}

## required columns per artifact ------------------------------------------------
# The stratum determines the grouping column; everything else follows the
# artifact (and, for auroc/counts, whether the file is horizon-level or a
# bootstrap file). Column sets taken from the verified OHSU round-two upload.

strata_group_col = function(strata) {
  fcase(
    is.na(strata),          NA_character_,
    strata == "ca",         "ca_01",
    strata == "liquid",     "liquid_01",
    strata == "mets",       "mets_01",
    strata == "components", "ca_01"
  )
}

required_cols_for = function(artifact, strata, horizon, variant) {

  g = strata_group_col(strata)

  if (artifact == "maxscores") {
    return(c("score_name", g, "max_value", "outcome", "n", "site"))
  }

  if (artifact == "auroc") {
    base = c("score_name", g, "n_obs", "n_events", "auroc", "auroc_se",
             "ci_lower", "ci_upper", "site", "estimable")
    return(c(base, if (is.na(horizon)) "metric" else "horizon"))
  }

  if (artifact == "counts") {
    base = c("score_name", g, "value", "outcome", "n", "site", "h")
    is_boot = !is.na(variant) && variant %in% BOOT_VARIANTS
    return(c(base, if (is_boot) "iter" else NULL))
  }

  if (artifact == "events") {
    return(c(g, "n_enc", "n_events", "site"))
  }

  if (artifact == "firstscore") {
    return(c("score_name", "ca_01", "value", "n", "site"))
  }

  if (artifact == "ever") {
    return(c("score_name", g, "ever_positive", "o_out", "n", "site"))
  }

  if (artifact == "sesp") {
    return(c("score_name", g, "n_total", "n_outcome", "n_no_outcome",
             "n_pos", "n_neg", "tp", "fp", "tn", "fn",
             "sensitivity", "specificity", "ppv", "npv", "site"))
  }

  if (artifact == "cuminc") {
    return(c("score", "ca_01", "o_out", "time_bin_start",
             "n_at_risk", "n_became_pos", "cum_inc", "site"))
  }

  if (artifact == "first") {
    return(c("score", "ca_01", "o_out", "n",
             "h_from_admit_sum", "h_from_admit_sumsq",
             "n_pos_day0", "n_pos_day1", "site"))
  }

  if (artifact == "upset" && !is.na(strata) && strata == "components") {
    return(c("ca_01", "o_out", "score", "component", "n", "n_encs", "site"))
  }

  if (artifact == "upset") {
    return(c("ca_01", "outcome", "sirs", "qsofa", "mews", "news", "mews_sf",
             "n", "site"))
  }

  # Round two: coarse bins plus a mean were replaced by poolable cumulative
  # counts at fixed thresholds, and by exact per-unit medians computed on
  # line-level data at the site. Both carry crossing_def (first / final_onset).
  if (artifact == "leadtime") {
    return(c("score_name", "ca_01", "outcome_key", "crossing_def",
             "threshold_h", "n_at_or_below", "n_total", "site"))
  }

  if (artifact == "leadtime_median") {
    return(c("score_name", "ca_01", "outcome_key", "crossing_def",
             "unit", "unit_id", "n_events", "median_h", "q25_h", "q75_h",
             "mean_h", "sd_h", "sum_hours", "sumsq_hours", "site"))
  }

  if (artifact == "crossclass") {
    return(c("score_name", "ca_01", "outcome_key", "crossed", "event",
             "n", "site"))
  }

  if (artifact == "monitoring") {
    return(c("ca_01", "measure", "n_enc", "sum_obs", "sumsq_obs",
             "sum_ward_hours", "sum_rate_per24h", "sumsq_rate_per24h", "site"))
  }

  if (artifact == "monitoring_bins") {
    return(c("ca_01", "measure", "rate_bin", "n_enc", "site"))
  }

  if (artifact == "missing_vlab") {
    return(c("variable", "ca_01", "n_total", "n_missing", "pct_missing", "site"))
  }

  if (artifact == "news_o2_resolution") {
    return(c("ca_01", "resolution_branch", "n", "site"))
  }

  character(0)
}

## file specifications: root files ----------------------------------------------

file_specs_root = list(

  flow = list(
    pattern       = "figure_s01_flow_{site}.csv",
    required_cols = c("step", "n_remaining_ca", "n_excluded_ca",
                      "n_remaining_no", "n_excluded_no")
  ),

  cat = list(
    pattern       = "table_02_cat_{site}.csv",
    required_cols = c("ca_01", "var", "category", "n")
  ),

  cont = list(
    pattern       = "table_02_cont_{site}.csv",
    required_cols = c("ca_01", "n", "age_sum", "age_sumsq", "vw_sum", "vw_sumsq",
                      "los_sum", "los_sumsq", "site")
  ),

  cancer_codes = list(
    pattern       = "cancer_codes_primary_{site}.csv",
    required_cols = c("ca_icd10_enc", "n", "site")
  ),

  hospital_types = list(
    pattern       = "hospital_types_{site}.csv",
    required_cols = c("hospital_type", "n_hospitals", "site")
  ),

  run_report = list(
    pattern       = "run_report_{site}.csv",
    required_cols = c("site", "n_encounters", "n_cancer", "n_ed_admit",
                      "pct_cancer", "outcome_rate_pct", "n_score_rows",
                      "n_hid_multi_pt", "completed")
  ),

  missing_demog = list(
    pattern       = "missing_demog_{site}.csv",
    required_cols = c("variable", "n_missing", "n_total", "pct_missing", "site")
  ),

  env_manifest = list(
    pattern       = "env_manifest_{site}.csv",
    required_cols = c("item", "value", "type", "site")
  ),

  admission_dx_chapter = list(
    pattern       = "admission_dx_chapter-ca-{site}.csv",
    required_cols = c("ca_01", "chapter", "n", "site")
  ),

  admission_dx_stem = list(
    pattern       = "admission_dx_stem-ca-{site}.csv",
    required_cols = c("ca_01", "code_stem", "n", "site")
  )
)

manifest_spec = list(
  pattern       = "MANIFEST-{site}.csv",
  required_cols = c("relative_path", "n_rows", "n_cols", "md5",
                    "site", "run_date", "pipeline_version")
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

  # composite rate identity vs run_report (run_report rounds to 0.1%)
  outcome_rate_tol = 0.0015
)

## score-specific valid ranges -------------------------------------------------

score_ranges_short = list(
  sirs    = c(min = 0L, max = 4L),
  qsofa   = c(min = 0L, max = 3L),
  mews    = c(min = 0L, max = 14L),
  news    = c(min = 0L, max = 20L),
  mews_sf = c(min = 0L, max = 17L)
)

## expected flow steps (round two: six rows, ending at the reconciliation) ------
# Alternate wordings retained for cross-site label drift; matching is partial.

flow_step_order = c(
  "Adult inpatient admissions during study period",
  "After excluding patients not admitted through the ED",
  "After excluding patients who were in the ICU before hitting the wards",
  "After excluding patients who were in the ICU before admission to the wards",
  "After excluding encounters with < 6h data",
  "After excluding encounters with outcomes too early",
  "After excluding encounters with no calculable score before the outcome"
)

FLOW_N_STEPS      = 6L
FLOW_FINAL_MARKER = "no calculable score"

## outlier detection threshold (SDs from pooled mean) --------------------------

outlier_sd_threshold = 2.5

# helper functions -------------------------------------------------------------

## null coalescing operator ----------------------------------------------------

`%||%` = function(a, b) if (is.null(a)) b else a

## discover sites from folder structure ----------------------------------------

discover_sites = function(main_folder) {

  folders = list.dirs(main_folder, recursive = FALSE, full.names = FALSE)

  folders = folders[!str_starts(folders, "\\.")]
  folders = folders[folders %in% expected_sites |
                      !folders %in% c("code", "config", "for_pooled_analysis",
                                      "proj_tables", "proj_output", "qc_output",
                                      "upload_to_box", "upload_to_box_v2",
                                      "main", "threshold", "horizon",
                                      "sensitivity", "diagnostics")]
  folders
}

## read a file with schema check ------------------------------------------------

read_checked_file = function(filepath, required_cols) {

  if (!file.exists(filepath)) {
    return(list(success = FALSE, error = paste("File not found:", filepath), data = NULL))
  }

  tryCatch({
    dt = fread(filepath)

    missing_cols = setdiff(required_cols, names(dt))
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

## read a root file for a site --------------------------------------------------

read_root_file = function(main_folder, site, file_spec) {

  filename = str_replace(file_spec$pattern, "\\{site\\}", site)
  filepath = file.path(main_folder, site, filename)

  read_checked_file(filepath, file_spec$required_cols)
}

## read an expected grammar file for a site -------------------------------------

read_grammar_file = function(main_folder, site, subdir, filename,
                             artifact, strata, horizon, variant) {

  filepath = file.path(main_folder, site, subdir, filename)
  read_checked_file(filepath, required_cols_for(artifact, strata, horizon, variant))
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
# Checks every expected file (grammar-derived plus root) for presence and
# schema, and flags every on-disk file that the grammar does not expect
# (stale round-one, dxgroup, or meta files).

qc_structure = function(main_folder, sites) {

  expected_results = list()
  unexpected_rows  = list()

  for (site in sites) {

    ## grammar-derived files
    expected = expand_expected_files(site)

    for (i in seq_len(nrow(expected))) {

      result = read_grammar_file(
        main_folder,
        site,
        expected$subdir[i],
        expected$filename[i],
        expected$artifact[i],
        expected$strata[i],
        expected$horizon[i],
        expected$variant[i]
      )

      expected_results[[length(expected_results) + 1]] = tidytable(
        site      = site,
        file_type = paste(expected$subdir[i], expected$artifact[i], sep = "/"),
        location  = expected$subdir[i],
        filename  = expected$filename[i],
        exists    = result$success,
        error     = result$error %||% NA_character_,
        n_rows    = if (result$success) nrow(result$data) else NA_integer_,
        n_cols    = if (result$success) ncol(result$data) else NA_integer_
      )
    }

    ## root files
    for (file_type in names(file_specs_root)) {

      spec   = file_specs_root[[file_type]]
      result = read_root_file(main_folder, site, spec)

      expected_results[[length(expected_results) + 1]] = tidytable(
        site      = site,
        file_type = file_type,
        location  = "root",
        filename  = str_replace(spec$pattern, "\\{site\\}", site),
        exists    = result$success,
        error     = result$error %||% NA_character_,
        n_rows    = if (result$success) nrow(result$data) else NA_integer_,
        n_cols    = if (result$success) ncol(result$data) else NA_integer_
      )
    }

    ## unexpected files: anything on disk the grammar and root specs omit
    expected_paths = c(
      file.path(expected$subdir, expected$filename),
      vapply(file_specs_root, function(s) {
        str_replace(s$pattern, "\\{site\\}", site)
      }, character(1)),
      str_replace(manifest_spec$pattern, "\\{site\\}", site)
    )

    on_disk = list.files(
      file.path(main_folder, site),
      pattern    = "\\.csv$",
      recursive  = TRUE,
      full.names = FALSE
    )

    extra = setdiff(on_disk, expected_paths)
    for (f in extra) {
      unexpected_rows[[length(unexpected_rows) + 1]] = tidytable(
        site   = site,
        file   = f,
        reason = "not in the round-two expected catalog (stale upload?)"
      )
    }
  }

  summary_dt = rbindlist(expected_results)
  unexpected = if (length(unexpected_rows) > 0) rbindlist(unexpected_rows) else NULL

  list(summary = summary_dt, unexpected = unexpected)
}

# QC module 2: sample sizes and cancer prevalence ------------------------------

qc_sample_sizes = function(main_folder, sites, bounds) {

  results = list()

  for (site in sites) {

    cont_result = read_root_file(main_folder, site, file_specs_root$cont)

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

qc_flow = function(main_folder, sites, bounds, outlier_threshold) {

  all_flow = list()

  for (site in sites) {
    result = read_root_file(main_folder, site, file_specs_root$flow)
    if (result$success) {
      dt      = copy(result$data)
      dt$site = site
      all_flow[[site]] = dt
    }
  }

  if (length(all_flow) == 0) {
    return(list(site_level = NULL, outliers = NULL, pooled = NULL, structure = NULL))
  }

  flow_combined = rbindlist(all_flow, fill = TRUE)

  # structural check: six rows per site, final row is the reconciliation step
  structure_check = flow_combined[, .(
    n_steps   = .N,
    has_final = any(grepl(FLOW_FINAL_MARKER, step, ignore.case = TRUE))
  ), by = site]
  structure_check[, flag_structure := n_steps != FLOW_N_STEPS | has_final == FALSE]

  # normalize step names and assign numeric order
  flow_combined[, step_clean := str_squish(tolower(step))]
  step_order_clean = str_squish(tolower(flow_step_order))

  flow_combined[, step_num := match(step_clean, step_order_clean)]

  # for any unmatched steps, try partial matching
  flow_combined[is.na(step_num), step_num := {
    sapply(step_clean[is.na(step_num)], function(s) {
      matches = which(sapply(step_order_clean, function(so) {
        grepl(substr(s, 1, 30), so, fixed = TRUE) |
          grepl(substr(so, 1, 30), s, fixed = TRUE)
      }))
      if (length(matches) > 0) matches[1] else NA_integer_
    })
  }]

  # fallback: assign order by remaining count (descending)
  flow_combined[is.na(step_num),
                step_num := frank(-n_remaining_ca, ties.method = "dense") + 100L,
                by = site]

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

  # sanity check: cap at 100%
  flow_combined[pct_excluded_ca > 100, pct_excluded_ca := NA_real_]
  flow_combined[pct_excluded_no > 100, pct_excluded_no := NA_real_]

  # internal identity: remaining + excluded at each step must equal the
  # previous step's remaining
  flow_combined[, `:=`(
    flag_arith_ca = !is.na(prev_remaining_ca) &
      (n_remaining_ca + n_excluded_ca != prev_remaining_ca),
    flag_arith_no = !is.na(prev_remaining_no) &
      (n_remaining_no + n_excluded_no != prev_remaining_no)
  )]

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
  flow_combined = merge(
    flow_combined,
    pooled[, .(step, mean_pct_ca, sd_pct_ca, mean_pct_no, sd_pct_no)],
    by = "step", all.x = TRUE
  )

  flow_combined[, `:=`(
    z_ca = (pct_excluded_ca - mean_pct_ca) / sd_pct_ca,
    z_no = (pct_excluded_no - mean_pct_no) / sd_pct_no
  )]

  # the ED-admission step removes most encounters by design at every site, so
  # it is exempt from the high-exclusion bound (the bound exists to catch
  # extraction errors at the other steps)
  flow_combined[, `:=`(
    flag_outlier_ca = abs(z_ca) > outlier_threshold & !is.na(z_ca),
    flag_outlier_no = abs(z_no) > outlier_threshold & !is.na(z_no),
    flag_high_excl  = pmax(pct_excluded_ca, pct_excluded_no, na.rm = TRUE) > bounds$excl_pct_max &
      !grepl("not admitted through the ED", step, ignore.case = TRUE)
  )]

  outliers = flow_combined[
    flag_outlier_ca == TRUE | flag_outlier_no == TRUE | flag_high_excl == TRUE |
      flag_arith_ca == TRUE | flag_arith_no == TRUE
  ]

  list(
    site_level = flow_combined,
    outliers   = outliers,
    pooled     = pooled,
    structure  = structure_check
  )
}

# QC module 4: outcome rates ---------------------------------------------------
# Round-two table_02_cat carries composite/nohospice/wardicu/warddeath/hospicedc
# flags plus dead, hospice, icu, va, imv.

qc_outcomes = function(main_folder, sites, bounds, outlier_threshold) {

  all_cat = list()

  for (site in sites) {
    result = read_root_file(main_folder, site, file_specs_root$cat)
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

  # outcome variables (round two)
  outcome_vars = c("composite", "nohospice", "wardicu", "warddeath", "hospicedc",
                   "dead", "hospice", "icu", "va", "imv")

  outcomes = cat_combined[var %in% outcome_vars & category == "1"]
  outcomes = merge(outcomes, totals, by = c("site", "ca_01"))
  outcomes[, rate := n / N]

  # pivot wider for easier comparison
  outcomes_wide = dcast(outcomes, site + ca_01 ~ var, value.var = "rate")

  # flag implausible rates against the retained clinical bounds
  outcomes_wide[, `:=`(
    flag_mort_low  = !is.na(dead) & dead < bounds$mortality_min,
    flag_mort_high = !is.na(dead) & dead > bounds$mortality_max,
    flag_icu_low   = !is.na(icu) & icu < bounds$icu_min,
    flag_icu_high  = !is.na(icu) & icu > bounds$icu_max,
    flag_hosp_low  = !is.na(hospice) & hospice < bounds$hospice_min,
    flag_hosp_high = !is.na(hospice) & hospice > bounds$hospice_max
  )]

  # component decomposition identity: composite counts equal the sum of the
  # three components (exact by construction in 02_scores.R)
  comp_counts = dcast(
    cat_combined[var %in% c("composite", "wardicu", "warddeath", "hospicedc") &
                   category == "1"],
    site + ca_01 ~ var, value.var = "n"
  )
  comp_counts[, flag_decomp := composite != (wardicu + warddeath + hospicedc)]

  # pooled stats for outlier detection
  pooled = outcomes[, .(
    mean_rate = mean(rate, na.rm = TRUE),
    sd_rate   = sd(rate, na.rm = TRUE)
  ), by = .(var, ca_01)]

  outcomes = merge(outcomes, pooled, by = c("var", "ca_01"))
  outcomes[, z := (rate - mean_rate) / sd_rate]
  outcomes[, flag_outlier := abs(z) > outlier_threshold & !is.na(z)]

  list(
    long          = outcomes,
    wide          = outcomes_wide,
    pooled        = pooled,
    outliers      = outcomes[flag_outlier == TRUE],
    decomposition = comp_counts
  )
}

# QC module 5: demographics ----------------------------------------------------

qc_demographics = function(main_folder, sites, bounds, outlier_threshold) {

  all_cont = list()
  all_cat  = list()

  for (site in sites) {

    cont_result = read_root_file(main_folder, site, file_specs_root$cont)
    if (cont_result$success) {
      dt = copy(cont_result$data)
      if (!"site" %in% names(dt)) dt$site = site
      all_cont[[site]] = dt
    }

    cat_result = read_root_file(main_folder, site, file_specs_root$cat)
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

qc_heterogeneity = function(main_folder, sites) {

  all_cat = list()

  for (site in sites) {
    result = read_root_file(main_folder, site, file_specs_root$cat)
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

  # key proportions to check (round-two outcome set)
  key_vars = c("composite", "wardicu", "warddeath", "hospicedc",
               "dead", "hospice", "icu", "va", "imv", "female")

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
# Checks the main composite maxscores and the composite horizon counts.

qc_score_ranges = function(main_folder, sites, score_ranges_short) {

  results = list()

  for (site in sites) {

    ## encounter-level maxscores
    max_fn = build_artifact_filename("maxscores", site, strata = "ca",
                                     outcome = "composite")
    result = read_grammar_file(main_folder, site, "main", max_fn,
                               "maxscores", "ca", NA_integer_, NA_character_)

    if (result$success) {

      dt = copy(result$data)
      dt[, score_clean := str_remove(score_name, "_total")]

      for (sc in names(score_ranges_short)) {

        range_info = score_ranges_short[[sc]]
        score_rows = dt[score_clean == sc]

        if (nrow(score_rows) == 0) next

        min_val = min(score_rows$max_value, na.rm = TRUE)
        max_val = max(score_rows$max_value, na.rm = TRUE)

        if (min_val < range_info["min"] || max_val > range_info["max"]) {
          results[[paste(site, sc, "enc", sep = "_")]] = tidytable(
            site         = site,
            score        = sc,
            source       = "maxscores (encounter)",
            observed_min = min_val,
            observed_max = max_val,
            expected_min = range_info["min"],
            expected_max = range_info["max"],
            flag         = TRUE
          )
        }
      }
    }

    ## horizon counts
    for (h in c(12L, 24L)) {

      cnt_fn = build_artifact_filename("counts", site, strata = "ca",
                                       outcome = "composite", horizon = h)
      result = read_grammar_file(main_folder, site, "horizon", cnt_fn,
                                 "counts", "ca", h, NA_character_)

      if (!result$success) next

      dt = copy(result$data)
      dt[, score_clean := str_remove(score_name, "_total")]

      for (sc in names(score_ranges_short)) {

        range_info = score_ranges_short[[sc]]
        score_rows = dt[score_clean == sc]

        if (nrow(score_rows) == 0) next

        min_val = min(score_rows$value, na.rm = TRUE)
        max_val = max(score_rows$value, na.rm = TRUE)

        if (min_val < range_info["min"] || max_val > range_info["max"]) {
          results[[paste(site, sc, "h", h, sep = "_")]] = tidytable(
            site         = site,
            score        = sc,
            source       = paste0("counts (h", h, ")"),
            observed_min = min_val,
            observed_max = max_val,
            expected_min = range_info["min"],
            expected_max = range_info["max"],
            flag         = TRUE
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
# Composite files only; the per-outcome twins share their generating code.

qc_threshold_outputs = function(main_folder, sites, bounds) {

  results = list()

  ## sensitivity/specificity: plausibility + count identities
  sesp_results = list()

  for (site in sites) {

    fn     = build_artifact_filename("sesp", site, strata = "ca", outcome = "composite")
    result = read_grammar_file(main_folder, site, "threshold", fn,
                               "sesp", "ca", NA_integer_, NA_character_)

    if (!result$success) next

    dt      = copy(result$data)
    dt$site = site

    dt[, `:=`(
      flag_sens_low  = sensitivity < bounds$sens_min,
      flag_sens_high = sensitivity > bounds$sens_max,
      flag_spec_low  = specificity < bounds$spec_min,
      flag_spec_high = specificity > bounds$spec_max,
      flag_ppv       = ppv < 0 | ppv > 1,
      flag_npv       = npv < 0 | npv > 1,
      flag_cells     = abs((tp + fn + fp + tn) - n_total) > 1,
      flag_margins   = (tp + fn != n_outcome) | (fp + tn != n_no_outcome) |
                       (tp + fp != n_pos)    | (fn + tn != n_neg)
    )]

    flagged = dt[
      flag_sens_low == TRUE | flag_sens_high == TRUE |
        flag_spec_low == TRUE | flag_spec_high == TRUE |
        flag_ppv == TRUE | flag_npv == TRUE |
        flag_cells == TRUE | flag_margins == TRUE
    ]

    if (nrow(flagged) > 0) {
      sesp_results[[site]] = flagged
    }
  }

  results$sesp = if (length(sesp_results) > 0) rbindlist(sesp_results, fill = TRUE) else NULL

  ## cumulative incidence: bounded within [0, 1]
  # No monotonicity check: the site estimator is bin-wise positivity among
  # encounters still followed at each bin start (03a_artifacts.R), so the
  # at-risk denominator shrinks over time and late-bin decreases are a
  # property of the estimator, not a data defect.
  cuminc_results = list()

  for (site in sites) {

    fn     = build_artifact_filename("cuminc", site, strata = "ca", outcome = "composite")
    result = read_grammar_file(main_folder, site, "threshold", fn,
                               "cuminc", "ca", NA_integer_, NA_character_)

    if (!result$success) next

    dt      = copy(result$data)
    dt$site = site

    dt[, flag_bounds := cum_inc < 0 | cum_inc > 1]

    flagged = dt[flag_bounds == TRUE]

    if (nrow(flagged) > 0) {
      cuminc_results[[site]] = flagged[, .(site, score, ca_01, o_out,
                                           time_bin_start, cum_inc,
                                           flag_bounds)]
    }
  }

  results$cuminc = if (length(cuminc_results) > 0) rbindlist(cuminc_results, fill = TRUE) else NULL

  ## time to first positive
  first_results = list()

  for (site in sites) {

    fn     = build_artifact_filename("first", site, strata = "ca", outcome = "composite")
    result = read_grammar_file(main_folder, site, "threshold", fn,
                               "first", "ca", NA_integer_, NA_character_)

    if (!result$success) next

    dt      = copy(result$data)
    dt$site = site

    dt[, mean_time := h_from_admit_sum / n]
    dt[, flag_time := mean_time < bounds$time_to_pos_min |
                        mean_time > bounds$time_to_pos_max]

    flagged = dt[flag_time == TRUE]

    if (nrow(flagged) > 0) {
      first_results[[site]] = flagged
    }
  }

  results$first = if (length(first_results) > 0) rbindlist(first_results, fill = TRUE) else NULL

  results
}

# QC module 9: sensitivity analysis consistency --------------------------------
# se variants follow direction rules (no_ed_req is a superset; the others are
# subsets); carry-forward variants score the identical cohort and must match
# the main encounter counts exactly.

qc_sensitivity_analyses = function(main_folder, sites) {

  results = list()

  for (site in sites) {

    main_fn = build_artifact_filename("maxscores", site, strata = "ca",
                                      outcome = "composite")
    main_result = read_grammar_file(main_folder, site, "main", main_fn,
                                    "maxscores", "ca", NA_integer_, NA_character_)

    if (!main_result$success) next

    main_n = main_result$data[, .(n_main = sum(n)), by = .(score_name, ca_01)]

    ## se variants
    for (variant in SE_VARIANTS) {

      var_fn = build_artifact_filename("maxscores", site, strata = "ca",
                                       outcome = "composite", variant = variant)
      var_result = read_grammar_file(main_folder, site, "sensitivity", var_fn,
                                     "maxscores", "ca", NA_integer_, variant)

      if (!var_result$success) next

      var_n = var_result$data[, .(n_var = sum(n)), by = .(score_name, ca_01)]

      comparison = merge(main_n, var_n, by = c("score_name", "ca_01"), all = TRUE)
      comparison[, ratio := n_var / n_main]
      comparison[, site := site]
      comparison[, variant := variant]

      if (variant == "se_no_ed_req") {
        comparison[, flag := n_var < n_main * 0.95]   # should be a superset
      } else {
        comparison[, flag := n_var > n_main * 1.05]   # should not be larger
      }

      results[[paste(site, variant, sep = "_")]] = comparison
    }

    ## carry-forward variants: exact cohort identity
    for (variant in CF_VARIANTS) {

      var_fn = build_artifact_filename("maxscores", site, strata = "ca",
                                       variant = variant)
      var_result = read_grammar_file(main_folder, site, "sensitivity", var_fn,
                                     "maxscores", "ca", NA_integer_, variant)

      if (!var_result$success) next

      var_n = var_result$data[, .(n_var = sum(n)), by = .(score_name, ca_01)]

      comparison = merge(main_n, var_n, by = c("score_name", "ca_01"), all = TRUE)
      comparison[, ratio := n_var / n_main]
      comparison[, site := site]
      comparison[, variant := variant]
      comparison[, flag := n_var != n_main]

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

# QC module 10: internal consistency checks ------------------------------------
# Encounter counts from table_02_cont, the flow final row, and run_report must
# agree exactly (02_scores.R sets all three from the reconciled cohort). The
# composite rate from table_02_cat must reproduce run_report outcome_rate_pct.

qc_internal_consistency = function(main_folder, sites, bounds) {

  results = list()

  for (site in sites) {

    # source 1: table_02_cont
    cont_result = read_root_file(main_folder, site, file_specs_root$cont)
    n_cont = if (cont_result$success) cont_result$data[, sum(n)] else NA_integer_

    # source 2: flow final row
    flow_result = read_root_file(main_folder, site, file_specs_root$flow)
    if (flow_result$success) {
      flow_dt   = flow_result$data
      final_row = flow_dt[grepl(FLOW_FINAL_MARKER, step, ignore.case = TRUE)]
      n_flow = if (nrow(final_row) == 1) {
        final_row$n_remaining_ca + final_row$n_remaining_no
      } else {
        NA_integer_
      }
    } else {
      n_flow = NA_integer_
    }

    # source 3: run_report n_ed_admit and outcome rate
    report_result = read_root_file(main_folder, site, file_specs_root$run_report)
    if (report_result$success) {
      n_report    = report_result$data$n_ed_admit[1]
      rate_report = report_result$data$outcome_rate_pct[1] / 100
    } else {
      n_report    = NA_integer_
      rate_report = NA_real_
    }

    # composite rate from table_02_cat
    cat_result = read_root_file(main_folder, site, file_specs_root$cat)
    if (cat_result$success && !is.na(n_cont)) {
      n_composite = cat_result$data[var == "composite" & category == "1", sum(n)]
      rate_cat    = n_composite / n_cont
    } else {
      rate_cat = NA_real_
    }

    n_values = c(cont = n_cont, flow = n_flow, report = n_report)
    n_values = n_values[!is.na(n_values)]

    flag_n = if (length(n_values) >= 2) {
      max(n_values) != min(n_values)
    } else {
      NA
    }

    flag_rate = if (!is.na(rate_cat) && !is.na(rate_report)) {
      abs(rate_cat - rate_report) > bounds$outcome_rate_tol
    } else {
      NA
    }

    results[[site]] = tidytable(
      site        = site,
      n_cont      = n_cont,
      n_flow      = n_flow,
      n_report    = n_report,
      rate_cat    = rate_cat,
      rate_report = rate_report,
      flag_n      = flag_n,
      flag_rate   = flag_rate,
      flag        = isTRUE(flag_n) | isTRUE(flag_rate)
    )
  }

  if (length(results) == 0) return(NULL)

  rbindlist(results, fill = TRUE)
}

# QC module 11: MANIFEST reconciliation -----------------------------------------
# Compares each site's shipped MANIFEST against the files actually present:
# listed-but-absent, present-but-unlisted, and row-count / md5 mismatches.

qc_manifest = function(main_folder, sites) {

  summaries = list()
  flag_rows = list()

  for (site in sites) {

    manifest_result = read_root_file(main_folder, site, manifest_spec)

    if (!manifest_result$success) {
      summaries[[site]] = tidytable(
        site  = site,
        error = manifest_result$error
      )
      next
    }

    manifest = copy(manifest_result$data)

    on_disk = list.files(
      file.path(main_folder, site),
      pattern    = "\\.csv$",
      recursive  = TRUE,
      full.names = FALSE
    )
    on_disk = on_disk[!grepl("^MANIFEST-", basename(on_disk))]

    listed_missing   = setdiff(manifest$relative_path, on_disk)
    present_unlisted = setdiff(on_disk, manifest$relative_path)

    for (f in listed_missing) {
      flag_rows[[length(flag_rows) + 1]] = tidytable(
        site = site, file = f, issue = "listed in MANIFEST but absent on disk"
      )
    }
    for (f in present_unlisted) {
      flag_rows[[length(flag_rows) + 1]] = tidytable(
        site = site, file = f, issue = "present on disk but absent from MANIFEST"
      )
    }

    # verify row counts and checksums for the intersection
    common = intersect(manifest$relative_path, on_disk)
    n_nrow_mismatch = 0L
    n_md5_mismatch  = 0L

    for (f in common) {

      fp       = file.path(main_folder, site, f)
      expected = manifest[relative_path == f]

      actual_md5 = unname(tools::md5sum(fp))
      if (!identical(actual_md5, expected$md5[1])) {
        n_md5_mismatch = n_md5_mismatch + 1L
        flag_rows[[length(flag_rows) + 1]] = tidytable(
          site = site, file = f, issue = "md5 checksum mismatch"
        )
        next   # a changed file will usually also change row counts; one flag suffices
      }

      actual_nrow = nrow(fread(fp, showProgress = FALSE))
      if (!is.na(expected$n_rows[1]) && actual_nrow != expected$n_rows[1]) {
        n_nrow_mismatch = n_nrow_mismatch + 1L
        flag_rows[[length(flag_rows) + 1]] = tidytable(
          site = site, file = f, issue = "row count differs from MANIFEST"
        )
      }
    }

    summaries[[site]] = tidytable(
      site             = site,
      n_listed         = nrow(manifest),
      n_on_disk        = length(on_disk),
      n_listed_missing = length(listed_missing),
      n_unlisted       = length(present_unlisted),
      n_md5_mismatch   = n_md5_mismatch,
      n_nrow_mismatch  = n_nrow_mismatch,
      pipeline_version = manifest$pipeline_version[1],
      run_date         = as.character(manifest$run_date[1]),
      error            = NA_character_
    )
  }

  list(
    summary = rbindlist(summaries, fill = TRUE),
    flags   = if (length(flag_rows) > 0) rbindlist(flag_rows) else NULL
  )
}

# main QC runner ---------------------------------------------------------------

run_qc = function(main_folder            = here(),
                  sites_expected         = expected_sites,
                  plausibility_bounds    = bounds,
                  score_ranges           = score_ranges_short,
                  outlier_threshold      = outlier_sd_threshold,
                  check_analysis_outputs = TRUE) {

  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("FEDERATED DATA QUALITY CONTROL REPORT (round two)\n")
  cat(strrep("=", 70), "\n", sep = "")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Main folder:", main_folder, "\n\n")

  # discover or use expected sites
  discovered_sites = discover_sites(main_folder)

  if (is.null(sites_expected) || length(sites_expected) == 0) {
    sites = discovered_sites
    cat("Discovered", length(sites), "site folders:", paste(sites, collapse = ", "), "\n\n")
  } else {
    missing_sites = setdiff(sites_expected, discovered_sites)
    extra_sites   = setdiff(discovered_sites, sites_expected)
    sites         = intersect(sites_expected, discovered_sites)

    if (length(missing_sites) > 0) {
      cat("⚠️  MISSING EXPECTED SITES:", paste(missing_sites, collapse = ", "), "\n")
    }
    if (length(extra_sites) > 0) {
      cat("ℹ️  Extra folders found (ignored):", paste(extra_sites, collapse = ", "), "\n")
    }
    cat("Processing", length(sites), "sites:", paste(sites, collapse = ", "), "\n\n")
  }

  if (length(sites) == 0) {
    stop("No sites to process!")
  }

  results  = list()
  n_issues = 0

  # ---------------------------------------------------------------------------
  print_header(1, "structural integrity")

  results$structure = qc_structure(main_folder, sites)

  missing_files = results$structure$summary[exists == FALSE]
  if (nrow(missing_files) > 0) {
    cat("⚠️  Missing or malformed expected files:\n")
    print(missing_files[, .(site, location, filename, error)][1:min(20, .N)])
    if (nrow(missing_files) > 20) cat("... and", nrow(missing_files) - 20, "more\n")
    n_issues = n_issues + nrow(missing_files)
  } else {
    cat("✓ All", nrow(results$structure$summary), "expected files present with",
        "required columns\n")
  }

  if (!is.null(results$structure$unexpected)) {
    cat("\n⚠️  Unexpected files on disk (stale uploads must be deleted and rerun):\n")
    print(results$structure$unexpected)
    n_issues = n_issues + nrow(results$structure$unexpected)
  } else {
    cat("✓ No unexpected files on disk\n")
  }

  # ---------------------------------------------------------------------------
  print_header(2, "MANIFEST reconciliation")

  results$manifest = qc_manifest(main_folder, sites)

  print(results$manifest$summary)

  if (!is.null(results$manifest$flags) && nrow(results$manifest$flags) > 0) {
    cat("\n⚠️  MANIFEST discrepancies:\n")
    print(results$manifest$flags[1:min(20, .N)])
    if (nrow(results$manifest$flags) > 20) {
      cat("... and", nrow(results$manifest$flags) - 20, "more\n")
    }
    n_issues = n_issues + nrow(results$manifest$flags)
  } else {
    cat("\n✓ Every MANIFEST reconciles with the files on disk\n")
  }

  # ---------------------------------------------------------------------------
  print_header(3, "sample sizes & cancer prevalence")

  results$sample_sizes = qc_sample_sizes(main_folder, sites, plausibility_bounds)

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

  # ---------------------------------------------------------------------------
  print_header(4, "exclusion flow")

  results$flow = qc_flow(main_folder, sites, plausibility_bounds, outlier_threshold)

  if (!is.null(results$flow$structure)) {
    bad_structure = results$flow$structure[flag_structure == TRUE]
    if (nrow(bad_structure) > 0) {
      cat("⚠️  Flow files without six rows ending at the reconciliation step:\n")
      print(bad_structure)
      n_issues = n_issues + nrow(bad_structure)
    } else {
      cat("✓ Six flow rows per site, reconciliation step present\n\n")
    }
  }

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
    cat("\n⚠️  Flow flags (outlier %, arithmetic identity):\n")
    print(results$flow$outliers[, .(
      site,
      step   = str_trunc(step, 30),
      pct_ca = fmt_dec(pct_excluded_ca, 1),
      pct_no = fmt_dec(pct_excluded_no, 1),
      flag_arith_ca, flag_arith_no
    )])
    n_issues = n_issues + nrow(results$flow$outliers)
  } else {
    cat("\n✓ No flow flags\n")
  }

  # ---------------------------------------------------------------------------
  print_header(5, "outcome rates")

  results$outcomes = qc_outcomes(main_folder, sites, plausibility_bounds, outlier_threshold)

  if (!is.null(results$outcomes$pooled)) {
    cat("Pooled outcome rates:\n")
    print(results$outcomes$pooled[, .(
      var, ca_01,
      mean = fmt_pct(mean_rate),
      sd   = fmt_pct(sd_rate)
    )])
  }

  if (!is.null(results$outcomes$decomposition)) {
    bad_decomp = results$outcomes$decomposition[flag_decomp == TRUE]
    if (nrow(bad_decomp) > 0) {
      cat("\n⚠️  Composite ≠ wardicu + warddeath + hospicedc:\n")
      print(bad_decomp)
      n_issues = n_issues + nrow(bad_decomp)
    } else {
      cat("\n✓ Composite decomposes exactly into its three components\n")
    }
  }

  if (!is.null(results$outcomes$outliers) && nrow(results$outcomes$outliers) > 0) {
    cat("\n⚠️  Outlier outcome rates:\n")
    print(results$outcomes$outliers[, .(site, ca_01, var, rate = fmt_pct(rate))])
    n_issues = n_issues + nrow(results$outcomes$outliers)
  } else {
    cat("\n✓ No outlier outcome rates\n")
  }

  # ---------------------------------------------------------------------------
  print_header(6, "demographics")

  results$demographics = qc_demographics(main_folder, sites,
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

  # ---------------------------------------------------------------------------
  print_header(7, "cross-site heterogeneity")

  results$heterogeneity = qc_heterogeneity(main_folder, sites)

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

  # ---------------------------------------------------------------------------
  print_header(8, "score value ranges")

  results$score_ranges = qc_score_ranges(main_folder, sites, score_ranges)

  if ("message" %in% names(results$score_ranges)) {
    cat("✓", results$score_ranges$message, "\n")
  } else if (nrow(results$score_ranges) > 0) {
    cat("⚠️  Score values outside expected ranges:\n")
    print(results$score_ranges)
    n_issues = n_issues + nrow(results$score_ranges)
  }

  if (check_analysis_outputs) {

    # -------------------------------------------------------------------------
    print_header(9, "threshold analysis outputs")

    results$threshold = qc_threshold_outputs(main_folder, sites, plausibility_bounds)

    if (!is.null(results$threshold$sesp) && nrow(results$threshold$sesp) > 0) {
      cat("⚠️  Sensitivity/specificity flags:\n")
      print(results$threshold$sesp[, .(site, score_name, ca_01,
                                       sensitivity, specificity,
                                       flag_cells, flag_margins)])
      n_issues = n_issues + nrow(results$threshold$sesp)
    } else {
      cat("✓ Sensitivity/specificity values plausible; count identities hold\n")
    }

    if (!is.null(results$threshold$cuminc) && nrow(results$threshold$cuminc) > 0) {
      cat("\n⚠️  Cumulative incidence flags:\n")
      print(results$threshold$cuminc)
      n_issues = n_issues + nrow(results$threshold$cuminc)
    } else {
      cat("✓ Cumulative incidence values bounded within [0, 1]\n")
    }

    if (!is.null(results$threshold$first) && nrow(results$threshold$first) > 0) {
      cat("\n⚠️  Time-to-first-positive flags:\n")
      print(results$threshold$first)
      n_issues = n_issues + nrow(results$threshold$first)
    } else {
      cat("✓ Time-to-first-positive values plausible\n")
    }

    # -------------------------------------------------------------------------
    print_header(10, "sensitivity analysis consistency")

    results$sensitivity = qc_sensitivity_analyses(main_folder, sites)

    if (!is.null(results$sensitivity$flagged) && nrow(results$sensitivity$flagged) > 0) {
      cat("⚠️  Unexpected sensitivity analysis sample sizes:\n")
      print(results$sensitivity$flagged[, .(site, variant, score_name, ca_01,
                                            n_main, n_var, ratio)])
      n_issues = n_issues + nrow(results$sensitivity$flagged)
    } else {
      cat("✓ Sensitivity variant sample sizes consistent",
          "(cf variants equal main exactly)\n")
    }

    # -------------------------------------------------------------------------
    print_header(11, "internal consistency")

    results$consistency = qc_internal_consistency(main_folder, sites,
                                                  plausibility_bounds)

    if (!is.null(results$consistency)) {
      print(results$consistency[, .(site, n_cont, n_flow, n_report,
                                    rate_cat  = fmt_pct(rate_cat),
                                    rate_rpt  = fmt_pct(rate_report),
                                    flag)])
      flagged_cons = results$consistency[flag == TRUE]
      if (nrow(flagged_cons) > 0) {
        cat("\n⚠️  Internal consistency failures (counts or composite rate)\n")
        n_issues = n_issues + nrow(flagged_cons)
      } else {
        cat("\n✓ Encounter counts and composite rates consistent across outputs\n")
      }
    }
  }

  # ---------------------------------------------------------------------------
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

# export flagged items ----------------------------------------------------------

export_flags = function(results, output_dir = here("qc_output")) {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if (!is.null(results$structure$summary)) {
    missing_files = results$structure$summary[exists == FALSE]
    if (nrow(missing_files) > 0) {
      fwrite(missing_files, file.path(output_dir, "flag_missing_files.csv"))
    }
  }

  if (!is.null(results$structure$unexpected) && nrow(results$structure$unexpected) > 0) {
    fwrite(results$structure$unexpected, file.path(output_dir, "flag_unexpected_files.csv"))
  }

  if (!is.null(results$manifest$flags) && nrow(results$manifest$flags) > 0) {
    fwrite(results$manifest$flags, file.path(output_dir, "flag_manifest.csv"))
  }

  if (!is.null(results$flow$outliers) && nrow(results$flow$outliers) > 0) {
    fwrite(results$flow$outliers, file.path(output_dir, "flag_flow_outliers.csv"))
  }

  if (!is.null(results$outcomes$outliers) && nrow(results$outcomes$outliers) > 0) {
    fwrite(results$outcomes$outliers, file.path(output_dir, "flag_outcome_outliers.csv"))
  }

  if (is.data.table(results$score_ranges) &&
      !("message" %in% names(results$score_ranges)) &&
      nrow(results$score_ranges) > 0) {
    fwrite(results$score_ranges, file.path(output_dir, "flag_score_ranges.csv"))
  }

  if (!is.null(results$sensitivity$flagged) && nrow(results$sensitivity$flagged) > 0) {
    fwrite(results$sensitivity$flagged, file.path(output_dir, "flag_sensitivity.csv"))
  }

  if (!is.null(results$consistency)) {
    flagged_cons = results$consistency[flag == TRUE]
    if (nrow(flagged_cons) > 0) {
      fwrite(flagged_cons, file.path(output_dir, "flag_internal_consistency.csv"))
    }
  }

  message("Flags exported to: ", output_dir)
}

# ==============================================================================
# HTML REPORT
# Requires: ggplot2, htmltools, base64enc
# ==============================================================================

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

pal_cancer  = c("0" = "#4575b4", "1" = "#d73027")
pal_outcome = c(
  "composite" = "#313695",
  "nohospice" = "#4575b4",
  "wardicu"   = "#91bfdb",
  "warddeath" = "#d73027",
  "hospicedc" = "#fc8d59",
  "dead"      = "#a50026",
  "hospice"   = "#fdae61",
  "icu"       = "#fee090",
  "va"        = "#74add1",
  "imv"       = "#abd9e9"
)
pal_flag = c("TRUE" = "#d73027", "FALSE" = "#1a9850")

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

  # focus on the composite and its components
  key_outcomes = c("composite", "wardicu", "warddeath", "hospicedc")
  dt = dt[var %in% key_outcomes]

  # nicer labels
  dt[, var_label := fcase(
    var == "composite", "Composite",
    var == "wardicu",   "Ward-to-ICU",
    var == "warddeath", "Ward death",
    var == "hospicedc", "Hospice discharge",
    default = var
  )]

  pooled = outcomes_results$pooled[var %in% key_outcomes]
  pooled[, var_label := fcase(
    var == "composite", "Composite",
    var == "wardicu",   "Ward-to-ICU",
    var == "warddeath", "Ward death",
    var == "hospicedc", "Hospice discharge",
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
      title = "Composite Outcome and Components by Site and Cancer Status",
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
    geom_errorbar(aes(xmin = age_mean - age_sd, xmax = age_mean + age_sd),
                  width = 0.3, alpha = 0.5, orientation = "y",
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
    geom_errorbar(aes(xmin = vw_mean - vw_sd, xmax = vw_mean + vw_sd),
                  width = 0.3, alpha = 0.5, orientation = "y",
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
  n_unexpected    = if (!is.null(qc_results$structure$unexpected)) {
    nrow(qc_results$structure$unexpected)
  } else 0
  n_manifest_flags = if (!is.null(qc_results$manifest$flags)) {
    nrow(qc_results$manifest$flags)
  } else 0
  n_sample_flags = nrow(qc_results$sample_sizes[
    flag_n_low == TRUE | flag_n_high == TRUE | flag_ca_low == TRUE | flag_ca_high == TRUE
  ])
  n_flow_outliers    = if (!is.null(qc_results$flow$outliers)) nrow(qc_results$flow$outliers) else 0
  n_outcome_outliers = if (!is.null(qc_results$outcomes$outliers)) nrow(qc_results$outcomes$outliers) else 0
  n_total_issues     = n_missing_files + n_unexpected + n_manifest_flags +
    n_sample_flags + n_flow_outliers + n_outcome_outliers

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
    tags$div(class = paste("summary-card", if (n_manifest_flags + n_unexpected > 0) "warn" else "ok"),
             tags$div(class = "number", n_manifest_flags + n_unexpected),
             tags$div(class = "label", "Manifest/Stale Flags")),
    tags$div(class = paste("summary-card", if (n_total_issues > 0) "warn" else "ok"),
             tags$div(class = "number", n_total_issues),
             tags$div(class = "label", "Total Flags"))
  )

  # section: manifest reconciliation
  section_manifest = collapsible_section(
    title = "MANIFEST Reconciliation",
    id = "manifest",
    content = tagList(
      status_badge(n_manifest_flags + n_unexpected, "discrepancies"),
      tags$h3("Per-Site Summary"),
      dt_to_html(qc_results$manifest$summary),
      if (n_manifest_flags > 0) tagList(
        tags$h3("Discrepancies"),
        dt_to_html(qc_results$manifest$flags, max_rows = 100)
      ),
      if (n_unexpected > 0) tagList(
        tags$h3("Unexpected Files on Disk"),
        dt_to_html(qc_results$structure$unexpected, max_rows = 100)
      )
    )
  )

  # section: sample sizes
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

  # section: exclusion flow
  cat("  Creating exclusion flow plot...\n")
  p_flow = plot_exclusion_flow(qc_results$flow)

  section_flow = collapsible_section(
    title = "Exclusion Flow",
    id = "exclusion-flow",
    content = tagList(
      status_badge(n_flow_outliers, "flags"),
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

  # section: outcomes
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

  # section: demographics
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

  # section: heterogeneity
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

  # section: score ranges (if issues)
  section_scores = NULL
  if (!is.null(qc_results$score_ranges) &&
      is.data.table(qc_results$score_ranges) &&
      !("message" %in% names(qc_results$score_ranges)) &&
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

  # section: file inventory
  section_files = collapsible_section(
    title = "File Inventory",
    id = "files",
    open = FALSE,
    content = tagList(
      status_badge(n_missing_files, "missing"),
      dt_to_html(qc_results$structure$summary[, .(
        site, file_type, location, filename, exists, n_rows
      )], max_rows = 200)
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
      section_manifest,
      section_sample,
      section_flow,
      section_outcomes,
      section_demo,
      section_het,
      section_scores,
      section_files,
      tags$div(
        class = "footnote",
        tags$p("Generated by qc.R (round two)"),
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

# run --------------------------------------------------------------------------

qc_results = run_qc(
  main_folder            = here(),
  sites_expected         = expected_sites,
  check_analysis_outputs = TRUE
)

export_flags(qc_results)

generate_qc_report(
  qc_results,
  output_file = here("qc_report.html")
)
