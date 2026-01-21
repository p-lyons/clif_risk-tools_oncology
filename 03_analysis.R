scores = read_parquet(here("proj_tables", "scores_full.parquet"))

# constants --------------------------------------------------------------------

THRESHOLDS = tidytable(
  score_name = c("sirs_total", "qsofa_total", "mews_total", "news_total", "mews_sf_total"),
  threshold  = c(2L, 2L, 5L, 5L, 7L)
)

HORIZONS = c(12L, 24L)

VARIANTS = c("main", "se_no_ed_req", "se_fullcode_only", "se_win0_96h", "se_one_enc_per_pt")

# cleaner allowed artifacts list
.allowed = list(
  main        = c("maxscores"),
  threshold   = c("ever", "sesp", "cuminc", "first", "upset"),
  sensitivity = c("maxscores", "counts"),
  horizon     = c("counts"),
  diagnostics = c("overall", "by_cancer", "max_scores"),
  meta        = c("coefficients", "score_sds")
)

# helper functions -------------------------------------------------------------

make_y = function(h_to_event, horizon) {
  as.integer(!is.na(h_to_event) & h_to_event >= 0 & h_to_event <= horizon)
}

#' Build consistent filename for artifacts
#' Pattern: {artifact}[-{strata}][-h{horizon}][-{variant}]-{site}.csv
#' All components except artifact and site are optional
#' 
.build_filename = function(artifact, site, strata = NULL, horizon = NULL, variant = NULL) {
  
  stopifnot(nzchar(artifact), nzchar(site))
  
  parts = artifact
  
  if (!is.null(strata) && nzchar(strata)) {
    parts = paste0(parts, "-", strata)
  }
  if (!is.null(horizon) && !is.na(horizon)) {
    parts = paste0(parts, "-h", horizon)
  }
  if (!is.null(variant) && nzchar(variant)) {
    parts = paste0(parts, "-", variant)
  }
  paste0(parts, "-", site, ".csv")
}

write_artifact = function(df, analysis, artifact, site,
                          strata = NULL, horizon = NULL, variant = NULL,
                          root = "upload_to_box") {
  
  stopifnot(analysis %in% names(.allowed), artifact %in% .allowed[[analysis]])
  
  dir = file.path(root, analysis)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  
  fn   = .build_filename(artifact, site, strata, horizon, variant)
  path = file.path(dir, fn)
  
  message("    Writing: ", path)
  fwrite(df, path)
  invisible(path)
}

sample_one_encounter_per_patient = function(enc_df) {
  dt = as.data.table(enc_df)
  set.seed(2025L)
  dt[, .SD[sample(.N, 1L)], by = patient_id]$joined_hosp_id
}

sample_one_idx = function(dt) {
  dt[, .I[sample(.N, 1L)], by = joined_hosp_id]$V1
}

# variant materializers --------------------------------------------------------

#' Materialize score-level (long) data for a variant
#' Returns data.table with one row per score measurement
materialize_variant_long = function(variant, scores_long_base, scores, fc) {
  
  dt = copy(scores_long_base)
  
  switch(
    variant,
    main = {
      dt = dt[ed_admit_01 == 1L]
    },
    se_no_ed_req = {
      dt = dt
    },
    se_fullcode_only = {
      dt = dt[ed_admit_01 == 1L & fullcode_01 == 1L]
    },
    se_win0_96h = {
      dt = dt[ed_admit_01 == 1L & h_from_admit >= 0 & h_from_admit <= 96]
    },
    se_one_enc_per_pt = {
      enc_tbl   = funique(select(scores, joined_hosp_id, patient_id, ed_admit_01))
      enc_tbl   = fsubset(enc_tbl, ed_admit_01 == 1L)
      keep_encs = sample_one_encounter_per_patient(enc_tbl)
      dt        = dt[joined_hosp_id %in% keep_encs]
    },
    stop("Unknown variant: ", variant)
  )
  
  dt[]
}

#' Materialize encounter-level max scores for a variant
#' Returns data.table with one row per encounter × score
materialize_variant_max = function(variant, scores, cohort, fc) {
  
  # start with appropriate subset of scores
  dt_scores = switch(
    variant,
    main = {
      fsubset(scores, ed_admit_01 == 1L)
    },
    se_no_ed_req = {
      copy(scores)
    },
    se_fullcode_only = {
      fc_encs = fsubset(cohort, tolower(initial_code_status) == "full")$joined_hosp_id
      fsubset(scores, ed_admit_01 == 1L & joined_hosp_id %in% fc_encs)
    },
    se_win0_96h = {
      # KEY FIX: filter to 0-96h BEFORE taking max
      fsubset(scores, ed_admit_01 == 1L & h_from_admit >= 0 & h_from_admit <= 96)
    },
    se_one_enc_per_pt = {
      enc_tbl   = funique(select(scores, joined_hosp_id, patient_id, ed_admit_01))
      enc_tbl   = fsubset(enc_tbl, ed_admit_01 == 1L)
      keep_encs = sample_one_encounter_per_patient(enc_tbl)
      fsubset(scores, joined_hosp_id %in% keep_encs)
    },
    stop("Unknown variant: ", variant)
  )
  
  if (nrow(dt_scores) == 0) {
    warning("No data for variant: ", variant)
    return(NULL)
  }
  
  # calculate max scores per encounter
  dt_max = fgroup_by(dt_scores, joined_hosp_id) |>
    fsummarize(
      patient_id  = ffirst(patient_id),
      hospital_id = ffirst(hospital_id),
      ca_01       = ffirst(ca_01),
      ed_admit_01 = ffirst(ed_admit_01),
      fullcode_01 = fifelse(ffirst(joined_hosp_id) %in% fc, 1L, 0L),
      sirs_max    = fmax(sirs_total,    na.rm = TRUE),
      qsofa_max   = fmax(qsofa_total,   na.rm = TRUE),
      mews_max    = fmax(mews_total,    na.rm = TRUE),
      news_max    = fmax(news_total,    na.rm = TRUE),
      mews_sf_max = fmax(mews_sf_total, na.rm = TRUE),
      outcome     = ffirst(o_primary_01)
    ) |>
    pivot_longer(
      cols      = ends_with("_max"),
      names_to  = "score_name",
      values_to = "max_value"
    ) |>
    ftransform(score_name = str_remove(score_name, "_max"))
  
  # handle -Inf from fmax on empty sets
  dt_max[is.infinite(max_value), max_value := NA_integer_]
  
  as.data.table(dt_max)[]
}

aggregate_maxscores = function(dt, site_lowercase) {
  
  dt = as.data.table(dt)
  need = c("score_name", "ca_01", "max_value", "outcome")
  miss = setdiff(need, names(dt))
  if (length(miss)) stop("aggregate_maxscores(): missing cols: ", paste(miss, collapse = ", "))
  
  dt[, .(n = .N), by = .(score_name, ca_01, max_value, outcome)
  ][, site := site_lowercase][]
}

run_horizon_counts = function(dt, horizons, site_lowercase) {
  
  counts_list = lapply(horizons, function(HH) {
    dt_copy = copy(dt)
    dt_copy[, outcome := make_y(h_to_event, HH)]
    dt_copy[, .(n = .N), by = .(score_name, ca_01, value, outcome)
    ][, `:=`(site = site_lowercase, h = HH)]
  })
  
  rbindlist(counts_list, use.names = TRUE)[]
}

run_horizon_counts_bootstrap = function(dt, horizons, site_lowercase, B = 400L) {
  
  set.seed(2025L)
  boot_list = vector("list", B)
  base_dt   = dt[, .(joined_hosp_id, score_name, ca_01, value, h_to_event)]
  setkey(base_dt, joined_hosp_id)
  
  for (b in seq_len(B)) {
    idx  = sample_one_idx(base_dt)
    samp = base_dt[idx]
    
    counts_b = rbindlist(lapply(horizons, function(HH) {
      samp[, outcome := make_y(h_to_event, HH)]
      samp[, .(n = .N), by = .(score_name, ca_01, value, outcome)
      ][, `:=`(site = site_lowercase, h = HH, iter = b)]
    }), use.names = TRUE)
    
    boot_list[[b]] = counts_b
    if (b %% 25 == 0) message("    Bootstrap iteration ", b, "/", B)
  }
  
  rbindlist(boot_list, use.names = TRUE)[]
}

fit_one_score = function(df, score, outcome_col = "outcome", cancer_col = "ca_01",
                         hosp_col = "hospital_id", score_name_col = "score_name",
                         max_value_col = "max_value") {
  
  df_sub = df[df[[score_name_col]] == score, , drop = FALSE]
  if (nrow(df_sub) == 0L) stop("No rows for score = ", score)
  
  has_multi_hosp = fnunique(df_sub[[hosp_col]]) > 1
  
  f = if (has_multi_hosp) {
    as.formula(paste0(outcome_col, " ~ ", max_value_col, " * ", cancer_col, " + (1|", hosp_col, ")"))
  } else {
    as.formula(paste0(outcome_col, " ~ ", max_value_col, " * ", cancer_col))
  }
  
  fit = tryCatch(
    glmmTMB(f, family = binomial(), data = df_sub),
    warning = function(w) {
      if (grepl("non-positive-definite", conditionMessage(w))) {
        glmmTMB(
          f,
          family  = binomial(),
          data    = df_sub,
          control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"))
        )
      } else {
        warning(w)
        glmmTMB(f, family = binomial(), data = df_sub)
      }
    }
  )
  
  co = summary(fit)$coefficients$cond
  
  int_rows = grep(paste0("^", max_value_col, ":"), rownames(co), value = TRUE)
  int_rows = int_rows[grepl(paste0("^", max_value_col, ":", cancer_col), int_rows)]
  
  if (length(int_rows) != 1L) {
    stop("Could not uniquely identify interaction term. Found: ", paste(int_rows, collapse = ", "))
  }
  
  est = co[int_rows, "Estimate"]
  se  = co[int_rows, "Std. Error"]
  
  tidytable(
    score       = score,
    beta_int    = est,
    se_int      = se,
    site_n      = nobs(fit),
    n_events    = fsum(df_sub[[outcome_col]] == 1L, na.rm = TRUE),
    n_hospitals = fnunique(df_sub[[hosp_col]]),
    converged   = isTRUE(fit$fit$convergence == 0),
    hess_pd     = !any(is.na(se))
  )
}

# ==============================================================================
# PREPARE BASE DATASETS
# ==============================================================================

message("\n== Preparing base datasets ==")

jp = select(cohort, patient_id, joined_hosp_id, hospital_id)
fc = fsubset(cohort, tolower(initial_code_status) == "full")$joined_hosp_id

# add patient/hospital info to scores
scores = join(scores, jp, how = "inner", multiple = FALSE) |>
  ftransform(h_from_admit = as.numeric(difftime(time, in_dttm, units = "hours")))

# long format for point-level analyses
scores_long_base =
  select(scores, joined_hosp_id, ends_with("01"), ends_with("dttm"), time, ends_with("total")) |>
  select(-outcome_nohospc_dttm, -end_dttm, -out_dttm) |>
  pivot_longer(
    cols      = ends_with("total"),
    names_to  = "score_name",
    values_to = "value"
  ) |>
  ftransform(h_to_event   = as.numeric(difftime(outcome_dttm, time,    units = "hours"))) |>
  ftransform(h_from_admit = as.numeric(difftime(time,         in_dttm, units = "hours")))

scores_long_base = join(scores_long_base, jp, how = "inner", multiple = FALSE) |>
  ftransform(fullcode_01 = fifelse(joined_hosp_id %in% fc, 1L, 0L)) |>
  select(ends_with("id"), ends_with("01"), score_name, value, starts_with("h_"))

setDT(scores_long_base)
setkey(scores_long_base, joined_hosp_id, patient_id)

message("  scores_long_base: ", format(nrow(scores_long_base), big.mark = ","), " rows")

# ==============================================================================
# THRESHOLD ANALYSES (main cohort only)
# ==============================================================================

message("\n== Threshold analyses ==")

# --- ever positive ------------------------------------------------------------

ever_positive =
  fsubset(scores, ed_admit_01 == 1) |>
  pivot_longer(
    cols      = ends_with("total"),
    names_to  = "score_name",
    values_to = "value"
  ) |>
  join(THRESHOLDS, how = "inner", multiple = FALSE) |>
  ftransform(positive = value >= threshold) |>
  fsubset(positive == TRUE) |>
  roworder(time) |>
  fgroup_by(joined_hosp_id, score_name) |>
  fsummarize(
    time_to_positive_h = ffirst(h_from_admit),
    first_positive_val = ffirst(value)
  )

all_encs =
  fsubset(cohort, ed_admit_01 == 1) |>
  ftransform(deadhospice_01 = pmin(dead_01 + hospice_01, 1L)) |>
  select(joined_hosp_id, ca_01, deadhospice_01)

ever_positive_complete =
  tidyr::expand_grid(
    joined_hosp_id = all_encs$joined_hosp_id,
    score_name     = THRESHOLDS$score_name
  ) |>
  as_tidytable() |>
  join(all_encs,      how = "left", multiple = FALSE) |>
  join(ever_positive, how = "left", multiple = FALSE) |>
  ftransform(ever_positive = as.integer(!is.na(time_to_positive_h))) |>
  select(joined_hosp_id, score_name, ca_01, ever_positive, deadhospice_01)

# export: ever positive counts
ever_positive_agg = as.data.table(ever_positive_complete)[
  , .(n = .N), by = .(score_name, ca_01, ever_positive, deadhospice_01)
][, site := site_lowercase][]

write_artifact(
  df       = ever_positive_agg,
  analysis = "threshold",
  artifact = "ever",
  site     = site_lowercase,
  strata   = "ca"
)

# export: sensitivity/specificity
ever_positive_sesp = as.data.table(ever_positive_complete)[
  , .(
    n_total      = .N,
    n_outcome    = sum(deadhospice_01 == 1L, na.rm = TRUE),
    n_no_outcome = sum(deadhospice_01 == 0L, na.rm = TRUE),
    n_pos        = sum(ever_positive == 1L, na.rm = TRUE),
    n_neg        = sum(ever_positive == 0L, na.rm = TRUE),
    tp           = sum(ever_positive == 1L & deadhospice_01 == 1L, na.rm = TRUE),
    fp           = sum(ever_positive == 1L & deadhospice_01 == 0L, na.rm = TRUE),
    tn           = sum(ever_positive == 0L & deadhospice_01 == 0L, na.rm = TRUE),
    fn           = sum(ever_positive == 0L & deadhospice_01 == 1L, na.rm = TRUE)
  ),
  by = .(score_name, ca_01)
][, `:=`(
  sensitivity = tp / (tp + fn),
  specificity = tn / (tn + fp),
  ppv         = tp / (tp + fp),
  npv         = tn / (tn + fn),
  site        = site_lowercase
)][]

write_artifact(
  df       = ever_positive_sesp,
  analysis = "threshold",
  artifact = "sesp",
  site     = site_lowercase,
  strata   = "ca"
)

# --- first positive (component analysis) --------------------------------------

message("  Computing first positive components...")

extract_first_positive = function(scores_dt, score_prefix, threshold_val, score_total_col) {
  
  scores_dt |>
    fsubset(get(score_total_col) >= threshold_val & ed_admit_01 == 1) |>
    select(joined_hosp_id, ca_01, o_primary_01, h_from_admit,
            matches(paste0("^", score_prefix, "_"))) |>
    select(-matches("_total$")) |>
    roworder(h_from_admit) |>
    fgroup_by(joined_hosp_id) |>
    ffirst() |>
    pivot_longer(
      cols         = matches(paste0("^", score_prefix, "_")),
      names_to     = "component",
      values_to    = "value",
      names_prefix = paste0(score_prefix, "_")
    ) |>
    ftransform(score = score_prefix) |>
    fsubset(value > 0)
}

first_sirs   = extract_first_positive(scores, "sirs",  2L, "sirs_total")
first_qsofa  = extract_first_positive(scores, "qsofa", 2L, "qsofa_total")
first_mews   = extract_first_positive(scores, "mews",  5L, "mews_total")
first_news   = extract_first_positive(scores, "news",  7L, "news_total")

# mews_sf needs special handling for the sf column
first_mewssf =
  fsubset(scores, mews_sf_total >= 7L & ed_admit_01 == 1) |>
  select(joined_hosp_id, ca_01, o_primary_01, h_from_admit, matches("^mews_"), sf) |>
  select(-matches("_total$")) |>
  roworder(h_from_admit) |>
  fgroup_by(joined_hosp_id) |>
  ffirst() |>
  rename(mews_sf = sf) |>
  pivot_longer(
    cols         = matches("^mews_"),
    names_to     = "component",
    values_to    = "value",
    names_prefix = "mews_"
  ) |>
  ftransform(score = "mews_sf") |>
  fsubset(value > 0)

firsts = rowbind(first_sirs, first_qsofa, first_mews, first_news, first_mewssf)

# export: time to first positive
time_to_positive =
  select(firsts, joined_hosp_id, score, ca_01, o_primary_01, h_from_admit) |>
  funique() |>
  ftransform(pos_day0 = fifelse(h_from_admit <= 24, 1L, 0L)) |>
  ftransform(pos_day1 = fifelse(h_from_admit > 24 & h_from_admit <= 48, 1L, 0L)) |>
  fgroup_by(score, ca_01, o_primary_01) |>
  fsummarize(
    n                  = fnobs(joined_hosp_id),
    h_from_admit_sum   = fsum(h_from_admit),
    h_from_admit_sumsq = fsum(h_from_admit^2),
    n_pos_day0         = fsum(pos_day0),
    n_pos_day1         = fsum(pos_day1)
  ) |>
  ftransform(site = site_lowercase)

write_artifact(
  df       = time_to_positive,
  analysis = "threshold",
  artifact = "first",
  site     = site_lowercase,
  strata   = "ca"
)

# export: component upset data
n_positive =
  select(firsts, joined_hosp_id, ca_01, o_primary_01, score) |>
  funique() |>
  fgroup_by(ca_01, o_primary_01, score) |>
  fsummarize(n_encs = GRPN())

upset_components =
  select(firsts, joined_hosp_id, ca_01, o_primary_01, score, component) |>
  funique() |>
  join(n_positive, how = "left", multiple = TRUE) |>
  fgroup_by(ca_01, o_primary_01, score, component) |>
  fsummarize(n = GRPN(), n_encs = fmax(n_encs)) |>
  ftransform(site = site_lowercase)

write_artifact(
  df       = upset_components,
  analysis = "threshold",
  artifact = "upset",
  site     = site_lowercase,
  strata   = "components"
)

# --- cumulative incidence of score positivity ---------------------------------

message("  Computing cumulative incidence...")

enc_windows =
  fsubset(scores, ed_admit_01 == 1) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(
    ca_01              = ffirst(ca_01),
    o_primary_01       = ffirst(o_primary_01),
    h_admit_to_outcome = ffirst(as.numeric(difftime(outcome_dttm, in_dttm, units = "hours"))),
    h_admit_to_dc      = ffirst(as.numeric(difftime(out_dttm, in_dttm, units = "hours"))),
    max_obs_time       = fmax(h_from_admit, na.rm = TRUE)
  ) |>
  ftransform(h_followup_end = pmin(h_admit_to_outcome, h_admit_to_dc, max_obs_time, na.rm = TRUE))

time_breaks = seq(0, 168, by = 8)

first_pos_times =
  select(firsts, joined_hosp_id, score, h_from_admit) |>
  funique() |>
  rename(h_first_positive = h_from_admit)

cuminc_data =
  tidyr::expand_grid(
    joined_hosp_id = enc_windows$joined_hosp_id,
    score          = unique(firsts$score),
    time_bin_start = time_breaks
  ) |>
  as_tidytable() |>
  join(enc_windows, how = "left", multiple = FALSE) |>
  join(first_pos_times, how = "left", multiple = FALSE) |>
  ftransform(
    at_risk         = time_bin_start < h_followup_end,
    became_positive = !is.na(h_first_positive) & h_first_positive <= time_bin_start,
    positive_in_bin = !is.na(h_first_positive) &
      h_first_positive >= time_bin_start &
      h_first_positive < time_bin_start + 8
  )

cuminc_summary =
  fsubset(cuminc_data, at_risk == TRUE) |>
  fgroup_by(score, ca_01, o_primary_01, time_bin_start) |>
  fsummarize(
    n_at_risk    = GRPN(),
    n_became_pos = fsum(became_positive, na.rm = TRUE),
    n_pos_in_bin = fsum(positive_in_bin, na.rm = TRUE),
    cum_inc      = fsum(became_positive, na.rm = TRUE) / GRPN()
  ) |>
  ftransform(site = site_lowercase)

write_artifact(
  df       = cuminc_summary,
  analysis = "threshold",
  artifact = "cuminc",
  site     = site_lowercase,
  strata   = "ca"
)

# ==============================================================================
# VARIANT ANALYSES
# ==============================================================================

message("\n== Variant analyses ==")

variant_diagnostics = list()

for (v in VARIANTS) {
  
  message("\n--- Variant: ", v, " ---")
  
  # --- point-level (long) data ------------------------------------------------
  
  dt_long = materialize_variant_long(v, scores_long_base, scores, fc)
  
  diag_overall = dt_long[, .(
    n_rows     = .N,
    n_enc      = uniqueN(joined_hosp_id),
    n_pat      = uniqueN(patient_id),
    n_enc_ca   = uniqueN(joined_hosp_id[ca_01 == 1L]),
    n_enc_noca = uniqueN(joined_hosp_id[ca_01 == 0L])
  )]
  
  diag_by_cancer = dt_long[, .SD[1], by = joined_hosp_id][, .(
    n_enc        = .N,
    n_events_12h = sum(make_y(h_to_event, 12L) == 1L, na.rm = TRUE),
    n_events_24h = sum(make_y(h_to_event, 24L) == 1L, na.rm = TRUE),
    evt_rate_12h = mean(make_y(h_to_event, 12L), na.rm = TRUE),
    evt_rate_24h = mean(make_y(h_to_event, 24L), na.rm = TRUE)
  ), by = ca_01]
  
  message("  Encounters: ", format(diag_overall$n_enc, big.mark = ","),
          " (CA: ", format(diag_overall$n_enc_ca, big.mark = ","),
          ", No CA: ", format(diag_overall$n_enc_noca, big.mark = ","), ")")
  
  # horizon counts
  message("  Computing horizon counts...")
  counts = run_horizon_counts(dt_long, HORIZONS, site_lowercase)
  
  for (HH in HORIZONS) {
    write_artifact(
      df       = counts[h == HH],
      analysis = if (v == "main") "horizon" else "sensitivity",
      artifact = "counts",
      site     = site_lowercase,
      strata   = "ca",
      horizon  = HH,
      variant  = if (v == "main") NULL else v
    )
  }
  
  # bootstrap counts
  message("  Computing bootstrap counts...")
  counts_boot = run_horizon_counts_bootstrap(dt_long, HORIZONS, site_lowercase, B = 400L)
  
  for (HH in HORIZONS) {
    write_artifact(
      df       = counts_boot[h == HH],
      analysis = if (v == "main") "horizon" else "sensitivity",
      artifact = "counts",
      site     = site_lowercase,
      strata   = "ca",
      horizon  = HH,
      variant  = if (v == "main") "boot" else paste0(v, "-boot")
    )
  }
  
  # --- encounter-level max scores ---------------------------------------------
  
  message("  Computing max scores...")
  dt_max = materialize_variant_max(v, scores, cohort, fc)
  
  if (!is.null(dt_max) && nrow(dt_max) > 0) {
    
    diag_max = dt_max[, .(
      n_enc      = .N / 5L,  # 5 score types per encounter
      n_events   = sum(outcome == 1L, na.rm = TRUE) / 5L,
      event_rate = mean(outcome, na.rm = TRUE)
    ), by = ca_01]
    
    message("  Max scores: ", format(diag_max[1, n_enc], big.mark = ","), " encounters")
    
    dt_max_agg = aggregate_maxscores(dt_max, site_lowercase)
    
    write_artifact(
      df       = dt_max_agg,
      analysis = if (v == "main") "main" else "sensitivity",
      artifact = "maxscores",
      site     = site_lowercase,
      strata   = "ca",
      variant  = if (v == "main") NULL else v
    )
    
    variant_diagnostics[[v]] = list(
      overall   = diag_overall,
      by_cancer = diag_by_cancer,
      max_scores = diag_max
    )
    
  } else {
    message("  No max scores for this variant")
    variant_diagnostics[[v]] = list(
      overall   = diag_overall,
      by_cancer = diag_by_cancer
    )
  }
}

# export diagnostics
message("\n== Writing diagnostics ==")

diag_overall_all = rbindlist(
  lapply(names(variant_diagnostics), function(v) {
    cbind(variant = v, variant_diagnostics[[v]]$overall)
  }),
  use.names = TRUE, fill = TRUE
)[, site := site_lowercase]

diag_by_cancer_all = rbindlist(
  lapply(names(variant_diagnostics), function(v) {
    cbind(variant = v, variant_diagnostics[[v]]$by_cancer)
  }),
  use.names = TRUE, fill = TRUE
)[, site := site_lowercase]

diag_max_all = rbindlist(
  lapply(names(variant_diagnostics), function(v) {
    if (!is.null(variant_diagnostics[[v]]$max_scores)) {
      cbind(variant = v, variant_diagnostics[[v]]$max_scores)
    }
  }),
  use.names = TRUE, fill = TRUE
)
if (nrow(diag_max_all) > 0) diag_max_all[, site := site_lowercase]

write_artifact(df = diag_overall_all,   analysis = "diagnostics", artifact = "overall",    site = site_lowercase)
write_artifact(df = diag_by_cancer_all, analysis = "diagnostics", artifact = "by_cancer",  site = site_lowercase)
if (nrow(diag_max_all) > 0) {
  write_artifact(df = diag_max_all, analysis = "diagnostics", artifact = "max_scores", site = site_lowercase)
}

# ==============================================================================
# SUBGROUP: CANCER TYPE (LIQUID VS SOLID)
# ==============================================================================

message("\n== Cancer type subgroup analysis ==")

liq = fsubset(cohort, liquid_01 == 1)$joined_hosp_id

# max scores by liquid/solid
dt_max_main = materialize_variant_max("main", scores, cohort, fc)

dt_max_liquid =
  fsubset(dt_max_main, ca_01 == 1) |>
  ftransform(liquid_01 = fifelse(joined_hosp_id %in% liq, 1L, 0L))

dt_max_liquid_agg = dt_max_liquid[
  , .(n = .N), by = .(score_name, liquid_01, max_value, outcome)
][, site := site_lowercase][]

write_artifact(
  df       = dt_max_liquid_agg,
  analysis = "main",
  artifact = "maxscores",
  site     = site_lowercase,
  strata   = "liquid"
)

# horizon counts by liquid/solid
scores_long_cancer = fsubset(scores_long_base, ca_01 == 1 & ed_admit_01 == 1) |>
  ftransform(liquid_01 = fifelse(joined_hosp_id %in% liq, 1L, 0L))

setDT(scores_long_cancer)
scores_long_cancer[, outcome := make_y(h_to_event, 24L)]

counts_liquid_24h = scores_long_cancer[
  , .(n = .N), by = .(score_name, liquid_01, value, outcome)
][, `:=`(site = site_lowercase, h = 24L)]

write_artifact(
  df       = counts_liquid_24h,
  analysis = "horizon",
  artifact = "counts",
  site     = site_lowercase,
  strata   = "liquid",
  horizon  = 24L
)

# ==============================================================================
# META-ANALYSIS INPUTS
# ==============================================================================

message("\n== Meta-analysis inputs ==")

library(glmmTMB)

df_model =
  fsubset(dt_max_main, ed_admit_01 == 1) |>
  fgroup_by(joined_hosp_id, score_name) |>
  fsummarize(
    hospital_id = ffirst(hospital_id),
    ca_01       = fmax(ca_01),
    max_value   = fmax(max_value),
    outcome     = fmax(outcome)
  )

scores_to_run = funique(df_model$score_name)

site_fit_tbl = rbindlist(
  lapply(scores_to_run, function(sc) {
    message("  Fitting: ", sc)
    fit_one_score(df_model, score = sc)
  }),
  use.names = TRUE, fill = TRUE
)[, site := site_lowercase][]

write_artifact(
  df       = site_fit_tbl,
  analysis = "meta",
  artifact = "coefficients",
  site     = site_lowercase
)

score_sds =
  fgroup_by(df_model, score_name) |>
  fsummarize(
    sd_score     = fsd(max_value, na.rm = TRUE),
    mean_score   = fmean(max_value, na.rm = TRUE),
    n_encounters = fnobs(max_value)
  ) |>
  rename(score = score_name) |>
  ftransform(site = site_lowercase)

write_artifact(
  df       = score_sds,
  analysis = "meta",
  artifact = "score_sds",
  site     = site_lowercase
)

# ==============================================================================
# UPSET PLOT DATA (encounter-level)
# ==============================================================================

message("\n== Upset plot data ==")

THRESHOLDS_short = copy(THRESHOLDS)
THRESHOLDS_short[, score_name := str_remove(score_name, "_total")]

upset_enc =
  fsubset(dt_max_main, ed_admit_01 == 1) |>
  join(THRESHOLDS_short, how = "left", multiple = TRUE) |>
  ftransform(positive = fifelse(max_value >= threshold, 1L, 0L)) |>
  select(joined_hosp_id, ca_01, outcome, score_name, positive) |>
  pivot_wider(names_from = score_name, values_from = positive)

upset_counts = upset_enc[
  , .(n = .N), by = .(ca_01, outcome, sirs, qsofa, mews, news, mews_sf)
][, site := site_lowercase]

setorder(upset_counts, -n)

write_artifact(
  df       = upset_counts,
  analysis = "threshold",
  artifact = "upset",
  site     = site_lowercase,
  strata   = "ca"
)

# ==============================================================================
# CLEANUP
# ==============================================================================

message("\n== Complete ==")

rm(dt_long, dt_max, dt_max_main, dt_max_liquid, counts, counts_boot)
rm(scores_long_base, scores_long_cancer, cuminc_data)
gc()

message("Files written to: ", normalizePath(file.path("upload_to_box")))
