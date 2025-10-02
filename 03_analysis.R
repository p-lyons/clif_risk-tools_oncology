# constants --------------------------------------------------------------------

THRESHOLDS = tidytable(
  score_name = c("sirs_total", "qsofa_total", "mews_total", "news_total"),
  threshold  = c(2L, 2L, 5L, 5L)
)

HORIZONS = c(12L, 24L)

VARIANTS = c("main", "se_no_ed_req", "se_fullcode_only", "se_win0_120h", "se_one_enc_per_pt")

.allowed = list(
  main        = c("maxscores"),
  threshold   = c("ever", "upset"),
  sensitivity = c("maxscores", "counts"),
  horizon     = c("counts")
)

# helper functions -------------------------------------------------------------

## make outcome for horizon analysis -------------------------------------------

make_y = function(h_to_event, horizon) {
  as.integer(!is.na(h_to_event) & h_to_event >= 0 & h_to_event <= horizon)
}

## function to write files -----------------------------------------------------

.build_filename = function(artifact, site, strata = NULL, horizon = NULL, variant = NULL) {
  stopifnot(nzchar(artifact), nzchar(site))
  nm = artifact
  if (!is.null(strata)  && nzchar(strata))  nm = paste0(nm, "_", strata)
  if (!is.null(horizon) && nzchar(horizon)) nm = paste0(nm, "_h", horizon)
  if (!is.null(variant) && nzchar(variant)) nm = paste0(nm, "-", variant)
  paste0(nm, "-", site, ".csv")
}

write_artifact = function(df, analysis, artifact, site,
                          strata = NULL, horizon = NULL, variant = NULL,
                          root = "proj_output") {
  stopifnot(analysis %in% names(.allowed), artifact %in% .allowed[[analysis]])
  dir = file.path(root, analysis)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  fn   = .build_filename(artifact, site, strata, horizon, variant)
  path = file.path(dir, fn)
  tidytable::fwrite(df, path)
  path
}

## deterministic one-encounter-per-patient sampler -----------------------------

sample_one_encounter_per_patient = function(enc_df) {
  dt = as.data.table(enc_df)
  set.seed(2025L)
  dt[, .SD[sample(.N, 1L)], by = patient_id]$joined_hosp_id
}

## sample one observation per encounter (for bootstrapping) --------------------

sample_one_idx = function(dt) {
  dt[, .I[sample(.N, 1L)], by = joined_hosp_id]$V1
}

## materialize a variant view of scores_long_base ------------------------------

materialize_variant = function(variant) {
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
      dt = dt[fullcode_01 == 1L]
    },
    se_win0_120h = {
      dt = dt[h_from_admit >= 0 & h_from_admit <= 120]
    },
    se_one_enc_per_pt = {
      enc_tbl   = funique(fselect(scores, joined_hosp_id, patient_id))
      keep_encs = sample_one_encounter_per_patient(enc_tbl)
      dt        = dt[joined_hosp_id %in% keep_encs]
    },
    stop("Unknown variant: ", variant)
  )
  
  dt[]
}

## materialize variant view of encounter-level max scores ----------------------

materialize_variant_max = function(variant, scores_max_enc) {
  dt = copy(scores_max_enc)
  
  switch(
    variant,
    main = {
      dt = dt[ed_admit_01 == 1L]
    },
    se_no_ed_req = {
      dt = dt
    },
    se_fullcode_only = {
      dt = dt[fullcode_01 == 1L]
    },
    se_win0_120h = {
      return(NULL)
    },
    se_one_enc_per_pt = {
      enc_tbl   = dt[, .(patient_id, joined_hosp_id)]
      keep_encs = sample_one_encounter_per_patient(enc_tbl)
      dt        = dt[joined_hosp_id %in% keep_encs]
    },
    stop("Unknown variant: ", variant)
  )
  
  # remove identifiers before returning
  dt = fselect(dt, -patient_id, -joined_hosp_id)
  as.data.table(dt)[]
}

## aggregate max scores to counts ----------------------------------------------

aggregate_maxscores = function(dt, site_lowercase) {
  dt   = as.data.table(dt)
  need = c("score_name","ca_01","max_value","outcome")
  miss = setdiff(need, names(dt))
  if (length(miss)) stop("aggregate_maxscores(): missing cols: ", paste(miss, collapse=", "))
  dt[, .(n = .N), by = .(score_name, ca_01, max_value, outcome)
  ][, site := site_lowercase][]
}

## collapse small cells --------------------------------------------------------

collapse_small = function(dt, grpvars, val_var = "value", n_var = "n", thresh = 5, max_pct_collapse = 0.01) {
  
  dt          = as.data.table(dt)
  total_n     = sum(dt[[n_var]], na.rm = TRUE)
  collapsed_n = 0L
  
  dt_out = dt[, {
    setorderv(.SD, val_var, -1L)
    vals = .SD[[val_var]]
    ns   = .SD[[n_var]]
    keep = rep(TRUE, length(ns))
    
    for (i in seq_along(ns)) {
      if (ns[i] < thresh && i < length(ns)) {
        collapsed_n <<- collapsed_n + ns[i]
        ns[i + 1] = ns[i + 1] + ns[i]
        keep[i]   = FALSE
      }
    }
    if (length(ns) > 1L && ns[length(ns)] < thresh) {
      collapsed_n <<- collapsed_n + ns[length(ns)]
      ns[length(ns) - 1L] = ns[length(ns) - 1L] + ns[length(ns)]
      keep[length(ns)]    = FALSE
    }
    .(value = vals[keep], n = ns[keep])
  }, by = grpvars]
  
  setnames(dt_out, c("value", "n"), c(val_var, n_var))
  pct = collapsed_n / total_n
  
  if (pct > max_pct_collapse) {
    stop(sprintf("❌ Collapsed %.4f%% of observations (>%.0f%% allowed). Review needed.", 
                 100*pct, 100*max_pct_collapse))
  } else {
    message(sprintf("✅ Collapsed safely: %.4f%% of observations collapsed (≤%.0f%% allowed).",
                    100*pct, 100*max_pct_collapse))
  }
  dt_out
}

## horizon counts (non-bootstrapped) -------------------------------------------

run_horizon_counts = function(dt, horizons, site_lowercase) {
  counts_list = lapply(horizons, function(HH) {
    dt[, outcome := make_y(h_to_event, HH)]
    dt[, .(n = .N), by = .(score_name, ca_01, value, outcome)
    ][, `:=`(site = site_lowercase, h = HH)]
  })
  
  out = rbindlist(counts_list, use.names = TRUE)
  
  if (exists("collapse_small")) {
    out = collapse_small(out, grpvars = c("h", "score_name", "ca_01", "outcome"))
  }
  out[]
}

## horizon counts (clustered bootstrapping) ------------------------------------

run_horizon_counts_bootstrap = function(dt, horizons, site_lowercase, B = 100L) {
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
    if (b %% 10 == 0) message("  Bootstrap iteration ", b, "/", B)
  }
  rbindlist(boot_list, use.names = TRUE)[]
}

# create long-form base with precomputed time features -------------------------

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

jp = select(cohort, patient_id, joined_hosp_id)
fc = fsubset(cohort, tolower(initial_code_status) == "full") |> pull(joined_hosp_id)

scores_long_base = 
  join(scores_long_base, jp, how = "left", multiple = F) |>
  ftransform(fullcode_01 = if_else(joined_hosp_id %in% fc, 1L, 0L)) |>
  select(ends_with("id"), ends_with("01"), score_name, value, starts_with("h_"))

setDTthreads(0)
setDT(scores_long_base)
setkey(scores_long_base, joined_hosp_id, patient_id)

# encounter-level maximum scores (for main AUROC analysis) --------------------

## add patient_id and h_from_admit to scores for downstream use ----------------

scores = 
  join(scores, jp, how = "left", multiple = FALSE) |>
  ftransform(h_from_admit = as.numeric(difftime(time, in_dttm, units = "hours")))

## create encounter-level max scores -------------------------------------------

scores_max_enc = 
  scores |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(
    patient_id  = ffirst(patient_id),
    ca_01       = ffirst(ca_01),
    ed_admit_01 = ffirst(ed_admit_01),
    fullcode_01 = if_else(ffirst(joined_hosp_id) %in% fc, 1L, 0L),
    sirs_max    = fmax(sirs_total,  na.rm = TRUE),
    qsofa_max   = fmax(qsofa_total, na.rm = TRUE),
    mews_max    = fmax(mews_total,  na.rm = TRUE),
    news_max    = fmax(news_total,  na.rm = TRUE),
    outcome     = ffirst(o_primary_01)
  ) |>
  pivot_longer(
    cols      = ends_with("max"),
    names_to  = "score_name",
    values_to = "max_value"
  ) |>
  ftransform(score_name = str_remove(score_name, "_max"))

# ever positive analysis (time to first threshold crossing) -------------------

ever_positive = 
  scores |>
  fsubset(ed_admit_01 == 1) |>
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

## add back encounters that never went positive --------------------------------

all_encs = 
  fsubset(cohort, ed_admit_01 == 1) |>
  fselect(joined_hosp_id, ca_01)

ever_positive_complete =
  tidyr::expand_grid(
    joined_hosp_id = all_encs$joined_hosp_id,
    score_name     = THRESHOLDS$score_name
  ) |>
  tidytable::as_tidytable() |>
  # explicit keys prevent ca_01.x / ca_01.y
  join(all_encs,      how = "left", multiple = FALSE) |>
  join(ever_positive, how = "left", multiple = FALSE) |>
  ftransform(ever_positive = as.integer(!is.na(time_to_positive_h))) |>
  fselect(joined_hosp_id, score_name, ca_01, ever_positive)

# sanity check
stopifnot(all(c("score_name","ca_01","ever_positive") %chin% names(ever_positive_complete)))

## aggregate to counts (no identifiers) ----------------------------------------

ever_positive_agg =
  as.data.table(ever_positive_complete)[
    , .(n = .N), by = .(score_name, ca_01, ever_positive)
  ][, site := site_lowercase][]

write_artifact(
  df       = ever_positive_agg,
  analysis = "threshold",
  artifact = "ever",
  site     = site_lowercase,
  strata   = "ca"
)

# run analyses across variants -------------------------------------------------

for (v in VARIANTS) {
  message("\n== Variant: ", v, " ==")
  
  ## time-varying analyses -----------------------------------------------------
  
  dt_v = materialize_variant(v)
  diag_tbl = dt_v[, .(
    n_rows       = .N,
    n_enc        = uniqueN(joined_hosp_id),
    n_pat        = uniqueN(patient_id),
    evt_rate_12h = mean(make_y(h_to_event, 12L), na.rm = TRUE),
    evt_rate_24h = mean(make_y(h_to_event, 24L), na.rm = TRUE)
  )]
  print(diag_tbl)
  
  ### standard horizon counts (all observations) -------------------------------
  
  message("  Computing horizon counts...")
  counts_by_point = run_horizon_counts(dt_v, HORIZONS, site_lowercase)
  
  for (HH in sort(unique(counts_by_point$h))) {
    write_artifact(
      df       = counts_by_point[h == HH],
      analysis = if (v == "main") "horizon" else "sensitivity",
      artifact = "counts",
      site     = site_lowercase,
      strata   = "ca",
      horizon  = HH,
      variant  = if (v == "main") NULL else v
    )
  }
  
  ### bootstrapped horizon counts (1 obs per encounter, 100 iterations) --------
  
  message("  Computing bootstrapped horizon counts...")
  counts_boot = run_horizon_counts_bootstrap(dt_v, HORIZONS, site_lowercase, B = 100L)
  
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
  
  ## encounter-level max score analyses ----------------------------------------
  
  dt_max = materialize_variant_max(v, scores_max_enc)
  
  if (!is.null(dt_max)) {
    message("  Writing encounter-level max scores...")
    dt_max_agg = aggregate_maxscores(dt_max, site_lowercase)
    
  if (exists("collapse_small")) {
    dt_max_agg = collapse_small(dt_max_agg, grpvars = c("score_name", "ca_01", "outcome"), val_var = "max_value")
  }
    
    write_artifact(
      df       = dt_max_agg,
      analysis = if (v == "main") "main" else "sensitivity",
      artifact = "maxscores",
      site     = site_lowercase,
      strata   = "ca",
      variant  = if (v == "main") NULL else v
    )
  } else {
    message("  Skipping encounter-level max scores (not applicable for this variant)")
  }
}

message("\n== All variants complete ==")

# subgroup analysis by cancer type ---------------------------------------------

message("\n== Cancer type subgroup analysis ==")

## ensure liquid_01 is in the base data ----------------------------------------

liq = 
  fsubset(cohort, liquid_01 == 1) |>
  pull(joined_hosp_id)
  
ever_positive_complete = 
  fsubset(ever_positive_complete, ca_01 == 1) |>
  ftransform(liquid_01 = if_else(joined_hosp_id %in% liq, 1L, 0L))
  
scores_long_base = 
  fsubset(scores_long_base, ca_01 == 1 & ed_admit_01 == 1) |>
  ftransform(liquid_01 = if_else(joined_hosp_id %in% liq, 1L, 0L))

## max scores by encounter (liquid subgroup) -----------------------------------

dt_max_liquid = 
  fsubset(scores_max_enc, ca_01 == 1 & ed_admit_01 == 1) |>
  ftransform(liquid_01 = if_else(joined_hosp_id %in% liq, 1L, 0L)) |>
  fselect(-patient_id, -joined_hosp_id, -ca_01, -ed_admit_01, -fullcode_01)

dt_max_liquid_agg = 
  as.data.table(dt_max_liquid)[
    , .(n = .N), by = .(score_name, liquid_01, max_value, outcome)
  ][, site := site_lowercase][]

if (exists("collapse_small")) {
  dt_max_liquid_agg = collapse_small(
    dt_max_liquid_agg, 
    grpvars = c("score_name", "liquid_01", "outcome"), 
    val_var = "max_value"
  )
}

write_artifact(
  df       = dt_max_liquid_agg,
  analysis = "main",
  artifact = "maxscores",
  site     = site_lowercase,
  strata   = "liquid"
)

## 24h horizon counts (liquid subgroup) ----------------------------------------

message("  Computing 24h horizon counts by liquid_01...")

dt_liquid_24h = 
  fsubset(scores_long_base, ca_01 == 1 & ed_admit_01 == 1)

dt_liquid_24h[, outcome := make_y(h_to_event, 24L)]

counts_liquid_24h = 
  dt_liquid_24h[, .(n = .N), by = .(score_name, liquid_01, value, outcome)
  ][, `:=`(site = site_lowercase, h = 24L)]

if (exists("collapse_small")) {
  counts_liquid_24h = collapse_small(
    counts_liquid_24h, 
    grpvars = c("score_name", "liquid_01", "outcome")
  )
}

write_artifact(
  df       = counts_liquid_24h,
  analysis = "horizon",
  artifact = "counts",
  site     = site_lowercase,
  strata   = "liquid",
  horizon  = 24L
)

# upset plot data --------------------------------------------------------------

rm(dt_liquid_24h, counts_liquid_24h, dt_max_liquid_agg, dt_max_liquid); gc()

THRESHOLDS$score_name = str_remove_all(THRESHOLDS$score_name, "_total")

upset = 
  fsubset(scores_max_enc, ed_admit_01 == 1) |>
  join(THRESHOLDS, how = "left", multiple = T) |>
  ftransform(positive = if_else(max_value >= threshold, 1L, 0L)) |>
  select(joined_hosp_id, ca_01, outcome, score_name, positive) |>
  pivot_wider(names_from = score_name, values_from = positive)

pooled_upset_counts = {
  x =
    fgroup_by(upset, ca_01, outcome, sirs, qsofa, mews, news) |>
    fsummarise(n = fnobs(joined_hosp_id))
  k = fsum(x$n < 5L)
  message(k, " rows set to 5 (n < 5).")
  x$n = pmax(x$n, 5L)
  roworder(x, -n)
}

write_artifact(
  df       = pooled_upset_counts,
  analysis = "threshold",
  artifact = "upset",
  site     = site_lowercase,
  strata   = "ca",
  horizon  = NULL
)

