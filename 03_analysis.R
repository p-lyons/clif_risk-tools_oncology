
# start main analysis script here

# setup ------------------------------------------------------------------------

## function to write files -----------------------------------------------------

.allowed <- list(
  main        = c("auc","cm"),
  threshold   = c("ever","first","cm","comps","leadtime"),
  sensitivity = c("auc","cm","ever","first","comps","leadtime"),
  subgroups   = c("auc","cm","counts"),
  horizon     = c("auc","cm","counts")
)

# filename: {artifact}{_strata?}{_h{hrs}?}{-variant?}-{site}.csv
.build_filename <- function(artifact, site, strata = NULL, horizon = NULL, variant = NULL) {
  stopifnot(nzchar(artifact), nzchar(site))
  nm <- artifact
  if (!is.null(strata)  && nzchar(strata))  nm <- paste0(nm, "_", strata)
  if (!is.null(horizon) && nzchar(horizon)) nm <- paste0(nm, "_h", horizon)
  if (!is.null(variant) && nzchar(variant)) nm <- paste0(nm, "-", variant)
  paste0(nm, "-", site, ".csv")
}

# write CSV to proj_output/{analysis}/...
write_artifact <- function(df, analysis, artifact, site,
                           strata = NULL, horizon = NULL, variant = NULL,
                           root = "proj_output") {
  stopifnot(analysis %in% names(.allowed),
            artifact %in% .allowed[[analysis]])
  dir <- file.path(root, analysis)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  fn   <- .build_filename(artifact, site, strata, horizon, variant)
  path <- file.path(dir, fn)
  data.table::fwrite(df, path)
  path
}

# max score before outcome -----------------------------------------------------

## get each encounter's maximal pre-outcome score ------------------------------

scores_max = 
  fgroup_by(scores, joined_hosp_id) |>
  fsummarize(
    sirs_max  = fmax(sirs_total),
    qsofa_max = fmax(qsofa_total),
    mews_max  = fmax(mews_total),
    news_max  = fmax(news_total)
  )

scores = join(scores, scores_max, how = "left", multiple = T) 

## data frame for primary analysis ---------------------------------------------

main_out = 
  select(scores, joined_hosp_id, ca_01, o_primary_01, ends_with("max")) |>
  funique() |>
  pivot_longer(
    cols      = ends_with("max"),
    names_to  = "score_name",
    values_to = "value"
  ) |>
  fmutate(score_name = str_remove(score_name, "_max")) |>
  fgroup_by(ca_01, o_primary_01, score_name, value) |>
  fsummarize(n = fnobs(joined_hosp_id))

### ok to collapse a very small number of cells to prevent n < 5 ---------------

#### collapsing function -------------------------------------------------------

collapse_small <- function(dt, grpvars, val_var = "value", n_var = "n", thresh = 5, max_pct_collapse = 0.01) {
  
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
        collapsed_n <<- collapsed_n + ns[i]   # track collapsed counts
        ns[i + 1] <- ns[i + 1] + ns[i]
        keep[i]   <- FALSE
      }
    }
    if (length(ns) > 1L && ns[length(ns)] < thresh) {
      collapsed_n <<- collapsed_n + ns[length(ns)]
      ns[length(ns) - 1L] <- ns[length(ns) - 1L] + ns[length(ns)]
      keep[length(ns)]    <- FALSE
    }
    .(value = vals[keep], n = ns[keep])
  }, by = grpvars]
  
  setnames(dt_out, c("value", "n"), c(val_var, n_var))
  
  pct = collapsed_n / total_n
  
  if (pct > max_pct_collapse) {
    stop(sprintf("❌ Collapsed %.4f%% of observations (>%.0f%% allowed). Review needed.", 
                 100 * pct, 100 * max_pct_collapse))
  } else {
    message(sprintf("✅ Collapsed safely: %.4f%% of observations collapsed (≤%.0f%% allowed).",
                    100 * pct, 100 * max_pct_collapse))
  }
  
  dt_out
}

#### run collapsing function ---------------------------------------------------

collapsed_out = 
  collapse_small(
    dt      = main_out,
    grpvars = c("ca_01","o_primary_01","score_name"),
    val_var = "value",
    n_var   = "n",
    thresh  = 5
  ) |>
  fmutate(site = site_lowercase)

write_artifact(
  df       = collapsed_out,
  analysis = "main",
  artifact = "auc",
  site     = site_lowercase,
  strata   = "ca",
  variant  = "counts"
)

rm(main_out, collapsed_out, scores_max); gc()

## when does the max score occur? ----------------------------------------------

### sirs -----------------------------------------------------------------------
# 
# sirs = 
#   select(scores, joined_hosp_id, time, starts_with("sirs")) |>
#   fsubset(sirs_total == sirs_max) |>
#   roworder(time) |>
#   fgroup_by(joined_hosp_id) |> 
#   ffirst()
#   
# ### qsofa -----------------------------------------------------------------------
# 
# qsofa = 
#   select(scores, joined_hosp_id, time, starts_with("qsofa")) |>
#   fsubset(qsofa_total == qsofa_max) |>
#   roworder(time) |>
#   fgroup_by(joined_hosp_id) |> 
#   ffirst()
# 
# ### mews -----------------------------------------------------------------------
# 
# mews = 
#   select(scores, joined_hosp_id, time, starts_with("mews")) |>
#   fsubset(mews_total == mews_max) |>
#   roworder(time) |>
#   fgroup_by(joined_hosp_id) |> 
#   ffirst()
# 
# ### news -----------------------------------------------------------------------
# 
# news = 
#   select(scores, joined_hosp_id, time, starts_with("news")) |>
#   fsubset(news_total == news_max) |>
#   roworder(time) |>
#   fgroup_by(joined_hosp_id) |> 
#   ffirst()

############

# threshold positivity ---------------------------------------------------------

## establish threshold definitions ---------------------------------------------

score_name = c("sirs_total","qsofa_total","mews_total","news_total")
threshold  = c(2L, 2L, 5L, 5L)
thresh_tbl = bind_cols(score_name = score_name, threshold = threshold)

## first time, if any, crossing each threshold ---------------------------------

scores_long = 
  select(scores, joined_hosp_id, ca_01, o_primary_01, time, ends_with("total")) |>
  pivot_longer(
    cols      = ends_with("total"),
    names_to  = "score_name",
    values_to = "value"
  ) 

scores_thresh = 
  join(scores_long, thresh_tbl, how = "left", multiple = T) |>
  fsubset(value >= threshold) |>
  roworder(time) |>
  fgroup_by(joined_hosp_id, score_name) |>
  ffirst() 

## data for confusion matrix ---------------------------------------------------

scores_long_thresh = 
  join(scores_long, scores_thresh, how = "left", multiple = T) |>
  fmutate(threshold = if_else(is.na(threshold), 0L, 1L)) |>
  roworder(time) |>
  fgroup_by(joined_hosp_id, score_name, threshold) |>
  ffirst()

thresh_counts = 
  fgroup_by(scores_long_thresh, score_name, ca_01, o_primary_01) |>
  fsummarize(
    n_pos = fsum(threshold),
    n_tot = fnobs(joined_hosp_id)
  ) |>
  fmutate(site = site_lowercase)

### save and clean up ----------------------------------------------------------

write_artifact(
  df       = thresh_counts,
  analysis = "threshold",
  artifact = "cm",
  site     = site_lowercase,
  strata   = "ca"
)

rm(threshold, thresh_counts, scores_long_thresh, scores_long); gc()

## time to threshold -----------------------------------------------------------

### add start and stop times to scores_thresh ----------------------------------

thresh_tbl   = select(thresh_tbl, -threshold)
score_thresh = select(scores_thresh, -value, -threshold)

time_df = 
  select(scores, joined_hosp_id, ca_01, in_dttm, end_dttm) |>
  funique() |>
  expand_grid(thresh_tbl) |>
  join(score_thresh, how = "left", multiple = T) |>
  fmutate(
    t_event_h  = as.numeric(difftime(time, in_dttm), "hours"),
    t_censor_h = pmin(as.numeric(difftime(end_dttm,  in_dttm), "hours"), 168),
    event      = if_else(!is.na(t_event_h) & t_event_h <= t_censor_h & t_event_h <= 168, 1L, 0L),
    t_event_h  = if_else(event == 1L, t_event_h, Inf),
    key        = 1L
  ) |>
  select(joined_hosp_id, ca_01, score_name, event, t_event_h, t_censor_h, key)

### time bins for x-axis -------------------------------------------------------

edges = c(
  seq(00L, 008L, 02L),      
  seq(12L, 024L, 04L),   
  seq(36L, 168L, 12L)  
)

bins = tidytable(
  bin_id = seq_along(edges[-1]),
  l      = edges[-length(edges)],
  u      = edges[-1],
  key    = 1L
)

### create shareable cif df ----------------------------------------------------

life_counts = 
  join(bins, time_df, how = "inner", multiple = T) |>
  fmutate(
    at_risk     = as.integer(t_event_h >= l & t_censor_h >= l),
    n_event     = as.integer(event == 1L & t_event_h >= l & t_event_h < u),
    n_discharge = as.integer(t_censor_h >= l & t_censor_h < u & t_censor_h < t_event_h)
  ) |>
  fgroup_by(score_name, ca_01, bin_id, l, u) |>
  fsummarize(
    n_risk      = fsum(at_risk),
    n_event     = fsum(n_event),
    n_discharge = fsum(n_discharge)
  ) |>
  roworder(score_name, ca_01, bin_id) |>
  fmutate(site = site_lowercase)

write_artifact(
  df       = life_counts,
  analysis = "threshold",
  artifact = "ever",
  site     = site_lowercase,
  strata   = "ca",
  variant  = "cif"
)

### summary stats for time to thresholds/outcomes ------------------------------

time_to_thresh = 
  select(scores, joined_hosp_id, in_dttm, outcome_dttm) |>
  join(scores_thresh, how = "inner", multiple = T) |>
  fmutate(time_to_threshold = as.numeric(difftime(time, in_dttm), "hours")) |> 
  fmutate(time_to_outcome   = as.numeric(difftime(outcome_dttm, time), "hours")) |> 
  select(joined_hosp_id, ca_01, o_primary_01, score_name, time_to_threshold, time_to_outcome)

time_to_thresh = 
  fgroup_by(time_to_thresh, score_name) |>
  fsummarize(
    n_thresh       = fnobs(time_to_threshold),
    sum_thresh     = fsum(time_to_threshold),
    sumsq_thresh   = fsum(time_to_threshold^2),
    min_thresh     = fmin(time_to_threshold),
    q25_thresh     = fquantile(time_to_threshold, 0.25),
    median_thresh  = fmedian(time_to_threshold),
    q75_thresh     = fquantile(time_to_threshold, 0.75),
    max_thresh     = fmax(time_to_threshold),
    n_outcome      = fnobs(time_to_outcome),
    sum_outcome    = fsum(time_to_outcome),
    sumsq_outcome  = fsum(time_to_outcome^2),
    min_outcome    = fmin(time_to_outcome),
    q25_outcome    = fquantile(time_to_outcome, 0.25),
    median_outcome = fmedian(time_to_outcome),
    q75_outcome    = fquantile(time_to_outcome, 0.75),
    max_outcome    = fmax(time_to_outcome)
  ) |>
  fmutate(site = site_lowercase)

write_artifact(
  df       = time_to_thresh,
  analysis = "threshold",
  artifact = "leadtime",
  site     = site_lowercase,
  strata   = "ca"
)

rm(time_to_thresh, life_counts, edges, bins, thresh_tbl, scores_thresh); gc()

# time based -------------------------------------------------------------------

## determine whether outcome occurred in subsequent n hours from each time -----

add_outcome_nh <- function(df, n) {
  col = paste0("outcome_", n, "h")
  mutate(df, !!col := if_else(as.numeric(difftime(outcome_dttm, time, units = "hours")) <= n, 1L, 0L, 0L))
}

df =
  select(scores, joined_hosp_id, ca_01, time, ends_with("total"), outcome_dttm) |>
  fsubset(time <= outcome_dttm | is.na(outcome_dttm)) |>
  add_outcome_nh(06) |> 
  add_outcome_nh(12) |> 
  add_outcome_nh(24) 

## put scores long -------------------------------------------------------------


df_long =
  funique(df) |>
  # pivot outcomes into long
  pivot_longer(
    cols          = starts_with("outcome_"),
    names_to      = "h",
    names_pattern = "outcome_(\\d+)h",
    values_to     = "outcome"
  ) |>
  # pivot predictors into long
  pivot_longer(
    cols      = ends_with("total"),
    names_to  = "score",
    values_to = "point"
  ) |>
  # make horizon an integer
  fmutate(h = as.integer(h))

rm(df, scores_long_h, outcomes_long_h); gc()

## analyze time-dependent aurocs -----------------------------------------------

### function for analysis ------------------------------------------------------

analyze_time_series =
  function(
    data, 
    predictors         = score_name, 
    outcomes           = c("outcome_6h", "outcome_12h", "outcome_24h"),
    se_target          = 0.6,
    threshold_override = NULL, # named vector, e.g. c("esm_1" = 0.3, "esm_2" = 0.2)
    ci_method          = "wilson" # method for binomial CIs
  ) {
    
    # Helper function to calculate binomial confidence intervals
    calc_binomial_ci = function(x, n, method = "wilson", conf_level = 0.95) {
      if (n == 0) return(c(NA_real_, NA_real_))
      
      if (method == "wilson") {
        # Wilson score interval
        z = qnorm(1 - (1 - conf_level) / 2)
        p = x / n
        denom = 1 + z^2 / n
        center = (p + z^2 / (2 * n)) / denom
        margin = z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
        return(c(max(0, center - margin), min(1, center + margin)))
      } else {
        # Clopper-Pearson exact method
        if (x == 0) {
          lower = 0
        } else {
          lower = qbeta((1 - conf_level) / 2, x, n - x + 1)
        }
        
        if (x == n) {
          upper = 1
        } else {
          upper = qbeta(1 - (1 - conf_level) / 2, x + 1, n - x)
        }
        return(c(lower, upper))
      }
    }
    
    # All combinations of predictors and outcomes
    analysis_grid = expand.grid(p = predictors, o = outcomes, stringsAsFactors = F)
    
    # Metrics for each combination
    results = map_dfr(1:nrow(analysis_grid), function(i) {
      
      pred    = analysis_grid$p[i]
      out     = analysis_grid$o[i]
      roc_obj = roc(data[[out]], data[[pred]])
      
      if (!is.null(threshold_override) && pred %in% names(threshold_override)) {
        best_thresh   = threshold_override[[pred]]
        achieved_sens = coords(roc_obj, x = best_thresh, input = "threshold", ret = "sensitivity")
      } else {
        
        # Thresholds at target sensitivity
        thresh_data = coords(
          roc_obj,
          x         = "all",
          ret       = c("threshold", "sensitivity"),
          transpose = FALSE
        ) |> 
          as.data.frame()
        
        thresh_data   = thresh_data[order(abs(thresh_data$sensitivity - se_target)), ]
        best_thresh   = thresh_data$threshold[1]
        achieved_sens = thresh_data$sensitivity[1]
      }
      
      pred_class    = if_else(data[[pred]] >= best_thresh, 1L, 0L)
      cm            = table(data[[out]], pred_class)
      
      # Ensure 2x2 matrix
      if(ncol(cm) == 1) {
        missing_class = setdiff(c(0,1), colnames(cm))
        cm            = cbind(cm, matrix(0, nrow = 2, ncol = 1, dimnames = list(NULL, missing_class)))
      }
      
      # Calculate metrics
      TP = cm["1", "1"]
      TN = cm["0", "0"]
      FP = cm["0", "1"]
      FN = cm["1", "0"]
      
      # Calculate point estimates
      sensitivity = TP / (TP + FN)
      specificity = TN / (TN + FP)
      ppv         = ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_)
      npv         = ifelse((TN + FN) > 0, TN / (TN + FN), NA_real_)
      
      # Calculate confidence intervals
      sens_ci = calc_binomial_ci(TP, TP + FN, method = ci_method)
      spec_ci = calc_binomial_ci(TN, TN + FP, method = ci_method)
      ppv_ci  = if((TP + FP) > 0) calc_binomial_ci(TP, TP + FP, method = ci_method) else c(NA_real_, NA_real_)
      npv_ci  = if((TN + FN) > 0) calc_binomial_ci(TN, TN + FN, method = ci_method) else c(NA_real_, NA_real_)
      
      tidytable(
        predictor     = pred,
        outcome       = out,
        auc           = as.numeric(auc(roc_obj)),
        auc_lci       = as.numeric(ci.auc(roc_obj)[1]),
        auc_uci       = as.numeric(ci.auc(roc_obj)[3]),
        sensitivity   = sensitivity,
        sensitivity_lci = sens_ci[1],
        sensitivity_uci = sens_ci[2],
        specificity   = specificity,
        specificity_lci = spec_ci[1],
        specificity_uci = spec_ci[2],
        ppv           = ppv,
        ppv_lci       = ppv_ci[1],
        ppv_uci       = ppv_ci[2],
        npv           = npv,
        npv_lci       = npv_ci[1],
        npv_uci       = npv_ci[2],
        threshold     = best_thresh,
        thresh_src    = ifelse(!is.null(threshold_override) && pred %in% names(threshold_override), "override", "sensitivity_target")
      )
    })
  }

## run function to get results -------------------------------------------------

thresholds = fread(here("output/performance_encounter-level_ohsu.csv"))
t_esm1     = fsubset(thresholds, predictor == "esm1_max") |> pull(threshold_used)
t_esm2     = fsubset(thresholds, predictor == "esm2_max") |> pull(threshold_used)

time_results = 
  analyze_time_series(
    data               = df,
    predictors         = c("esm_1", "esm_2"),
    outcomes           = c("outcome_4h", "outcome_12h", "outcome_1860h"),
    se_target          = 0.6,
    threshold_override = c("esm_1" = t_esm1, "esm_2" = t_esm2),
    ci_method          = "wilson"  # or "exact" for Clopper-Pearson
  ) |>
  select(-thresh_src) |>
  fmutate(outcome = if_else(outcome == "outcome_1860h", "outcome_infinity", outcome))

