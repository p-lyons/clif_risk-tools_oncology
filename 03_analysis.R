
# start main analysis script here

# setup ------------------------------------------------------------------------

## function to write files -----------------------------------------------------

.allowed <- list(
  main       = c("auc","cm"),
  threshold  = c("ever","first","cm","comps","leadtime"),
  sensitivity= c("auc","cm","ever","first","comps","leadtime"),
  subgroups  = c("auc","cm","counts"),
  horizon    = c("auc","cm","counts")
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
  fn  <- .build_filename(artifact, site, strata, horizon, variant)
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

fwrite(collapsed_out, here("proj_output", "main_output.csv"))

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

fwrite(thresh_counts, here("proj_output", paste0("main_thresh_", site_lowercase, ".csv")))

rm(thresh_counts, scores_long_thresh, scores_long); gc()

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

fwrite(life_counts, here("proj_output", paste0("ci_thresh_", site_lowercase, ".csv")))

### summary stats for time to thresholds/outcomes ------------------------------

time_to_thresh = 
  select(scores, joined_hosp_id, in_dttm, end_dttm) |>
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

fwrite(time_to_thresh, here("proj_output", paste0("time_to_thresh_", site_lowercase, ".csv")))

rm(time_to_thresh, life_counts, edges, bins, thresh_tbl, scores_thresh); gc()

# time based -------------------------------------------------------------------



