# 03a_artifacts.R — generalized artifact engine (P5)
#
# Generalizes the round-one 03_analysis.R so OUTCOME and STRATUM are parameters
# rather than hard-coded values, and adds an estimability rule to the AUROC
# artifacts. Discharges the outcome-decomposition (R1/R4), calibration (R4.3),
# operating-characteristics (R2.6/R4.2), subgroup (R4.5), and heterogeneity
# (R4.8) requests. The round-one score x cancer interaction export (glmmTMB
# meta/coefficients, meta/score_sds) was removed in this revision: the
# threshold-equivalence analyses supersede it in the manuscript, and no
# reported result depends on it.
#
# ==============================================================================
# PART 1 of 3 — constants, specification objects, and self-contained helpers.
#   Part 2: data prep + parameterized materializers + round-one reproduction.
#   Part 3: run matrix (threshold block, variant loop, STRATA iteration,
#           event tally, bootstrap) + cleanup.
# This part contains definitions only; sourcing it does no work.
# ==============================================================================

# ------------------------------------------------------------------------------
# Meeting decisions baked into this engine:
#   D2 — se_fullcode_only uses initial_code_status %in% c("full","presume_full")
#        in every path. The full-code encounter set `fc` is defined ONCE that way
#        (Part 2, prepare block) and used for both fullcode_01 and the
#        se_fullcode_only subset, so the point-level and encounter-max halves
#        agree (they did not in round one).
#   D3 — the threshold sensitivity/specificity block uses the outcome_times-
#        derived flag (o_composite_01, ward-to-ICU based) — the same flag the
#        AUROC/maxscores path uses — not the cohort-flag composite
#        pmin(wicu_01 + d_noicu_01 + h_noicu_01, 1). (Part 2/3.)
# ------------------------------------------------------------------------------

# constants --------------------------------------------------------------------

THRESHOLDS = tidytable(
  score_name = c("sirs_total", "qsofa_total", "mews_total", "news_total", "mews_sf_total"),
  threshold  = c(2L, 2L, 5L, 5L, 7L)
)

HORIZONS = c(12L, 24L)

VARIANTS = c("main", "se_no_ed_req", "se_fullcode_only", "se_win0_96h", "se_one_enc_per_pt")

# minimum events AND non-events for a site-level AUROC to be reported (change 5)
MIN_ESTIMABLE = 10L

# outcome specification (data, not code). Column names match scores_full.parquet
# as produced by P1 (flags/events) and P2 (end columns).
OUTCOMES = tidytable(
  key       = c("composite", "nohospice", "wardicu", "warddeath", "hospicedc"),
  flag_col  = c("o_composite_01", "o_nohospice_01", "o_wardicu_01",
                "o_warddeath_01", "o_hospicedc_01"),
  event_col = c("event_composite_dttm", "event_nohospice_dttm", "event_wardicu_dttm",
                "event_warddeath_dttm", "event_hospicedc_dttm"),
  end_col   = c("end_composite_dttm", "end_nohospice_dttm", "end_wardicu_dttm",
                "end_warddeath_dttm", "end_hospicedc_dttm")
)

# stratum specification. group_col is the column the stratum splits on; subset is
# the encounter filter applied before splitting. liquid_01/mets_01 are joined
# from cohort.parquet in Part 2 (scores_full carries only ca_01/ed_admit).
STRATA = tidytable(
  key       = c("ca", "liquid", "mets"),
  group_col = c("ca_01", "liquid_01", "mets_01"),
  subset    = c("all", "ca_01 == 1", "ca_01 == 1 & liquid_01 == 0")
)

# allowed artifacts per analysis subdirectory. "events" (Part 3) is the
# per-subgroup event tally for the fine-stratum × all-outcome cells.
.allowed = list(
  main        = c("maxscores", "auroc", "events", "firstscore"),
  threshold   = c("ever", "sesp", "cuminc", "first", "upset"),
  sensitivity = c("maxscores", "counts", "auroc"),
  horizon     = c("counts", "auroc"),
  diagnostics = c("overall", "by_cancer", "max_scores")
)

# filename build + parse -------------------------------------------------------

#' Build an artifact filename.
#' Pattern: {artifact}[-{strata}][-{outcome}][-h{horizon}][-{variant}]-{site}.csv
#' All components except artifact and site are optional. NA is treated as absent.
.build_filename = function(artifact, site, strata = NULL, outcome = NULL,
                           horizon = NULL, variant = NULL) {

  stopifnot(nzchar(artifact), nzchar(site))
  ok = function(x) !is.null(x) && !is.na(x) && nzchar(as.character(x))

  parts = artifact
  if (ok(strata))  parts = paste0(parts, "-", strata)
  if (ok(outcome)) parts = paste0(parts, "-", outcome)
  if (!is.null(horizon) && !is.na(horizon)) parts = paste0(parts, "-h", horizon)
  if (ok(variant)) parts = paste0(parts, "-", variant)

  paste0(parts, "-", site, ".csv")
}

#' Invert .build_filename(): recover the generating parameters from a filename.
#' The token vocabularies (strata, outcome, horizon) are disjoint and the token
#' order is fixed, so any middle token that is not a stratum, outcome, or horizon
#' is the variant.
.parse_filename = function(fn) {

  strata_vocab  = c(STRATA$key, "components")
  outcome_vocab = OUTCOMES$key

  stem = sub("\\.csv$", "", fn)
  toks = strsplit(stem, "-", fixed = TRUE)[[1]]
  n    = length(toks)

  out = list(artifact = toks[1], strata = NA_character_, outcome = NA_character_,
             horizon = NA_integer_, variant = NA_character_, site = toks[n])

  mid = if (n > 2) toks[2:(n - 1)] else character(0)
  for (tk in mid) {
    if (tk %in% strata_vocab) {
      out$strata = tk
    } else if (tk %in% outcome_vocab) {
      out$outcome = tk
    } else if (grepl("^h[0-9]+$", tk)) {
      out$horizon = as.integer(sub("^h", "", tk))
    } else {
      out$variant = tk
    }
  }
  out
}

#' QC helper: a filename must round-trip through the parser to the same string.
#' Used by the run loop to satisfy the "every emitted filename parses back to its
#' generating parameters" assertion.
.assert_filename_roundtrip = function(fn) {
  p = .parse_filename(fn)
  rebuilt = .build_filename(p$artifact, p$site, strata = p$strata, outcome = p$outcome,
                            horizon = p$horizon, variant = p$variant)
  if (!identical(rebuilt, fn)) {
    stop(sprintf("Filename does not round-trip: '%s' -> '%s'", fn, rebuilt), call. = FALSE)
  }
  invisible(TRUE)
}

# make_y (unchanged from round one) --------------------------------------------

make_y = function(h_to_event, horizon) {
  as.integer(!is.na(h_to_event) & h_to_event >= 0 & h_to_event <= horizon)
}

# site-level AUROC with DeLong CI + estimability (change 5) ---------------------
# Generalized over the grouping column so the same routine serves every stratum
# (group_col = "ca_01" reproduces round one). A cell with fewer than
# MIN_ESTIMABLE events OR non-events returns NA for the estimate and interval but
# retains n_obs / n_events and is flagged estimable = FALSE, so the coordinating
# center sees the denominators and excludes only the estimate from pooling.
compute_site_auroc = function(dt, score_col = "value", site_lowercase, group_col = "ca_01") {

  library(pROC)

  dt = as.data.table(dt)
  combos = unique(dt[, .(score_name, .grp = get(group_col))])

  mk_row = function(sc, gval, n_obs, n_events, auroc, se, lo, hi, estimable) {
    r = data.table(score_name = sc)
    r[, (group_col) := gval]
    r[, `:=`(n_obs = n_obs, n_events = n_events, auroc = auroc, auroc_se = se,
             ci_lower = lo, ci_upper = hi, site = site_lowercase, estimable = estimable)]
    r
  }

  results = lapply(seq_len(nrow(combos)), function(i) {

    sc = combos$score_name[i]
    g  = combos$.grp[i]

    sub = dt[score_name == sc & get(group_col) == g]

    n_obs    = nrow(sub)
    n_events = sum(sub$outcome == 1L, na.rm = TRUE)
    n_nonev  = n_obs - n_events

    estimable = (n_events >= MIN_ESTIMABLE) && (n_nonev >= MIN_ESTIMABLE)

    if (!estimable) {
      return(mk_row(sc, g, n_obs, n_events, NA_real_, NA_real_, NA_real_, NA_real_, FALSE))
    }

    roc_obj = tryCatch(
      roc(sub$outcome, sub[[score_col]], levels = c(0, 1), direction = "<", quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(roc_obj)) {
      return(mk_row(sc, g, n_obs, n_events, NA_real_, NA_real_, NA_real_, NA_real_, FALSE))
    }

    auc_val = as.numeric(auc(roc_obj))
    ci_obj  = tryCatch(ci.auc(roc_obj, method = "delong"), error = function(e) NULL)

    if (!is.null(ci_obj)) {
      ci_lo    = as.numeric(ci_obj[1])
      ci_hi    = as.numeric(ci_obj[3])
      auroc_se = (ci_hi - ci_lo) / (2 * 1.96)
    } else {
      ci_lo = ci_hi = auroc_se = NA_real_
    }

    mk_row(sc, g, n_obs, n_events, auc_val, auroc_se, ci_lo, ci_hi, TRUE)
  })

  rbindlist(results, use.names = TRUE)
}

# aggregate encounter-max scores to counts, generalized over grouping column ----
# The `outcome` column here is the 0/1 event indicator; the outcome KEY
# (composite, wardicu, …) is carried in the filename, not the table.
aggregate_maxscores = function(dt, site_lowercase, group_col = "ca_01") {

  dt = as.data.table(dt)
  need = c("score_name", group_col, "max_value", "outcome")
  miss = setdiff(need, names(dt))
  if (length(miss)) stop("aggregate_maxscores(): missing cols: ", paste(miss, collapse = ", "))

  dt[, .(n = .N), by = c("score_name", group_col, "max_value", "outcome")
  ][, site := site_lowercase][]
}

# ==============================================================================
# END PART 1
# ==============================================================================

# ==============================================================================
# PART 2 of 3 — setup, prepare, parameterized materializers, reproduction check.
# The materializers below are written in data.table so outcome/stratum column
# names can be resolved dynamically; because scores_full is zero-filled (no NA in
# the score totals), the per-encounter maxima are identical to the round-one
# collapse fmax() results.
# ==============================================================================

if (!exists("BOX_DIR")) {
  stop("BOX_DIR not found. Did you run 00_setup first?", call. = FALSE)
}

# read inputs ------------------------------------------------------------------

scores = read_parquet(here("proj_tables", "scores_full.parquet"))
cohort = read_parquet(here("proj_tables", "cohort.parquet"))

if (!exists("site_lowercase")) {
  site_lowercase = as.character(read_parquet(here("proj_tables", "site_lowercase.parquet"))$site_lowercase)
}

# clean up only the artifacts THIS script regenerates --------------------------
# Leaves outputs written earlier in the run by 02 (news_o2_resolution), 02b
# (*-cf{H}), and 02c (monitoring*, missing_vlab-ca) in place. In sensitivity/,
# where this script and 02b share the maxscores/counts base names, deletion is
# restricted to the four analysis-variant tokens so the carry-forward files stay.

cleanup_artifact_dirs = function(root = BOX_DIR) {

  variant_tokens = c("se_no_ed_req", "se_fullcode_only", "se_win0_96h", "se_one_enc_per_pt")

  for (subdir in names(.allowed)) {
    dir_path = here(root, subdir)
    if (!dir.exists(dir_path)) next

    bases = .allowed[[subdir]]

    pattern = if (subdir == "sensitivity") {
      paste0("^(", paste(bases, collapse = "|"), ")-.*-(",
             paste(variant_tokens, collapse = "|"), ")-.*\\.csv$")
    } else {
      paste0("^(", paste(bases, collapse = "|"), ")-.*\\.csv$")
    }

    old_files = list.files(dir_path, pattern = pattern, full.names = TRUE)
    if (length(old_files) > 0) {
      message("  Removing ", length(old_files), " regenerated file(s) from ", subdir, "/")
      file.remove(old_files)
    }
  }
  invisible(NULL)
}

message("\n== 03a: preparing datasets ==")
cleanup_artifact_dirs()

# prepare ----------------------------------------------------------------------
# D2: the full-code encounter set is full + presume_full, defined once here and
# used for both fullcode_01 and the se_fullcode_only subset. The cohort join adds
# the subgroup columns (liquid_01, mets_01, rank_enc) that scores_full lacks.

jp = select(cohort, patient_id, joined_hosp_id, hospital_id, liquid_01, mets_01, rank_enc)

fc = fsubset(cohort, tolower(initial_code_status) %in% c("full", "presume_full"))$joined_hosp_id

scores = join(scores, jp, how = "inner", multiple = FALSE) |>
  ftransform(h_from_admit = as.numeric(difftime(time, in_dttm, units = "hours")))

setDT(scores)

# sampling helpers (unchanged from round one) ----------------------------------

sample_one_encounter_per_patient = function(enc_df) {
  dt = as.data.table(enc_df)
  set.seed(2025L)
  dt[, .SD[sample(.N, 1L)], by = patient_id]$joined_hosp_id
}

# horizon counts, generalized over grouping column -----------------------------
# group_col = "ca_01" reproduces round one.

run_horizon_counts = function(dt, horizons, site_lowercase, group_col = "ca_01") {
  counts_list = lapply(horizons, function(HH) {
    dt_copy = copy(dt)
    dt_copy[, outcome := make_y(h_to_event, HH)]
    dt_copy[, .(n = .N), by = c("score_name", group_col, "value", "outcome")
    ][, `:=`(site = site_lowercase, h = HH)]
  })
  rbindlist(counts_list, use.names = TRUE)[]
}

# variant row-subset shared by both materializers ------------------------------

.apply_variant = function(dt, variant, fc, flag, evc) {
  switch(
    variant,
    main             = dt[ed_admit_01 == 1L],
    se_no_ed_req     = copy(dt),
    se_fullcode_only = dt[ed_admit_01 == 1L & joined_hosp_id %in% fc],
    se_win0_96h      = {
      d = dt[ed_admit_01 == 1L & h_from_admit >= 0 & h_from_admit <= 96]
      # censor events occurring beyond 96h from ward admission
      h_ev = as.numeric(difftime(d[[evc]], d$in_dttm, units = "hours"))
      d[[evc]]  = fifelse(!is.na(h_ev) & h_ev <= 96, d[[evc]], as.POSIXct(NA))
      d[[flag]] = fifelse(!is.na(h_ev) & h_ev <= 96, d[[flag]], 0L)
      d
    },
    se_one_enc_per_pt = {
      enc  = unique(dt[ed_admit_01 == 1L, .(joined_hosp_id, patient_id)])
      keep = sample_one_encounter_per_patient(enc)
      dt[joined_hosp_id %in% keep]
    },
    stop("Unknown variant: ", variant)
  )
}

.apply_stratum = function(dt, subset_expr) {
  if (subset_expr == "all") return(dt)
  keep = eval(parse(text = subset_expr), envir = dt)
  dt[keep]
}

# materialize encounter-level max scores for (variant, outcome, stratum) -------

materialize_variant_max = function(variant, scores, fc, outcome, stratum) {

  flag = outcome$flag_col; endc = outcome$end_col; evc = outcome$event_col
  gcol = stratum$group_col; subset_expr = stratum$subset

  dt = .apply_variant(as.data.table(scores), variant, fc, flag, evc)
  dt = .apply_stratum(dt, subset_expr)
  dt = dt[time < get(endc)]                       # per-outcome truncation
  if (nrow(dt) == 0L) return(NULL)

  dt_max = dt[, .(
    patient_id    = patient_id[1L],
    hospital_id   = hospital_id[1L],
    grp           = get(gcol)[1L],
    ed_admit_01   = ed_admit_01[1L],
    fullcode_01   = as.integer(joined_hosp_id[1L] %in% fc),
    sirs_max      = max(sirs_total,    na.rm = TRUE),
    qsofa_max     = max(qsofa_total,   na.rm = TRUE),
    mews_max      = max(mews_total,    na.rm = TRUE),
    news_max      = max(news_total,    na.rm = TRUE),
    mews_sf_max   = max(mews_sf_total, na.rm = TRUE),
    enc_news_any3 = max(news_any3,     na.rm = TRUE),
    outcome       = get(flag)[1L]
  ), by = joined_hosp_id]

  dt_max = melt(dt_max,
    measure.vars  = c("sirs_max", "qsofa_max", "mews_max", "news_max", "mews_sf_max"),
    variable.name = "score_name", value.name = "max_value")
  dt_max[, score_name := sub("_max$", "", as.character(score_name))]
  dt_max[is.infinite(max_value),    max_value := NA_integer_]
  dt_max[is.infinite(enc_news_any3), enc_news_any3 := 0L]
  setnames(dt_max, "grp", gcol)
  dt_max[]
}

# materialize point-level (long) rows for (variant, outcome, stratum) ----------

materialize_variant_long = function(variant, scores, fc, outcome, stratum) {

  flag = outcome$flag_col; endc = outcome$end_col; evc = outcome$event_col
  gcol = stratum$group_col; subset_expr = stratum$subset

  dt = .apply_variant(as.data.table(scores), variant, fc, flag, evc)
  dt = .apply_stratum(dt, subset_expr)
  dt = dt[time < get(endc)]
  if (nrow(dt) == 0L) return(NULL)

  dt[, h_to_event   := as.numeric(difftime(get(evc), time,    units = "hours"))]
  dt[, h_from_admit := as.numeric(difftime(time,     in_dttm, units = "hours"))]

  long = melt(dt,
    id.vars       = c("joined_hosp_id", "patient_id", gcol, "h_to_event", "h_from_admit"),
    measure.vars  = c("sirs_total", "qsofa_total", "mews_total", "news_total", "mews_sf_total"),
    variable.name = "score_name", value.name = "value")
  long[, score_name := as.character(score_name)]
  long[]
}

# reproduction self-check ------------------------------------------------------
# Confirms the generalized materializer, run for outcome = composite / stratum =
# ca / variant = main, reproduces a hard-coded round-one-style computation
# (which keys on o_primary_01 and truncates at end_dttm). This proves the
# parameterization is faithful before the engine runs on the new outcomes and
# subgroups. NOTE: the definitive comparison against the April submission is the
# separate byte-identical check against commit 37897cf.

.verify_composite_reproduction = function(scores, fc, site_lowercase) {

  oc = OUTCOMES[key == "composite"]
  st = STRATA[key == "ca"]

  gen = aggregate_maxscores(
    materialize_variant_max("main", scores, fc, oc, st),
    site_lowercase, group_col = "ca_01"
  )

  ref_scores = as.data.table(scores)[ed_admit_01 == 1L][time < end_dttm]
  ref_max = ref_scores[, .(
    ca_01       = ca_01[1L],
    sirs_max    = max(sirs_total,    na.rm = TRUE),
    qsofa_max   = max(qsofa_total,   na.rm = TRUE),
    mews_max    = max(mews_total,    na.rm = TRUE),
    news_max    = max(news_total,    na.rm = TRUE),
    mews_sf_max = max(mews_sf_total, na.rm = TRUE),
    outcome     = o_primary_01[1L]
  ), by = joined_hosp_id]
  ref_max = melt(ref_max,
    measure.vars  = c("sirs_max", "qsofa_max", "mews_max", "news_max", "mews_sf_max"),
    variable.name = "score_name", value.name = "max_value")
  ref_max[, score_name := sub("_max$", "", as.character(score_name))]
  ref_max[is.infinite(max_value), max_value := NA_integer_]
  ref = ref_max[, .(n = .N), by = .(score_name, ca_01, max_value, outcome)][, site := site_lowercase][]

  if (!fsetequal(gen, ref)) {
    stop("Generalized composite/ca/main maxscores do not reproduce the round-one reference.",
         call. = FALSE)
  }
  message("✅ Reproduction self-check passed: generalized composite/ca/main == round-one reference.")
  invisible(TRUE)
}

.verify_composite_reproduction(scores, fc, site_lowercase)

# ==============================================================================
# END PART 2  (Part 3: run matrix + threshold block + event tally + bootstrap)
# ==============================================================================

# ==============================================================================
# PART 3 of 3 (3a) — run matrix: maxscores / auroc / counts / events.
# The threshold block (ever/sesp/cuminc/first/upset) and the bootstrap are
# ports of the round-one blocks with the per-outcome flag swap, and follow in
# Part 3b.
# ==============================================================================

# artifact writer: builds the name, checks it round-trips, writes ---------------

write_artifact = function(df, analysis, artifact, site,
                          strata = NULL, outcome = NULL, horizon = NULL, variant = NULL,
                          root = BOX_DIR) {

  stopifnot(analysis %in% names(.allowed), artifact %in% .allowed[[analysis]])

  dir = here(root, analysis)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  fn = .build_filename(artifact, site, strata = strata, outcome = outcome,
                       horizon = horizon, variant = variant)
  .assert_filename_roundtrip(fn)   # QC: every emitted filename parses back to its params

  fwrite(df, file.path(dir, fn))
  invisible(fn)
}

# emitters ---------------------------------------------------------------------

emit_maxscores_auroc = function(scores, fc, variant, outcome, stratum, analysis, site_lowercase) {

  dt_max = materialize_variant_max(variant, scores, fc, outcome, stratum)
  if (is.null(dt_max)) return(invisible())

  gcol = stratum$group_col
  vtok = if (variant == "main") NULL else variant

  agg = aggregate_maxscores(dt_max, site_lowercase, group_col = gcol)
  write_artifact(agg, analysis, "maxscores", site_lowercase,
                 strata = stratum$key, outcome = outcome$key, variant = vtok)

  au = compute_site_auroc(dt_max, score_col = "max_value", site_lowercase, group_col = gcol)
  au[, metric := "encounter_max"]
  write_artifact(au, analysis, "auroc", site_lowercase,
                 strata = stratum$key, outcome = outcome$key, variant = vtok)

  invisible()
}

emit_counts_auroc_horizon = function(scores, fc, variant, outcome, stratum, horizons, analysis, site_lowercase) {

  long = materialize_variant_long(variant, scores, fc, outcome, stratum)
  if (is.null(long)) return(invisible())

  gcol = stratum$group_col
  vtok = if (variant == "main") NULL else variant

  counts = run_horizon_counts(long, horizons, site_lowercase, group_col = gcol)

  for (HH in horizons) {
    write_artifact(counts[h == HH], analysis, "counts", site_lowercase,
                   strata = stratum$key, outcome = outcome$key, horizon = HH, variant = vtok)

    dt_h = copy(long)
    dt_h[, outcome := make_y(h_to_event, HH)]
    au_h = compute_site_auroc(dt_h, score_col = "value", site_lowercase, group_col = gcol)
    au_h[, horizon := HH]
    write_artifact(au_h, analysis, "auroc", site_lowercase,
                   strata = stratum$key, outcome = outcome$key, horizon = HH, variant = vtok)
  }
  invisible()
}

# event tally (row 4): per-subgroup n_enc / n_events, no discrimination estimate

make_event_tally = function(scores, fc, outcome, stratum, site_lowercase) {

  flag = outcome$flag_col; endc = outcome$end_col
  gcol = stratum$group_col; subset_expr = stratum$subset

  dt = as.data.table(scores)[ed_admit_01 == 1L]
  dt = .apply_stratum(dt, subset_expr)
  dt = dt[time < get(endc)]
  if (nrow(dt) == 0L) return(NULL)

  enc = dt[, .(grp = get(gcol)[1L], ev = get(flag)[1L]), by = joined_hosp_id]
  out = enc[, .(n_enc = .N, n_events = sum(ev == 1L, na.rm = TRUE)), by = grp]
  setnames(out, "grp", gcol)
  out[, site := site_lowercase][]
}

# first-ward-observation score distribution (Reviewer 1, Major 5) --------------
# Each score's value at the FIRST ward observation per ED-admit encounter, as a
# count distribution by cancer status. "First ward observation" is the earliest
# scores_full timestamp for the encounter (scores_full is already ward-restricted
# and outcome-truncated, and carries one row per encounter-time). Fully poolable:
# the coordinating center sums n across sites per (score_name, ca_01, value) and
# recovers medians, IQRs, and means. Descriptive only -- no within-stratum
# discrimination is computed here, which would be circular.

make_firstscore_distribution = function(scores, site_lowercase) {

  score_cols = c("sirs_total", "qsofa_total", "mews_total", "news_total", "mews_sf_total")

  dt = as.data.table(scores)[ed_admit_01 == 1L]
  if (nrow(dt) == 0L) return(NULL)

  setorderv(dt, c("joined_hosp_id", "time"))
  first = dt[, .SD[1L], by = joined_hosp_id, .SDcols = c("ca_01", score_cols)]

  long = melt(first,
              id.vars       = c("joined_hosp_id", "ca_01"),
              measure.vars  = score_cols,
              variable.name = "score_name", value.name = "value")
  long[, score_name := as.character(score_name)]

  dist = long[!is.na(value), .(n = .N), by = .(score_name, ca_01, value)]
  setorder(dist, score_name, ca_01, value)
  dist[, site := site_lowercase][]
}

# run matrix -------------------------------------------------------------------

run_artifact_matrix = function(scores, fc, site_lowercase) {

  st_ca = STRATA[key == "ca"]

  # Row 1 — ca, all five outcomes, main variant
  for (k in OUTCOMES$key) {
    oc = OUTCOMES[key == k]
    emit_maxscores_auroc(scores, fc, "main", oc, st_ca, "main", site_lowercase)
    emit_counts_auroc_horizon(scores, fc, "main", oc, st_ca, HORIZONS, "horizon", site_lowercase)
  }

  # Row 2 — ca, composite only, four sensitivity variants
  oc_comp = OUTCOMES[key == "composite"]
  for (v in setdiff(VARIANTS, "main")) {
    emit_maxscores_auroc(scores, fc, v, oc_comp, st_ca, "sensitivity", site_lowercase)
    emit_counts_auroc_horizon(scores, fc, v, oc_comp, st_ca, HORIZONS, "sensitivity", site_lowercase)
  }

  # Row 3 — subgroups, composite + nohospice, main variant (h24 counts only)
  for (sk in c("liquid", "mets")) {
    st = STRATA[key == sk]
    for (k in c("composite", "nohospice")) {
      oc = OUTCOMES[key == k]
      emit_maxscores_auroc(scores, fc, "main", oc, st, "main", site_lowercase)
      emit_counts_auroc_horizon(scores, fc, "main", oc, st, 24L, "horizon", site_lowercase)
    }
  }

  # Row 4 — subgroups, all five outcomes, event tally only (no AUROC)
  for (sk in c("liquid", "mets")) {
    st = STRATA[key == sk]
    for (k in OUTCOMES$key) {
      oc = OUTCOMES[key == k]
      ev = make_event_tally(scores, fc, oc, st, site_lowercase)
      if (!is.null(ev)) {
        write_artifact(ev, "main", "events", site_lowercase, strata = st$key, outcome = oc$key)
      }
    }
  }
  invisible()
}

# cross-outcome QC (also flags exact-time ties, D7) ----------------------------

cross_outcome_qc = function(scores) {

  dt0 = as.data.table(scores)[ed_admit_01 == 1L]

  encN = sapply(OUTCOMES$key, function(k) {
    oc = OUTCOMES[key == k]
    uniqueN(dt0[time < get(oc$end_col)]$joined_hosp_id)
  })
  evN = sapply(OUTCOMES$key, function(k) {
    oc = OUTCOMES[key == k]
    d  = dt0[time < get(oc$end_col)]
    e  = d[, .(v = get(oc$flag_col)[1L]), by = joined_hosp_id]
    sum(e$v == 1L, na.rm = TRUE)
  })

  if (length(unique(encN)) != 1L) {
    stop("Encounter counts disagree across outcomes.", call. = FALSE)
  }
  if (evN[["composite"]] != evN[["wardicu"]] + evN[["warddeath"]] + evN[["hospicedc"]]) {
    stop("composite events != wardicu + warddeath + hospicedc (check exact-time ties, D7).", call. = FALSE)
  }
  if (evN[["nohospice"]] != evN[["wardicu"]] + evN[["warddeath"]]) {
    stop("nohospice events != wardicu + warddeath (check exact-time ties, D7).", call. = FALSE)
  }
  message("✅ Cross-outcome QC passed (encounter counts agree; composite = components; nohospice = wardicu + warddeath).")
  invisible()
}

# run --------------------------------------------------------------------------

message("\n== 03a: running artifact matrix ==")
run_artifact_matrix(scores, fc, site_lowercase)
cross_outcome_qc(scores)

# first-ward score distribution by cancer status (Reviewer 1, Major 5) ---------
firstscore = make_firstscore_distribution(scores, site_lowercase)
if (!is.null(firstscore)) {
  write_artifact(firstscore, "main", "firstscore", site_lowercase, strata = "ca")
  message(sprintf("  first-ward score distribution: %d rows (score x cancer x value)", nrow(firstscore)))
}

# Part 3b (threshold block + bootstrap) continues below.

# ==============================================================================
# END PART 3a
# ==============================================================================

# ==============================================================================
# PART 3b — threshold block (ca and liquid, per outcome) + bootstrap.
# Ported from round-one 03_analysis.R. The score-crossing computations
# (ever_positive, firsts) are outcome-independent and computed once; the
# outcome-dependent aggregations loop over the five outcomes. D3: the outcome
# flag is the outcome_times-derived o_{k}_01 (ward-to-ICU composite), taken from
# outcome_times so the sesp/ever denominators include the full cohort ED-admit
# set (matching round one), not just encounters that carry score rows.
# ==============================================================================

# per-encounter outcome flags for the full cohort ED-admit set (all five)
outcome_tt = as.data.table(read_parquet(here("proj_tables", "outcome_times.parquet")))

all_encs_base =
  fsubset(cohort, ed_admit_01 == 1) |>
  select(joined_hosp_id, ca_01, liquid_01) |>
  join(
    select(outcome_tt, joined_hosp_id,
           o_composite_01, o_nohospice_01, o_wardicu_01, o_warddeath_01, o_hospicedc_01),
    how = "left", multiple = FALSE
  )
setDT(all_encs_base)
# encounters with no recorded outcome are non-events for every outcome
for (fcol in OUTCOMES$flag_col) all_encs_base[is.na(get(fcol)), (fcol) := 0L]

# --- ever positive: first threshold crossing per encounter x score (once) -----

message("\n== 03a: threshold block ==")

ever_positive =
  fsubset(scores, ed_admit_01 == 1) |>
  pivot_longer(cols = ends_with("total"), names_to = "score_name", values_to = "value") |>
  join(THRESHOLDS, how = "inner", multiple = FALSE) |>
  ftransform(positive = value >= threshold) |>
  fsubset(positive == TRUE) |>
  roworder(time) |>
  fgroup_by(joined_hosp_id, score_name) |>
  fsummarize(
    time_to_positive_h = ffirst(h_from_admit),
    first_positive_val = ffirst(value)
  )

news_any3_pos =
  fsubset(scores, ed_admit_01 == 1 & news_any3 == 1L) |>
  roworder(time) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(
    time_to_positive_h = ffirst(h_from_admit),
    first_positive_val = ffirst(news_total)
  ) |>
  ftransform(score_name = "news_total")

ever_positive =
  rowbind(ever_positive, news_any3_pos) |>
  roworder(time_to_positive_h) |>
  fgroup_by(joined_hosp_id, score_name) |>
  fsummarize(
    time_to_positive_h = ffirst(time_to_positive_h),
    first_positive_val = ffirst(first_positive_val)
  )

rm(news_any3_pos)

# --- first-positive components per encounter x score (once) -------------------

extract_first_positive = function(scores_dt, score_prefix, threshold_val, score_total_col) {
  scores_dt |>
    fsubset(get(score_total_col) >= threshold_val & ed_admit_01 == 1) |>
    select(joined_hosp_id, ca_01, h_from_admit, matches(paste0("^", score_prefix, "_"))) |>
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

first_sirs  = extract_first_positive(scores, "sirs",  2L, "sirs_total")
first_qsofa = extract_first_positive(scores, "qsofa", 2L, "qsofa_total")
first_mews  = extract_first_positive(scores, "mews",  5L, "mews_total")

first_news =
  scores |>
  fsubset((news_total >= 5L | news_any3 == 1L) & ed_admit_01 == 1) |>
  select(joined_hosp_id, ca_01, h_from_admit, matches("^news_")) |>
  select(-matches("_total$"), -news_any3) |>
  roworder(h_from_admit) |>
  fgroup_by(joined_hosp_id) |>
  ffirst() |>
  pivot_longer(cols = matches("^news_"), names_to = "component", values_to = "value", names_prefix = "news_") |>
  ftransform(score = "news") |>
  fsubset(value > 0)

first_mewssf =
  fsubset(scores, mews_sf_total >= 7L & ed_admit_01 == 1) |>
  select(joined_hosp_id, ca_01, h_from_admit, matches("^mews_"), sf) |>
  select(-matches("_total$")) |>
  roworder(h_from_admit) |>
  fgroup_by(joined_hosp_id) |>
  ffirst() |>
  rename(mews_sf = sf) |>
  pivot_longer(cols = matches("^mews_"), names_to = "component", values_to = "value", names_prefix = "mews_") |>
  ftransform(score = "mews_sf") |>
  fsubset(value > 0)

firsts = rowbind(first_sirs, first_qsofa, first_mews, first_news, first_mewssf)
rm(first_sirs, first_qsofa, first_mews, first_news, first_mewssf)

# encounter follow-up windows for cumulative incidence (outcome-independent)
enc_windows =
  fsubset(scores, ed_admit_01 == 1) |>
  fgroup_by(joined_hosp_id) |>
  fsummarize(
    ca_01              = ffirst(ca_01),
    h_admit_to_dc      = ffirst(as.numeric(difftime(out_dttm, in_dttm, units = "hours"))),
    max_obs_time       = fmax(h_from_admit, na.rm = TRUE)
  )
setDT(enc_windows)

time_breaks = seq(0, 168, by = 8)

first_pos_times =
  select(firsts, joined_hosp_id, score, h_from_admit) |>
  funique() |>
  rename(h_first_positive = h_from_admit)

# --- per-outcome threshold artifacts ------------------------------------------

emit_threshold_for_outcome = function(outcome_key, flag_col) {

  # attach this outcome's flag as o_out to the full ED-admit denominator
  ae = all_encs_base[, .(joined_hosp_id, ca_01, o_out = get(flag_col))]

  # ever positive (complete grid over all encounters x scores)
  epc =
    tidyr::expand_grid(joined_hosp_id = ae$joined_hosp_id, score_name = THRESHOLDS$score_name) |>
    as_tidytable() |>
    join(ae,            how = "left", multiple = FALSE) |>
    join(ever_positive, how = "left", multiple = FALSE) |>
    ftransform(ever_positive = as.integer(!is.na(time_to_positive_h))) |>
    select(joined_hosp_id, score_name, ca_01, ever_positive, o_out)
  epc = as.data.table(epc)

  ever_agg = epc[, .(n = .N), by = .(score_name, ca_01, ever_positive, o_out)][, site := site_lowercase][]
  write_artifact(ever_agg, "threshold", "ever", site_lowercase, strata = "ca", outcome = outcome_key)

  sesp = epc[, .(
    n_total      = .N,
    n_outcome    = sum(o_out == 1L, na.rm = TRUE),
    n_no_outcome = sum(o_out == 0L, na.rm = TRUE),
    n_pos        = sum(ever_positive == 1L, na.rm = TRUE),
    n_neg        = sum(ever_positive == 0L, na.rm = TRUE),
    tp           = sum(ever_positive == 1L & o_out == 1L, na.rm = TRUE),
    fp           = sum(ever_positive == 1L & o_out == 0L, na.rm = TRUE),
    tn           = sum(ever_positive == 0L & o_out == 0L, na.rm = TRUE),
    fn           = sum(ever_positive == 0L & o_out == 1L, na.rm = TRUE)
  ), by = .(score_name, ca_01)
  ][, `:=`(
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    ppv         = tp / (tp + fp),
    npv         = tn / (tn + fn),
    site        = site_lowercase
  )][]
  write_artifact(sesp, "threshold", "sesp", site_lowercase, strata = "ca", outcome = outcome_key)

  # firsts joined to this outcome
  firsts_o = join(firsts, ae[, .(joined_hosp_id, o_out)], how = "left", multiple = FALSE)

  time_to_positive =
    select(firsts_o, joined_hosp_id, score, ca_01, o_out, h_from_admit) |>
    funique() |>
    ftransform(pos_day0 = fifelse(h_from_admit <= 24, 1L, 0L)) |>
    ftransform(pos_day1 = fifelse(h_from_admit > 24 & h_from_admit <= 48, 1L, 0L)) |>
    fgroup_by(score, ca_01, o_out) |>
    fsummarize(
      n                  = fnobs(joined_hosp_id),
      h_from_admit_sum   = fsum(h_from_admit),
      h_from_admit_sumsq = fsum(h_from_admit^2),
      n_pos_day0         = fsum(pos_day0),
      n_pos_day1         = fsum(pos_day1)
    ) |>
    ftransform(site = site_lowercase)
  write_artifact(time_to_positive, "threshold", "first", site_lowercase, strata = "ca", outcome = outcome_key)

  # upset components
  n_positive =
    select(firsts_o, joined_hosp_id, ca_01, o_out, score) |>
    funique() |>
    fgroup_by(ca_01, o_out, score) |>
    fsummarize(n_encs = GRPN())
  upset_components =
    select(firsts_o, joined_hosp_id, ca_01, o_out, score, component) |>
    funique() |>
    join(n_positive, how = "left", multiple = TRUE) |>
    fgroup_by(ca_01, o_out, score, component) |>
    fsummarize(n = GRPN(), n_encs = fmax(n_encs)) |>
    ftransform(site = site_lowercase)
  write_artifact(upset_components, "threshold", "upset", site_lowercase, strata = "components", outcome = outcome_key)

  # cumulative incidence of positivity
  ew = merge(enc_windows, ae[, .(joined_hosp_id, o_out)], by = "joined_hosp_id", all.x = TRUE)
  ew[, h_followup_end := pmin(h_admit_to_dc, max_obs_time, na.rm = TRUE)]

  cuminc_data =
    tidyr::expand_grid(
      joined_hosp_id = ew$joined_hosp_id,
      score          = unique(firsts$score),
      time_bin_start = time_breaks
    ) |>
    as_tidytable() |>
    join(ew,              how = "left", multiple = FALSE) |>
    join(first_pos_times, how = "left", multiple = FALSE) |>
    ftransform(
      at_risk         = time_bin_start < h_followup_end,
      became_positive = !is.na(h_first_positive) & h_first_positive <= time_bin_start
    )
  cuminc_summary =
    fsubset(cuminc_data, at_risk == TRUE) |>
    fgroup_by(score, ca_01, o_out, time_bin_start) |>
    fsummarize(
      n_at_risk    = GRPN(),
      n_became_pos = fsum(became_positive, na.rm = TRUE),
      cum_inc      = fsum(became_positive, na.rm = TRUE) / GRPN()
    ) |>
    ftransform(site = site_lowercase)
  write_artifact(cuminc_summary, "threshold", "cuminc", site_lowercase, strata = "ca", outcome = outcome_key)

  # upset co-positivity (encounter-level), per outcome
  THRESHOLDS_short = copy(as.data.table(THRESHOLDS))
  THRESHOLDS_short[, score_name := sub("_total$", "", score_name)]
  dtm = materialize_variant_max("main", scores, fc, OUTCOMES[key == outcome_key], STRATA[key == "ca"])
  upset_enc =
    join(dtm, THRESHOLDS_short, how = "left", multiple = TRUE) |>
    ftransform(positive = fifelse(
      score_name == "news",
      fifelse(max_value >= threshold | enc_news_any3 == 1L, 1L, 0L),
      fifelse(max_value >= threshold, 1L, 0L)
    )) |>
    select(joined_hosp_id, ca_01, outcome, score_name, positive) |>
    pivot_wider(names_from = score_name, values_from = positive)
  upset_enc = as.data.table(upset_enc)
  upset_counts = upset_enc[, .(n = .N), by = .(ca_01, outcome, sirs, qsofa, mews, news, mews_sf)][, site := site_lowercase]
  setorder(upset_counts, -n)
  write_artifact(upset_counts, "threshold", "upset", site_lowercase, strata = "ca", outcome = outcome_key)

  invisible()
}

for (k in OUTCOMES$key) {
  emit_threshold_for_outcome(k, OUTCOMES[key == k]$flag_col)
}

# ==============================================================================
# THRESHOLD BLOCK, LIQUID STRATUM (cancer only: hematologic vs solid)
# Reviewer 4 asked for threshold-specific operating characteristics for
# hematologic versus solid malignancies. The denominator is the cancer ED-admit
# encounter set, grouped on liquid_01 (1 = hematologic, 0 = solid), and the
# positivity definition is identical to the ca block above (standard thresholds,
# including the NEWS single-parameter rule via ever_positive). Composite and
# nohospice outcomes, matching the subgroup rows of the run matrix.
# ==============================================================================

message("\n== 03a: threshold block, liquid stratum ==")

for (outcome_key_l in c("composite", "nohospice")) {

  flag_col_l = OUTCOMES[key == outcome_key_l]$flag_col

  ae_l = all_encs_base[ca_01 == 1L, .(joined_hosp_id, liquid_01, o_out = get(flag_col_l))]

  epc_l =
    tidyr::expand_grid(joined_hosp_id = ae_l$joined_hosp_id, score_name = THRESHOLDS$score_name) |>
    as_tidytable() |>
    join(ae_l,          how = "left", multiple = FALSE) |>
    join(ever_positive, how = "left", multiple = FALSE) |>
    ftransform(ever_positive = as.integer(!is.na(time_to_positive_h))) |>
    select(joined_hosp_id, score_name, liquid_01, ever_positive, o_out)
  epc_l = as.data.table(epc_l)

  ever_agg_l = epc_l[, .(n = .N), by = .(score_name, liquid_01, ever_positive, o_out)][, site := site_lowercase][]
  write_artifact(ever_agg_l, "threshold", "ever", site_lowercase, strata = "liquid", outcome = outcome_key_l)

  sesp_l = epc_l[, .(
    n_total      = .N,
    n_outcome    = sum(o_out == 1L, na.rm = TRUE),
    n_no_outcome = sum(o_out == 0L, na.rm = TRUE),
    n_pos        = sum(ever_positive == 1L, na.rm = TRUE),
    n_neg        = sum(ever_positive == 0L, na.rm = TRUE),
    tp           = sum(ever_positive == 1L & o_out == 1L, na.rm = TRUE),
    fp           = sum(ever_positive == 1L & o_out == 0L, na.rm = TRUE),
    tn           = sum(ever_positive == 0L & o_out == 0L, na.rm = TRUE),
    fn           = sum(ever_positive == 0L & o_out == 1L, na.rm = TRUE)
  ), by = .(score_name, liquid_01)
  ][, `:=`(
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    ppv         = tp / (tp + fp),
    npv         = tn / (tn + fn),
    site        = site_lowercase
  )][]
  write_artifact(sesp_l, "threshold", "sesp", site_lowercase, strata = "liquid", outcome = outcome_key_l)

  rm(ae_l, epc_l, ever_agg_l, sesp_l)
}


# --- bootstrap horizon counts (composite / ca / main) -------------------------
# Cluster bootstrap over encounters. Encounters are drawn with replacement and
# EVERY score observation belonging to a drawn encounter enters the resample, so
# the bootstrap distribution belongs to the same point-level quantity the
# accompanying horizon counts report. (The prior version drew a single random
# observation per encounter-score, once, outside the iteration loop; each
# iteration therefore resampled the same fixed one-row-per-encounter table and
# summed to the encounter count rather than the observation count.)
#
# Expanding all observations of every drawn encounter is done arithmetically
# rather than physically. Each encounter contributes a fixed count to each
# (score_name, ca_01, value, outcome, h) cell, so multiplying those cell counts
# by the encounter's resampled frequency and summing is identical to materializing
# the expanded table, and keeps the per-iteration join to the collapsed cell set.

message("\n== 03a: bootstrap horizon counts ==")

run_horizon_counts_bootstrap = function(dt, horizons, site_lowercase, B = 400L) {

  set.seed(2025L)

  base_dt = as.data.table(dt)[, .(joined_hosp_id, score_name, ca_01, value, h_to_event)]
  for (HH in horizons) base_dt[, (paste0("y_", HH)) := make_y(h_to_event, HH)]

  ## collapse to one row per encounter x cell, per horizon
  cell_list = vector("list", length(horizons))

  for (i in seq_along(horizons)) {
    HH   = horizons[i]
    ycol = paste0("y_", HH)
    cc = base_dt[
      , .(n_rows = .N),
      by = c("joined_hosp_id", "score_name", "ca_01", "value", ycol)
    ]
    setnames(cc, ycol, "outcome")
    cc[, h := HH]
    cell_list[[i]] = cc
  }

  cells = rbindlist(cell_list, use.names = TRUE)
  setkey(cells, joined_hosp_id)

  enc_ids = unique(base_dt$joined_hosp_id)
  n_enc   = length(enc_ids)

  ## QC: the unweighted cell set must reproduce the point-level observation count
  n_obs_cells = cells[h == horizons[1L], sum(n_rows)]
  n_obs_base  = nrow(base_dt)

  if (n_obs_cells != n_obs_base) {
    stop(
      sprintf("Bootstrap cell collapse lost rows: %d cells vs %d observations.",
              n_obs_cells, n_obs_base),
      call. = FALSE
    )
  }

  boot_list = vector("list", B)

  for (b in seq_len(B)) {
    boot_encs = data.table(joined_hosp_id = sample(enc_ids, n_enc, replace = TRUE))
    enc_freq  = boot_encs[, .(weight = .N), by = joined_hosp_id]
    setkey(enc_freq, joined_hosp_id)

    samp = cells[enc_freq, nomatch = NULL]

    boot_list[[b]] = samp[
      , .(n = sum(weight * n_rows)),
      by = .(score_name, ca_01, value, outcome, h)
    ][, `:=`(site = site_lowercase, iter = b)]

    if (b %% 50 == 0) message("    bootstrap ", b, "/", B)
  }

  rbindlist(boot_list, use.names = TRUE)[]
}


long_main = materialize_variant_long("main", scores, fc, OUTCOMES[key == "composite"], STRATA[key == "ca"])
counts_boot = run_horizon_counts_bootstrap(long_main, HORIZONS, site_lowercase, B = 400L)
for (HH in HORIZONS) {
  write_artifact(counts_boot[h == HH], "horizon", "counts", site_lowercase,
                 strata = "ca", outcome = "composite", horizon = HH, variant = "boot")
}

message("\n== 03a complete ==")
message("Files written to: ", here(BOX_DIR))

# ==============================================================================
# END PART 3b  — 03a complete
# ==============================================================================
