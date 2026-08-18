# 03c_models.R — longitudinal (population-averaged) model, per site (P7)
#
# Reviewer 2 asked for a longitudinal method that accounts for the repeated,
# within-encounter measurements. For each score this fits a GEE with an
# exchangeable working correlation clustered on the encounter and exports the
# score x cancer interaction (with robust SE) for central meta-analysis. This
# coexists with the glmmTMB interaction export in 03a (meta/coefficients); the
# threshold-equivalence layer supersedes both in the manuscript, but the
# interaction meta-analysis remains available.
#
#   Inputs (proj_tables/):  scores_full.parquet, cohort.parquet
#   Reads (upload_to_box_v2/meta/): 03a's glmmTMB coefficients, for a sign QC
#   Model: outcome_24h ~ value * ca_01 (+ hospital_id where >1 hospital),
#          family = binomial, id = joined_hosp_id, corstr = "exchangeable".
#   Outcome: composite at the 24-hour horizon (D3: ward-to-ICU o_composite).

if (!exists("BOX_DIR")) {
  stop("BOX_DIR not found. Did you run 00_setup first?", call. = FALSE)
}
if (!exists("site_lowercase")) {
  site_lowercase = as.character(read_parquet(here("proj_tables", "site_lowercase.parquet"))$site_lowercase)
}

library(geepack)

# constants --------------------------------------------------------------------

SCORE_COLS   = c("sirs_total", "qsofa_total", "mews_total", "news_total", "mews_sf_total")
GEE_HORIZON  = 24L
MAX_CLUSTERS = 100000L
SAMPLE_SEED  = 2025L

make_y = function(h_to_event, horizon) {
  as.integer(!is.na(h_to_event) & h_to_event >= 0 & h_to_event <= horizon)
}

# build the long, 24h-horizon dataset (composite / main) -----------------------

message("\n== 03c: longitudinal GEE models ==")

scores = as.data.table(read_parquet(here("proj_tables", "scores_full.parquet")))
cohort = as.data.table(read_parquet(here("proj_tables", "cohort.parquet")))

jp = cohort[, .(joined_hosp_id, hospital_id)]

sc = merge(scores[ed_admit_01 == 1L], jp, by = "joined_hosp_id")
sc[, h_to_event := as.numeric(difftime(event_composite_dttm, time, units = "hours"))]
sc[, outcome_24h := make_y(h_to_event, GEE_HORIZON)]

long = melt(
  sc,
  id.vars       = c("joined_hosp_id", "hospital_id", "ca_01", "outcome_24h"),
  measure.vars  = SCORE_COLS,
  variable.name = "score_name", value.name = "value"
)
long[, score_name := as.character(score_name)]

# sample encounters (not observations), keeping all obs of chosen encounters ---

enc_ids           = unique(long$joined_hosp_id)
n_clusters_total  = length(enc_ids)

set.seed(SAMPLE_SEED)
keep_enc = if (n_clusters_total > MAX_CLUSTERS) sample(enc_ids, MAX_CLUSTERS) else enc_ids
sampling_fraction = length(keep_enc) / n_clusters_total

long = long[joined_hosp_id %in% keep_enc]

message(sprintf("  clusters: %s of %s sampled (fraction %.3f)",
                format(length(keep_enc), big.mark = ","),
                format(n_clusters_total, big.mark = ","), sampling_fraction))

# fit one GEE per score --------------------------------------------------------

fit_gee_one = function(d, score, n_clusters_total, sampling_fraction) {

  # geeglm requires clusters to be contiguous
  setorder(d, joined_hosp_id)

  multi_hosp = uniqueN(d$hospital_id) > 1L
  f = if (multi_hosp) {
    as.formula("outcome_24h ~ value * ca_01 + hospital_id")
  } else {
    as.formula("outcome_24h ~ value * ca_01")
  }

  na_row = function(converged) {
    data.table(
      score = score, beta_value = NA_real_, se_value = NA_real_,
      beta_ca = NA_real_, se_ca = NA_real_, beta_int = NA_real_, se_int = NA_real_,
      alpha_working = NA_real_, n_obs = nrow(d), n_clusters = uniqueN(d$joined_hosp_id),
      n_clusters_total = n_clusters_total, sampling_fraction = sampling_fraction,
      converged = converged
    )
  }

  fit = tryCatch(
    geeglm(f, family = binomial(), id = joined_hosp_id, corstr = "exchangeable", data = d),
    error = function(e) { message("    ", score, ": fit error — ", conditionMessage(e)); NULL }
  )
  if (is.null(fit)) return(na_row(FALSE))

  co   = summary(fit)$coefficients
  getc = function(term, col) if (term %in% rownames(co)) co[term, col] else NA_real_
  conv = tryCatch(isTRUE(fit$geese$error == 0), error = function(e) NA)
  alpha = tryCatch(as.numeric(fit$geese$alpha), error = function(e) NA_real_)

  data.table(
    score            = score,
    beta_value       = getc("value", "Estimate"),
    se_value         = getc("value", "Std.err"),
    beta_ca          = getc("ca_01", "Estimate"),
    se_ca            = getc("ca_01", "Std.err"),
    beta_int         = getc("value:ca_01", "Estimate"),
    se_int           = getc("value:ca_01", "Std.err"),
    alpha_working    = if (length(alpha)) alpha[1] else NA_real_,
    n_obs            = nrow(d),
    n_clusters       = uniqueN(d$joined_hosp_id),
    n_clusters_total = n_clusters_total,
    sampling_fraction = sampling_fraction,
    converged        = isTRUE(conv)
  )
}

gee_coefficients = rbindlist(lapply(SCORE_COLS, function(sc_name) {
  message("  fitting: ", sc_name)
  fit_gee_one(long[score_name == sc_name], sc_name, n_clusters_total, sampling_fraction)
}), use.names = TRUE)

gee_coefficients[, site := site_lowercase]

# QC ---------------------------------------------------------------------------

# convergence: report, do not fail (rows already carry converged)
n_nonconv = sum(!gee_coefficients$converged)
if (n_nonconv > 0L) {
  message("  ⚠ ", n_nonconv, " score(s) did not converge (emitted with converged = FALSE).")
} else {
  message("✅ all scores converged.")
}

# sign of interaction vs round-one glmmTMB (meta/coefficients); flag, do not fail
glmm_path = here(BOX_DIR, "meta", paste0("coefficients-", site_lowercase, ".csv"))
if (file.exists(glmm_path)) {
  glmm = fread(glmm_path)[, .(score, glmm_beta_int = beta_int)]
  cmp  = merge(gee_coefficients[, .(score, gee_beta_int = beta_int)], glmm, by = "score")
  cmp[, sign_agree := sign(gee_beta_int) == sign(glmm_beta_int)]
  disagree = cmp[sign_agree == FALSE & !is.na(gee_beta_int) & !is.na(glmm_beta_int)]$score
  if (length(disagree) > 0L) {
    message("  ⚠ GEE/glmmTMB interaction sign disagrees for: ", paste(disagree, collapse = ", "))
  } else {
    message("✅ GEE interaction signs agree with round-one glmmTMB.")
  }
} else {
  message("  (glmmTMB coefficients not found; skipping sign cross-check)")
}

# write ------------------------------------------------------------------------

meta_dir = here(BOX_DIR, "meta")
if (!dir.exists(meta_dir)) dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(gee_coefficients, file.path(meta_dir, paste0("gee_coefficients-", site_lowercase, ".csv")))

message("\n== 03c complete ==")
message("Files written to: ", meta_dir)

################################################################################
