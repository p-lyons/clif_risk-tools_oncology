# ==============================================================================
# 04_meta.R
# Meta-analysis of score × cancer interaction terms
# Tests whether score-outcome relationship differs by cancer status
# ==============================================================================

# setup ------------------------------------------------------------------------

library(metafor)
library(ggplot2)

# LOAD AND CHECK DATA ----------------------------------------------------------

message("\n== Loading meta-analysis inputs ==")

coef_data = coef_data_raw
score_sds = score_sds_raw

SCORE_LABS = c(
  "sirs"    = "SIRS",
  "qsofa"   = "QSOFA",
  "mews"    = "MEWS",
  "mews_sf" = "MEWS-SF",
  "news"    = "NEWS"
)

## check convergence -----------------------------------------------------------

if (any(!coef_data$converged)) {
  message("  WARNING: Some models did not converge:")
  print(coef_data[converged == FALSE, .(site, score)])
}

message("  Loaded coefficients from ", uniqueN(coef_data$site), " sites")
message("  Scores: ", paste(unique(coef_data$score), collapse = ", "))

# POOL STANDARD DEVIATIONS -----------------------------------------------------

message("\n== Pooling score SDs across sites ==")

pooled_sds = score_sds[, .(
  pooled_sd = sqrt(sum(sd_score^2 * (n_encounters - 1)) / (sum(n_encounters) - .N)),
  mean_score = sum(mean_score * n_encounters) / sum(n_encounters),
  total_n    = sum(n_encounters)
), by = score]

message("  Pooled SDs calculated:")
print(pooled_sds)

# RANDOM-EFFECTS META-ANALYSIS -------------------------------------------------

message("\n== Running random-effects meta-analysis ==")

meta_results = lapply(unique(coef_data$score), function(s) {
  
  score_data = coef_data[score == s]
  
  message("  Fitting meta-analysis for: ", s, " (", nrow(score_data), " sites)")
  
  # Fit REML with Knapp-Hartung adjustment
  mf = rma(
    yi     = score_data$beta_int,
    sei    = score_data$se_int,
    method = "REML",
    test   = "knha"
  )
  
  pr = predict(mf)
  
  data.table(
    score           = s,
    k_sites         = nrow(score_data),
    pooled_beta     = as.numeric(mf$beta),
    pooled_se       = mf$se,
    pooled_ci_lower = mf$ci.lb,
    pooled_ci_upper = mf$ci.ub,
    pooled_z        = mf$zval,
    pooled_p        = mf$pval,
    tau2            = mf$tau2,
    I2              = mf$I2,
    pred_lower      = pr$pi.lb,
    pred_upper      = pr$pi.ub
  )
}) |> rbindlist()

# MERGE POOLED SDS AND CALCULATE ORs -------------------------------------------

message("\n== Calculating odds ratios ==")

meta_results = merge(meta_results, pooled_sds, by = "score")

meta_results[, `:=`(
  # Per-point ORs
  pooled_or     = exp(pooled_beta),
  or_ci_lower   = exp(pooled_ci_lower),
  or_ci_upper   = exp(pooled_ci_upper),
  or_pred_lower = exp(pred_lower),
  or_pred_upper = exp(pred_upper)
)]

## per-SD translations ---------------------------------------------------------

meta_results[, `:=`(
  beta_per_sd     = pooled_beta * pooled_sd,
  or_per_sd       = exp(pooled_beta * pooled_sd),
  or_per_sd_lower = exp(pooled_ci_lower * pooled_sd),
  or_per_sd_upper = exp(pooled_ci_upper * pooled_sd)
)]

## prespecified equivalence band: OR per 1-SD within [0.95, 1.05] --------------

eq_lo_sd = 0.95
eq_hi_sd = 1.05
meta_results[, equivalent := (or_per_sd_lower > eq_lo_sd & or_per_sd_upper < eq_hi_sd)]

# FORMAT RESULTS TABLE ---------------------------------------------------------

message("\n== Formatting results ==")

meta_table = meta_results[, .(
  score,
  k_sites,
  `Pooled β (95% CI)`     = sprintf("%.3f (%.3f, %.3f)", pooled_beta, pooled_ci_lower, pooled_ci_upper),
  `OR per point (95% CI)` = sprintf("%.3f (%.3f, %.3f)", pooled_or, or_ci_lower, or_ci_upper),
  `OR per 1 SD (95% CI)`  = sprintf("%.3f (%.3f, %.3f)", or_per_sd, or_per_sd_lower, or_per_sd_upper),
  `Prediction Interval`   = sprintf("%.3f to %.3f", or_pred_lower, or_pred_upper),
  `P-value (KH)`          = fifelse(pooled_p < 0.001, "<0.001", sprintf("%.3f", pooled_p)),
  `τ²`                    = round(tau2, 4),
  `I²`                    = sprintf("%.1f%%", I2),
  `Equivalent (per SD)?`  = fifelse(equivalent, "Yes", "No")
)]

print(meta_table)

# INTERPRETATION GUIDE ---------------------------------------------------------

# Negative β: Score has WEAKER association with outcome in cancer patients
#   (each point increase predicts less additional risk in cancer vs non-cancer)
#
# Positive β: Score has STRONGER association with outcome in cancer patients
#
# Equivalence = TRUE: 95% CI falls entirely within equivalence band
#   → Score performs equivalently; AUROC differences reflect baseline risk, not score behavior
#
# Equivalence = FALSE: 95% CI extends outside equivalence band
#   → Score-outcome relationship fundamentally differs by cancer status
#   → Threshold adjustment alone may be insufficient

# SITE-LEVEL COEFFICIENTS FOR FOREST PLOT --------------------------------------

message("\n== Preparing forest plot data ==")

coef_for_plot = merge(coef_data, pooled_sds[, .(score, pooled_sd)], by = "score")

coef_for_plot[, `:=`(
  or_int       = exp(beta_int),
  or_int_lower = exp(beta_int - 1.96 * se_int),
  or_int_upper = exp(beta_int + 1.96 * se_int),
  score_lab    = factor(score, levels = names(SCORE_LABS), labels = SCORE_LABS)
)]

# Add pooled estimates for forest plot
pooled_for_plot = meta_results[, .(
  score,
  site         = "Pooled",
  beta_int     = pooled_beta,
  se_int       = pooled_se,
  or_int       = pooled_or,
  or_int_lower = or_ci_lower,
  or_int_upper = or_ci_upper,
  site_n       = total_n,
  n_events     = NA_integer_,
  converged    = TRUE
)]
pooled_for_plot[, score_lab := factor(score, levels = names(SCORE_LABS), labels = SCORE_LABS)]

forest_data = rbindlist(list(coef_for_plot, pooled_for_plot), fill = TRUE)
forest_data[, is_pooled := site == "Pooled"]

# SLOPE COMPARISON PLOT DATA ---------------------------------------------------

message("\n== Preparing slope comparison data ==")

# Use meta-analyzed coefficients for prediction curves when available
# Falls back to naive pooled glm if site-level coefficients are not exported

has_full_coefs = all(c("beta_intercept", "beta_score", "beta_cancer") %in% names(coef_data))

if (has_full_coefs && nrow(coef_data) > 0) {
  
  message("  Using meta-analyzed coefficients for slope visualization")
  
  score_list   = unique(coef_data$score)
  pred_results = list()
  
  for (s in score_list) {
    
    s_coefs = coef_data[score == s]
    
    # meta-analyze each coefficient across sites
    ma_intercept = tryCatch(
      rma(yi = s_coefs$beta_intercept, sei = s_coefs$se_intercept, method = "REML"),
      error = function(e) NULL
    )
    ma_score = tryCatch(
      rma(yi = s_coefs$beta_score, sei = s_coefs$se_score, method = "REML"),
      error = function(e) NULL
    )
    ma_cancer = tryCatch(
      rma(yi = s_coefs$beta_cancer, sei = s_coefs$se_cancer, method = "REML"),
      error = function(e) NULL
    )
    
    if (is.null(ma_intercept) || is.null(ma_score) || is.null(ma_cancer)) next
    
    b0  = as.numeric(ma_intercept$beta)
    b_s = as.numeric(ma_score$beta)
    b_c = as.numeric(ma_cancer$beta)
    b_i = meta_results[score == s]$pooled_beta
    
    # use maxscores to determine score range
    score_data_s = maxscores_ca_raw[analysis == "main" & score_name == s]
    score_range  = seq(min(score_data_s$max_value), max(score_data_s$max_value), by = 0.5)
    
    pred_no_ca = plogis(b0 + b_s * score_range)
    pred_ca    = plogis(b0 + b_c + (b_s + b_i) * score_range)
    
    pred_results[[s]] = data.table(
      score          = s,
      score_value    = rep(score_range, 2),
      cancer_status  = rep(c("Non-cancer", "Cancer"), each = length(score_range)),
      predicted_prob = c(pred_no_ca, pred_ca)
    )
  }
  
  all_preds = rbindlist(pred_results)
  
  # add interaction OR from meta-analysis for labeling
  all_preds = merge(
    all_preds,
    meta_results[, .(score, or_interaction = pooled_or)],
    by = "score"
  )
  all_preds[, score_label := sprintf("%s\nInteraction OR = %.3f", toupper(score), or_interaction)]
  
} else if (nrow(maxscores_ca_raw[analysis == "main"]) > 0) {
  
  message("  Full coefficients not available; falling back to pooled glm for visualization")
  
  score_data   = maxscores_ca_raw[analysis == "main"]
  score_list   = unique(score_data$score_name)
  pred_results = list()
  
  for (s in score_list) {
    
    s_data = score_data[score_name == s]
    s_data[, score_x_cancer := max_value * ca_01]
    
    mod   = glm(
      cbind(outcome, n - outcome) ~ ca_01 + max_value + score_x_cancer,
      data = s_data, family = binomial()
    )
    coefs = coef(mod)
    
    score_range = seq(min(s_data$max_value), max(s_data$max_value), by = 0.5)
    
    pred_no_ca = plogis(coefs["(Intercept)"] + coefs["max_value"] * score_range)
    pred_ca    = plogis(coefs["(Intercept)"] + coefs["ca_01"] +
                          (coefs["max_value"] + coefs["score_x_cancer"]) * score_range)
    
    pred_results[[s]] = data.table(
      score          = s,
      score_value    = rep(score_range, 2),
      cancer_status  = rep(c("Non-cancer", "Cancer"), each = length(score_range)),
      predicted_prob = c(pred_no_ca, pred_ca)
    )
  }
  
  all_preds = rbindlist(pred_results)
  all_preds = merge(
    all_preds,
    meta_results[, .(score, or_interaction = pooled_or)],
    by = "score"
  )
  all_preds[, score_label := sprintf("%s\nInteraction OR = %.3f", toupper(score), or_interaction)]
  
} else {
  message("  WARNING: No data for slope plots")
  all_preds = data.table()
}

# EXPORTS ----------------------------------------------------------------------

message("\n== Meta-analysis complete ==")

# Objects for figures
meta_results_final = meta_results
meta_table_final   = meta_table
coef_data_final    = coef_for_plot
forest_data_final  = forest_data
slope_preds_final  = all_preds
slope_models_final = all_models
