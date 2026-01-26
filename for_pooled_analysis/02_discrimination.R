# ==============================================================================
# 02_discrimination.R
# AUROC calculations using two complementary approaches:
#   1. Site-level meta-analysis (proper DeLong CIs at each site, then pooled)
#   2. Pooled weighted AUC (from aggregated counts, Hanley-McNeil CIs)
# ==============================================================================

# setup ------------------------------------------------------------------------

library(pROC)
library(WeightedROC)
library(metafor)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Calculate weighted AUC from aggregated count data (no expansion needed)
#' Uses WeightedROC package for memory-efficient computation
#' @param df data.table with columns: value (score), outcome (0/1), n (count)
#' @return list with auc, n_total, n_events
calculate_weighted_auc = function(df) {
  
  if (nrow(df) == 0) return(list(auc = NA_real_, n_total = 0L, n_events = 0L))
  
  n_total  = sum(df$n, na.rm = TRUE)
  n_events = sum(df$n[df$outcome == 1], na.rm = TRUE)
  
  if (n_events == 0 || n_events == n_total) {
    return(list(auc = NA_real_, n_total = n_total, n_events = n_events))
  }
  
  # WeightedROC expects: predictor scores, true labels (0/1), weights
  # We have aggregated data, so we pass value as score, outcome as label, n as weight
  tp_fp = tryCatch(
    WeightedROC(df$value, df$outcome, df$n),
    error = function(e) NULL
  )
  
  if (is.null(tp_fp)) {
    return(list(auc = NA_real_, n_total = n_total, n_events = n_events))
  }
  
  auc_val = WeightedAUC(tp_fp)
  
  list(auc = auc_val, n_total = n_total, n_events = n_events)
}

#' Calculate Hanley-McNeil standard error for AUC
#' Analytical approximation that doesn't require individual observations
#' @param auc the AUC value
#' @param n_pos number of positive cases (events)
#' @param n_neg number of negative cases (non-events)
#' @return standard error of AUC
hanley_mcneil_se = function(auc, n_pos, n_neg) {
  
  if (is.na(auc) || n_pos <= 0 || n_neg <= 0) return(NA_real_)
  
  # Hanley & McNeil (1982) approximation
  q1 = auc / (2 - auc)
  q2 = (2 * auc^2) / (1 + auc)
  
  var_auc = (auc * (1 - auc) + (n_pos - 1) * (q1 - auc^2) + (n_neg - 1) * (q2 - auc^2)) / 
            (n_pos * n_neg)
  
  sqrt(var_auc)
}

#' Meta-analyze site-level AUROCs using random-effects model
#' @param auroc_dt data.table with site, auroc, auroc_se columns
#' @return data.table with pooled estimate, CI, heterogeneity stats
meta_analyze_aurocs = function(auroc_dt) {
  
  # Filter to valid estimates
  dt = auroc_dt[!is.na(auroc) & !is.na(auroc_se) & auroc_se > 0]
  
  if (nrow(dt) < 2) {
    # Can't meta-analyze with <2 studies
    if (nrow(dt) == 1) {
      return(data.table(
        k_sites     = 1L,
        auroc       = dt$auroc[1],
        auroc_se    = dt$auroc_se[1],
        ci_lower    = dt$auroc[1] - 1.96 * dt$auroc_se[1],
        ci_upper    = dt$auroc[1] + 1.96 * dt$auroc_se[1],
        tau2        = NA_real_,
        I2          = NA_real_,
        method      = "single_site"
      ))
    }
    return(data.table(
      k_sites = 0L, auroc = NA_real_, auroc_se = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      tau2 = NA_real_, I2 = NA_real_, method = "insufficient_data"
    ))
  }
  
  # Random-effects meta-analysis
  ma = tryCatch(
    rma(yi = dt$auroc, sei = dt$auroc_se, method = "REML"),
    error = function(e) NULL
  )
  
  if (is.null(ma)) {
    return(data.table(
      k_sites = nrow(dt), auroc = NA_real_, auroc_se = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      tau2 = NA_real_, I2 = NA_real_, method = "meta_failed"
    ))
  }
  
  data.table(
    k_sites   = nrow(dt),
    auroc     = as.numeric(ma$beta),
    auroc_se  = ma$se,
    ci_lower  = ma$ci.lb,
    ci_upper  = ma$ci.ub,
    tau2      = ma$tau2,
    I2        = ma$I2,
    method    = "random_effects_meta"
  )
}

# ==============================================================================
# APPROACH 1: SITE-LEVEL META-ANALYSIS
# Uses DeLong CIs computed at each site, then pools via random-effects
# ==============================================================================

message("\n== Approach 1: Site-level meta-analysis ==")

## Encounter-level AUROCs ------------------------------------------------------

message("\n  -- Encounter-level AUROCs --")

if (exists("auroc_enc_raw") && nrow(auroc_enc_raw) > 0) {
  
  # Meta-analyze by score × cancer × analysis
  enc_meta_list = list()
  
  combos = unique(auroc_enc_raw[, .(score_name, ca_01, analysis)])
  
  for (i in seq_len(nrow(combos))) {
    
    sub = auroc_enc_raw[
      score_name == combos$score_name[i] & 
      ca_01 == combos$ca_01[i] & 
      analysis == combos$analysis[i]
    ]
    
    ma_result = meta_analyze_aurocs(sub)
    ma_result[, `:=`(
      score_name = combos$score_name[i],
      ca_01      = combos$ca_01[i],
      analysis   = combos$analysis[i],
      metric     = "Encounter max",
      n_total    = sum(sub$n_obs, na.rm = TRUE),
      n_events   = sum(sub$n_events, na.rm = TRUE)
    )]
    
    enc_meta_list[[i]] = ma_result
  }
  
  auroc_enc_meta = rbindlist(enc_meta_list, use.names = TRUE, fill = TRUE)
  
  message("  Encounter-level meta-analysis: ", nrow(auroc_enc_meta), " estimates")
  
} else {
  message("  WARNING: No site-level encounter AUROCs available")
  auroc_enc_meta = data.table()
}

## 24-hour horizon AUROCs ------------------------------------------------------

message("\n  -- 24-hour horizon AUROCs --")

if (exists("auroc_h24_raw") && nrow(auroc_h24_raw) > 0) {
  
  h24_meta_list = list()
  
  combos = unique(auroc_h24_raw[, .(score_name, ca_01, analysis)])
  
  for (i in seq_len(nrow(combos))) {
    
    sub = auroc_h24_raw[
      score_name == combos$score_name[i] & 
      ca_01 == combos$ca_01[i] & 
      analysis == combos$analysis[i]
    ]
    
    ma_result = meta_analyze_aurocs(sub)
    ma_result[, `:=`(
      score_name = combos$score_name[i],
      ca_01      = combos$ca_01[i],
      analysis   = combos$analysis[i],
      metric     = "24-hour horizon",
      n_total    = sum(sub$n_obs, na.rm = TRUE),
      n_events   = sum(sub$n_events, na.rm = TRUE)
    )]
    
    h24_meta_list[[i]] = ma_result
  }
  
  auroc_h24_meta = rbindlist(h24_meta_list, use.names = TRUE, fill = TRUE)
  
  message("  24h horizon meta-analysis: ", nrow(auroc_h24_meta), " estimates")
  
} else {
  message("  WARNING: No site-level 24h AUROCs available")
  auroc_h24_meta = data.table()
}

## Combine site-level meta-analysis results ------------------------------------

auroc_meta_all = rbindlist(list(auroc_enc_meta, auroc_h24_meta), use.names = TRUE, fill = TRUE)

if (nrow(auroc_meta_all) > 0) {
  auroc_meta_all[, approach := "site_meta"]
}

# ==============================================================================
# APPROACH 2: POOLED WEIGHTED AUC
# Uses aggregated counts with WeightedROC, Hanley-McNeil CIs
# ==============================================================================

message("\n== Approach 2: Pooled weighted AUC ==")

## Encounter-level (from maxscores) --------------------------------------------

message("\n  -- Encounter-level AUROCs (weighted) --")

analysis_variants = c("main", "fullcode_only", "no_ed_req", "win0_96h", "one_enc_per_pt")

enc_weighted_list = list()

for (variant in analysis_variants) {
  
  df_variant = maxscores_ca_raw[analysis == variant]
  
  if (nrow(df_variant) == 0) {
    message("    ", variant, ": 0 rows, skipping")
    next
  }
  
  combos = unique(df_variant[, .(score_name, ca_01)])
  
  for (j in seq_len(nrow(combos))) {
    
    sub = df_variant[score_name == combos$score_name[j] & ca_01 == combos$ca_01[j]]
    
    # Rename for calculate_weighted_auc
    sub_renamed = sub[, .(value = max_value, outcome, n)]
    
    result = calculate_weighted_auc(sub_renamed)
    
    # Hanley-McNeil SE
    n_pos = result$n_events
    n_neg = result$n_total - result$n_events
    se_hm = hanley_mcneil_se(result$auc, n_pos, n_neg)
    
    enc_weighted_list[[length(enc_weighted_list) + 1]] = data.table(
      score_name = combos$score_name[j],
      ca_01      = combos$ca_01[j],
      analysis   = variant,
      metric     = "Encounter max",
      auroc      = result$auc,
      auroc_se   = se_hm,
      ci_lower   = result$auc - 1.96 * se_hm,
      ci_upper   = result$auc + 1.96 * se_hm,
      n_total    = result$n_total,
      n_events   = result$n_events,
      method     = "weighted_hanley_mcneil"
    )
  }
  
  message("    ", variant, ": ", nrow(df_variant), " rows")
}

auroc_enc_weighted = rbindlist(enc_weighted_list, use.names = TRUE, fill = TRUE)
message("  Encounter-level weighted: ", nrow(auroc_enc_weighted), " estimates")

## 24-hour horizon (from counts) -----------------------------------------------

message("\n  -- 24-hour horizon AUROCs (weighted) --")

h24_weighted_list = list()

for (variant in analysis_variants) {
  
  df_variant = counts_h24_raw[analysis == variant]
  
  if (nrow(df_variant) == 0) {
    message("    ", variant, ": 0 rows, skipping")
    next
  }
  
  combos = unique(df_variant[, .(score_name, ca_01)])
  
  for (j in seq_len(nrow(combos))) {
    
    sub = df_variant[score_name == combos$score_name[j] & ca_01 == combos$ca_01[j]]
    
    result = calculate_weighted_auc(sub)
    
    n_pos = result$n_events
    n_neg = result$n_total - result$n_events
    se_hm = hanley_mcneil_se(result$auc, n_pos, n_neg)
    
    h24_weighted_list[[length(h24_weighted_list) + 1]] = data.table(
      score_name = combos$score_name[j],
      ca_01      = combos$ca_01[j],
      analysis   = variant,
      metric     = "24-hour horizon",
      auroc      = result$auc,
      auroc_se   = se_hm,
      ci_lower   = result$auc - 1.96 * se_hm,
      ci_upper   = result$auc + 1.96 * se_hm,
      n_total    = result$n_total,
      n_events   = result$n_events,
      method     = "weighted_hanley_mcneil"
    )
  }
  
  message("    ", variant, ": ", nrow(df_variant), " rows")
}

auroc_h24_weighted = rbindlist(h24_weighted_list, use.names = TRUE, fill = TRUE)
message("  24h horizon weighted: ", nrow(auroc_h24_weighted), " estimates")

## Combine weighted results ----------------------------------------------------

auroc_weighted_all = rbindlist(list(auroc_enc_weighted, auroc_h24_weighted), use.names = TRUE, fill = TRUE)

if (nrow(auroc_weighted_all) > 0) {
  auroc_weighted_all[, approach := "pooled_weighted"]
}

# ==============================================================================
# COMPARE APPROACHES
# ==============================================================================

message("\n== Comparing approaches ==")

if (nrow(auroc_meta_all) > 0 && nrow(auroc_weighted_all) > 0) {
  
  # Merge for comparison
  comparison = merge(
    auroc_meta_all[, .(score_name, ca_01, analysis, metric, 
                       auroc_meta = auroc, se_meta = auroc_se, ci_lo_meta = ci_lower, ci_hi_meta = ci_upper)],
    auroc_weighted_all[, .(score_name, ca_01, analysis, metric,
                           auroc_weighted = auroc, se_weighted = auroc_se, ci_lo_weighted = ci_lower, ci_hi_weighted = ci_upper)],
    by = c("score_name", "ca_01", "analysis", "metric"),
    all = TRUE
  )
  
  comparison[, `:=`(
    diff_auroc = auroc_meta - auroc_weighted,
    diff_se    = se_meta - se_weighted
  )]
  
  message("\n  Comparison of main analysis (meta vs weighted):")
  print(comparison[analysis == "main", .(metric, score_name, ca_01, 
                                          auroc_meta = round(auroc_meta, 3),
                                          auroc_weighted = round(auroc_weighted, 3),
                                          diff = round(diff_auroc, 4))])
  
  auroc_comparison = comparison
  
} else {
  message("  Cannot compare - one or both approaches have no data")
  auroc_comparison = data.table()
}

# ==============================================================================
# PRIMARY RESULTS (use site-level meta if available, else weighted)
# ==============================================================================

message("\n== Assembling primary results ==")

# Prefer site-level meta-analysis (proper CIs), fall back to weighted
if (nrow(auroc_meta_all) > 0) {
  auroc_primary = copy(auroc_meta_all)
  message("  Using site-level meta-analysis as primary")
} else {
  auroc_primary = copy(auroc_weighted_all)
  message("  Using pooled weighted AUC as primary (no site-level data)")
}

# Add labels
auroc_primary[, `:=`(
  ca_lab       = fifelse(ca_01 == 1L, "Cancer", "Non-cancer"),
  analysis_lab = factor(analysis, levels = names(ANALYSIS_LABS), labels = ANALYSIS_LABS),
  score_lab    = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS),
  metric_lab   = factor(metric, levels = c("Encounter max", "24-hour horizon"))
)]

# ==============================================================================
# AUROC DIFFERENCES (Cancer vs Non-cancer)
# ==============================================================================

message("\n== Computing AUROC differences ==")

auroc_diff = auroc_primary[, {
  
  auc_ca    = auroc[ca_01 == 1]
  auc_nonca = auroc[ca_01 == 0]
  se_ca     = auroc_se[ca_01 == 1]
  se_nonca  = auroc_se[ca_01 == 0]
  
  if (length(auc_ca) == 1 && length(auc_nonca) == 1 && 
      !is.na(auc_ca) && !is.na(auc_nonca)) {
    
    diff_val = auc_ca - auc_nonca
    
    # SE of difference (assuming independence between groups)
    se_diff = sqrt(se_ca^2 + se_nonca^2)
    
    # Z-test for difference
    z_val = diff_val / se_diff
    p_val = 2 * pnorm(-abs(z_val))
    
    list(
      auroc_ca     = auc_ca,
      auroc_nonca  = auc_nonca,
      diff_auc     = diff_val,
      se_diff      = se_diff,
      p_value      = p_val
    )
  } else {
    list(
      auroc_ca    = NA_real_,
      auroc_nonca = NA_real_,
      diff_auc    = NA_real_,
      se_diff     = NA_real_,
      p_value     = NA_real_
    )
  }
}, by = .(score_name, analysis, metric)]

auroc_diff[, `:=`(
  sig = fcase(
    is.na(p_value), "",
    p_value < 0.001, "***",
    p_value < 0.01,  "**",
    p_value < 0.05,  "*",
    default = ""
  ),
  score_lab    = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS),
  analysis_lab = factor(analysis, levels = names(ANALYSIS_LABS), labels = ANALYSIS_LABS),
  metric_lab   = factor(metric, levels = c("Encounter max", "24-hour horizon"))
)]

# ==============================================================================
# SUMMARY TABLES
# ==============================================================================

message("\n== Creating summary tables ==")

## Main analysis table ---------------------------------------------------------

auroc_main_table = auroc_primary[analysis == "main", .(
  metric_lab, score_lab, ca_lab,
  auroc_fmt = sprintf("%.3f (%.3f-%.3f)", auroc, ci_lower, ci_upper),
  n         = format_n(n_total),
  n_events  = format_n(n_events)
)]

auroc_main_wide = dcast(
  auroc_main_table,
  metric_lab + score_lab ~ ca_lab,
  value.var = "auroc_fmt"
)

## Difference table ------------------------------------------------------------

auroc_diff_table = auroc_diff[analysis == "main", .(
  metric_lab, score_lab,
  diff_fmt = sprintf("%.3f", diff_auc),
  p_fmt    = fifelse(!is.na(p_value) & p_value < 0.001, "<0.001", 
                     fifelse(!is.na(p_value), sprintf("%.3f", p_value), "—")),
  sig
)]

# ==============================================================================
# EXPORTS
# ==============================================================================

message("\n== Discrimination analysis complete ==")

if (nrow(auroc_meta_all) > 0) {
  message("  Site-level meta-analysis estimates: ", nrow(auroc_meta_all))
}
if (nrow(auroc_weighted_all) > 0) {
  message("  Pooled weighted estimates: ", nrow(auroc_weighted_all))
}
message("  Primary results: ", nrow(auroc_primary))

# Export objects for downstream scripts
auroc_results_final    = auroc_primary
auroc_meta_final       = auroc_meta_all
auroc_weighted_final   = auroc_weighted_all
auroc_comparison_final = auroc_comparison
auroc_diff_final       = auroc_diff
auroc_main_table_final = auroc_main_wide
auroc_diff_table_final = auroc_diff_table

# Also export site-level data for forest plots
if (exists("auroc_enc_raw")) auroc_site_enc = auroc_enc_raw
if (exists("auroc_h24_raw")) auroc_site_h24 = auroc_h24_raw

gc()
