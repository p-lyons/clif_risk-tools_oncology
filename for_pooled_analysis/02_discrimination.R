# 02_discrimination.R

# setup ------------------------------------------------------------------------

library(pROC)

# helper functions -------------------------------------------------------------

#' Calculate AUROC with DeLong CIs from aggregated count data
#' For encounter-level data where independence assumption holds
calculate_auroc_delong = function(df_input, analysis_name, score_col = "max_value") {
  
  if (nrow(df_input) == 0) {
    message("  Warning: No data for analysis: ", analysis_name)
    return(data.table())
  }
  
  combos = unique(df_input[, .(score_name, ca_01)])
  auroc_results = list()
  
  for (i in seq_len(nrow(combos))) {
    score  = combos$score_name[i]
    cancer = combos$ca_01[i]
    
    score_data = df_input[score_name == score & ca_01 == cancer]
    
    if (nrow(score_data) == 0) next
    
    # Check outcome variation
    n_pos = score_data[outcome == 1, sum(n, na.rm = TRUE)]
    n_neg = score_data[outcome == 0, sum(n, na.rm = TRUE)]
    total_n = n_pos + n_neg
    
    if (n_pos == 0 | n_neg == 0) {
      message("  Warning: No outcome variation for ", analysis_name, " ", score, " ca_01=", cancer)
      next
    }
    
    # Expand to individual observations
    df_expanded = score_data[rep(1:.N, n)]
    
    # ROC on expanded data
    roc_obj = roc(
      response  = df_expanded$outcome,
      predictor = df_expanded[[score_col]],
      quiet     = TRUE,
      levels    = c(0, 1),
      direction = "<"
    )
    
    ci_delong = ci.auc(roc_obj, method = "delong")
    
    auroc_results[[length(auroc_results) + 1]] = data.table(
      analysis   = analysis_name,
      score_name = score,
      ca_01      = cancer,
      n          = total_n,
      n_events   = n_pos,
      auroc      = as.numeric(auc(roc_obj)),
      ci_lower   = ci_delong[1],
      ci_upper   = ci_delong[3]
    )
    
    rm(df_expanded)
  }
  
  gc()
  rbindlist(auroc_results, use.names = TRUE)
}

delong_comparison = function(df_input, analysis_name, score_col = "max_value") {
  
  if (nrow(df_input) == 0) {
    return(data.table())
  }
  
  scores_list = unique(df_input$score_name)
  results = list()
  
  for (score in scores_list) {
    
    score_data = df_input[score_name == score]
    
    # Check both groups exist with outcome variation
    data_c0 = score_data[ca_01 == 0]
    data_c1 = score_data[ca_01 == 1]
    
    if (nrow(data_c0) == 0 | nrow(data_c1) == 0) {
      message("  Warning: Missing cancer group for ", analysis_name, " ", score)
      next
    }
    
    n_pos_c0 = data_c0[outcome == 1, sum(n, na.rm = TRUE)]
    n_neg_c0 = data_c0[outcome == 0, sum(n, na.rm = TRUE)]
    n_pos_c1 = data_c1[outcome == 1, sum(n, na.rm = TRUE)]
    n_neg_c1 = data_c1[outcome == 0, sum(n, na.rm = TRUE)]
    
    if (n_pos_c0 == 0 | n_neg_c0 == 0 | n_pos_c1 == 0 | n_neg_c1 == 0) {
      message("  Warning: No outcome variation for ", analysis_name, " ", score)
      next
    }
    
    # Expand to individual observations
    exp_c0 = data_c0[rep(1:.N, n)]
    exp_c1 = data_c1[rep(1:.N, n)]
    
    roc_c0 = roc(exp_c0$outcome, exp_c0[[score_col]], quiet = TRUE, levels = c(0, 1), direction = "<")
    roc_c1 = roc(exp_c1$outcome, exp_c1[[score_col]], quiet = TRUE, levels = c(0, 1), direction = "<")
    
    tst = roc.test(roc_c0, roc_c1, method = "delong", paired = FALSE)
    
    results[[length(results) + 1]] = data.table(
      analysis   = analysis_name,
      score_name = score,
      auroc_c0   = as.numeric(auc(roc_c0)),
      auroc_c1   = as.numeric(auc(roc_c1)),
      diff_auc   = as.numeric(auc(roc_c1) - auc(roc_c0)),
      p_delong   = unname(tst$p.value)
    )
    
    rm(exp_c0, exp_c1)
  }
  
  gc()
  rbindlist(results, use.names = TRUE)
}

#' Calculate AUROC from bootstrap iterations
#' For 24h horizon data where repeated observations violate independence
calculate_auroc_bootstrap = function(boot_data, point_data, analysis_name, score_col = "value") {
  
  if (nrow(boot_data) == 0 | nrow(point_data) == 0) {
    message("  Warning: No data for bootstrap analysis: ", analysis_name)
    return(data.table())
  }
  
  combos = unique(point_data[, .(score_name, ca_01)])
  auroc_results = list()
  
  for (i in seq_len(nrow(combos))) {
    score  = combos$score_name[i]
    cancer = combos$ca_01[i]
    
    # Point estimate from raw counts (for n and n_events)
    point_sub = point_data[score_name == score & ca_01 == cancer]
    n_total  = point_sub[, sum(n, na.rm = TRUE)]
    n_events = point_sub[outcome == 1, sum(n, na.rm = TRUE)]
    
    if (n_total == 0 | n_events == 0) next
    
    # Bootstrap iterations
    boot_sub = boot_data[score_name == score & ca_01 == cancer]
    
    if (nrow(boot_sub) == 0) {
      message("  Warning: No bootstrap data for ", score, " ca_01=", cancer)
      next
    }
    
    # Calculate AUROC for each bootstrap iteration
    boot_aurocs = boot_sub[, {
      # Expand this iteration's counts
      if (sum(n) == 0) return(data.table(auroc = NA_real_))
      
      df_iter = .SD[rep(1:.N, n)]
      
      # Check outcome variation
      if (length(unique(df_iter$outcome)) < 2) return(data.table(auroc = NA_real_))
      
      roc_obj = tryCatch(
        roc(df_iter$outcome, df_iter[[score_col]], quiet = TRUE, levels = c(0, 1), direction = "<"),
        error = function(e) NULL
      )
      
      if (is.null(roc_obj)) return(data.table(auroc = NA_real_))
      
      data.table(auroc = as.numeric(auc(roc_obj)))
    }, by = iter]
    
    # Remove failed iterations
    boot_aurocs = boot_aurocs[!is.na(auroc)]
    
    if (nrow(boot_aurocs) < 10) {
      message("  Warning: Too few valid bootstrap iterations for ", score, " ca_01=", cancer)
      next
    }
    
    # Point estimate (median of bootstrap) and percentile CIs
    auroc_point = median(boot_aurocs$auroc)
    ci_lower    = quantile(boot_aurocs$auroc, 0.025)
    ci_upper    = quantile(boot_aurocs$auroc, 0.975)
    
    auroc_results[[length(auroc_results) + 1]] = data.table(
      analysis   = analysis_name,
      score_name = score,
      ca_01      = cancer,
      n          = n_total,
      n_events   = n_events,
      auroc      = auroc_point,
      ci_lower   = ci_lower,
      ci_upper   = ci_upper,
      n_boot     = nrow(boot_aurocs)
    )
  }
  
  gc()
  rbindlist(auroc_results, use.names = TRUE)
}

#' Bootstrap-based comparison of AUROCs between cancer and non-cancer
#' Returns difference and p-value based on bootstrap distribution
bootstrap_comparison = function(boot_data, point_data, analysis_name, score_col = "value") {
  
  if (nrow(boot_data) == 0) {
    return(data.table())
  }
  
  scores_list = unique(point_data$score_name)
  results = list()
  
  for (score in scores_list) {
    
    boot_sub = boot_data[score_name == score]
    
    # Check both groups exist
    if (!all(c(0, 1) %in% boot_sub$ca_01)) {
      message("  Warning: Missing cancer group in bootstrap for ", score)
      next
    }
    
    # Calculate AUROC for each iteration × cancer group
    iter_aurocs = boot_sub[, {
      if (sum(n) == 0) return(data.table(auroc = NA_real_))
      
      df_iter = .SD[rep(1:.N, n)]
      
      if (length(unique(df_iter$outcome)) < 2) return(data.table(auroc = NA_real_))
      
      roc_obj = tryCatch(
        roc(df_iter$outcome, df_iter[[score_col]], quiet = TRUE, levels = c(0, 1), direction = "<"),
        error = function(e) NULL
      )
      
      if (is.null(roc_obj)) return(data.table(auroc = NA_real_))
      
      data.table(auroc = as.numeric(auc(roc_obj)))
    }, by = .(iter, ca_01)]
    
    # Pivot to get paired differences
    iter_wide = dcast(iter_aurocs, iter ~ ca_01, value.var = "auroc")
    setnames(iter_wide, c("0", "1"), c("auroc_c0", "auroc_c1"))
    iter_wide = iter_wide[!is.na(auroc_c0) & !is.na(auroc_c1)]
    
    if (nrow(iter_wide) < 10) {
      message("  Warning: Too few paired bootstrap iterations for ", score)
      next
    }
    
    iter_wide[, diff := auroc_c1 - auroc_c0]
    
    # Point estimates and CI for difference
    auroc_c0_est = median(iter_wide$auroc_c0)
    auroc_c1_est = median(iter_wide$auroc_c1)
    diff_est     = median(iter_wide$diff)
    
    # Two-sided p-value: proportion of bootstrap diffs on opposite side of 0
    p_boot = 2 * min(
      mean(iter_wide$diff <= 0),
      mean(iter_wide$diff >= 0)
    )
    
    results[[length(results) + 1]] = data.table(
      analysis   = analysis_name,
      score_name = score,
      auroc_c0   = auroc_c0_est,
      auroc_c1   = auroc_c1_est,
      diff_auc   = diff_est,
      p_boot     = p_boot,
      n_boot     = nrow(iter_wide)
    )
  }
  
  gc()
  rbindlist(results, use.names = TRUE)
}

# define analysis variants -----------------------------------------------------

analysis_variants = c("main", "fullcode_only", "no_ed_req", "win0_96h", "one_enc_per_pt")

# ==============================================================================
# ENCOUNTER-LEVEL AUROCs (DeLong - independence assumption holds)
# ==============================================================================

message("\n== Encounter-level AUROCs (DeLong method) ==")

enc_auroc_list  = list()
enc_delong_list = list()

for (variant in analysis_variants) {
  
  df_subset = maxscores_ca_raw[analysis == variant]
  
  n_rows = nrow(df_subset)
  n_enc  = df_subset[, sum(n, na.rm = TRUE)]
  message("  ", variant, ": ", format_n(n_rows), " rows, ~", format_n(n_enc), " encounters")
  
  if (n_rows > 0) {
    enc_auroc_list[[variant]]  = calculate_auroc_delong(df_subset, variant, score_col = "max_value")
    enc_delong_list[[variant]] = delong_comparison(df_subset, variant, score_col = "max_value")
  }
  
  gc()
}

auroc_encounter  = rbindlist(enc_auroc_list,  use.names = TRUE, fill = TRUE)
delong_encounter = rbindlist(enc_delong_list, use.names = TRUE, fill = TRUE)

# Merge DeLong comparison results
if (nrow(delong_encounter) > 0) {
  auroc_encounter = merge(
    auroc_encounter,
    delong_encounter[, .(analysis, score_name, p_delong, auroc_c0, auroc_c1, diff_auc)],
    by = c("analysis", "score_name"),
    all.x = TRUE
  )
}

setorder(auroc_encounter, analysis, score_name, ca_01)

rm(enc_auroc_list, enc_delong_list)
gc()


# ==============================================================================
# 24-HOUR HORIZON AUROCs
# Primary: All observations with DeLong CIs (note: treats repeated obs as independent)
# Sensitivity: Bootstrap with 1 obs/encounter (accounts for within-encounter correlation)
# ==============================================================================

message("\n== 24-hour horizon AUROCs ==")

h24_auroc_list      = list()
h24_comparison_list = list()

# ------------------------------------------------------------------------------
# PRIMARY ANALYSIS: All observations (DeLong CIs)
# Note: CIs are anti-conservative due to within-encounter correlation
# ------------------------------------------------------------------------------

message("\n  -- Primary analysis (all observations, DeLong CIs) --")
message("  Note: CIs treat repeated observations as independent")

for (variant in analysis_variants) {
  
  df_subset = counts_h24_raw[analysis == variant]
  
  n_rows = nrow(df_subset)
  n_obs  = df_subset[, sum(n, na.rm = TRUE)]
  message("  ", variant, ": ", format_n(n_rows), " rows, ~", format_n(n_obs), " observations")
  
  if (n_rows > 0) {
    h24_auroc_list[[variant]]      = calculate_auroc_delong(df_subset, variant, score_col = "value")
    h24_comparison_list[[variant]] = delong_comparison(df_subset, variant, score_col = "value")
  }
  
  gc()
}

auroc_24h_primary      = rbindlist(h24_auroc_list, use.names = TRUE, fill = TRUE)
comparison_24h_primary = rbindlist(h24_comparison_list, use.names = TRUE, fill = TRUE)

# Merge comparison results
if (nrow(comparison_24h_primary) > 0) {
  auroc_24h_primary = merge(
    auroc_24h_primary,
    comparison_24h_primary[, .(analysis, score_name, p_delong, auroc_c0, auroc_c1, diff_auc)],
    by = c("analysis", "score_name"),
    all.x = TRUE
  )
}

auroc_24h_primary[, ci_note := "DeLong (anti-conservative)"]

rm(h24_auroc_list, h24_comparison_list)
gc()

# ------------------------------------------------------------------------------
# SENSITIVITY ANALYSIS: Bootstrap (1 obs/encounter, 400 iterations)
# Accounts for within-encounter correlation
# ------------------------------------------------------------------------------

message("\n  -- Sensitivity analysis (clustered bootstrap, 1 obs/encounter) --")

h24_boot_auroc_list      = list()
h24_boot_comparison_list = list()

# Main analysis with bootstrap
boot_main  = boot_h24_raw[!grepl("se_", site)]
point_main = counts_h24_raw[analysis == "main"]
boot_main[, score_name := gsub("_total$", "", score_name)]

if (nrow(boot_main) > 0 & nrow(point_main) > 0) {
  
  n_iters = uniqueN(boot_main$iter)
  message("  main: ", format_n(n_iters), " bootstrap iterations")
  
  h24_boot_auroc_list[["main"]]      = calculate_auroc_bootstrap(boot_main, point_main, "main", score_col = "value")
  h24_boot_comparison_list[["main"]] = bootstrap_comparison(boot_main, point_main, "main", score_col = "value")
  
} else {
  message("  Warning: No bootstrap data available for main analysis")
}

gc()

auroc_24h_bootstrap      = rbindlist(h24_boot_auroc_list, use.names = TRUE, fill = TRUE)
comparison_24h_bootstrap = rbindlist(h24_boot_comparison_list, use.names = TRUE, fill = TRUE)

# Merge comparison results
if (nrow(comparison_24h_bootstrap) > 0) {
  comparison_24h_bootstrap[, p_delong := p_boot]  # Rename for consistency
  
  auroc_24h_bootstrap = merge(
    auroc_24h_bootstrap,
    comparison_24h_bootstrap[, .(analysis, score_name, p_delong, auroc_c0, auroc_c1, diff_auc)],
    by = c("analysis", "score_name"),
    all.x = TRUE
  )
}

auroc_24h_bootstrap[, ci_note := "Bootstrap (accounts for clustering)"]

rm(h24_boot_auroc_list, h24_boot_comparison_list, boot_main, point_main)
gc()

# ------------------------------------------------------------------------------
# Combine primary and sensitivity 24h results
# ------------------------------------------------------------------------------

# Tag analysis type
auroc_24h_primary[, h24_method := "primary"]
auroc_24h_bootstrap[, h24_method := "sensitivity_bootstrap"]

auroc_24h = rbindlist(list(auroc_24h_primary, auroc_24h_bootstrap), use.names = TRUE, fill = TRUE)

message("\n  24h primary estimates: ", nrow(auroc_24h_primary))
message("  24h bootstrap sensitivity estimates: ", nrow(auroc_24h_bootstrap))

setorder(auroc_24h, h24_method, analysis, score_name, ca_01)

rm(auroc_24h_primary, auroc_24h_bootstrap, comparison_24h_primary, comparison_24h_bootstrap)

# ==============================================================================
# COMBINE RESULTS
# ==============================================================================

message("\n== Combining results ==")

auroc_encounter[, metric := "Encounter max"]
auroc_encounter[, h24_method := NA_character_]

auroc_24h[, metric := "24-hour horizon"]

auroc_all = rbindlist(list(auroc_encounter, auroc_24h), use.names = TRUE, fill = TRUE)

# Add labels
auroc_all[, `:=`(
  ca_lab       = fifelse(ca_01 == 1L, "Cancer", "Non-cancer"),
  analysis_lab = factor(analysis, levels = names(ANALYSIS_LABS), labels = ANALYSIS_LABS),
  score_lab    = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS),
  metric_lab   = factor(metric, levels = c("Encounter max", "24-hour horizon"))
)]

# ==============================================================================
# SUMMARY TABLES
# ==============================================================================

message("\n== Creating summary tables ==")

## main analysis results (primary 24h, not bootstrap sensitivity) --------------

auroc_main = auroc_all[
  analysis == "main" & (is.na(h24_method) | h24_method == "primary"), 
  .(
    metric, score_lab, ca_lab, 
    auroc_fmt = sprintf("%.3f (%.3f-%.3f)", auroc, ci_lower, ci_upper),
    n         = format_n(n),
    n_events  = format_n(n_events)
  )
]

auroc_main_wide = dcast(
  auroc_main, 
  metric + score_lab ~ ca_lab, 
  value.var = "auroc_fmt"
)

## 24h bootstrap sensitivity comparison ----------------------------------------

auroc_24h_comparison = auroc_all[
  metric == "24-hour horizon" & analysis == "main",
  .(h24_method, score_lab, ca_lab, auroc, ci_lower, ci_upper)
]

if (nrow(auroc_24h_comparison) > 0) {
  auroc_24h_comparison[, auroc_fmt := sprintf("%.3f (%.3f-%.3f)", auroc, ci_lower, ci_upper)]
  
  auroc_24h_sens_table = dcast(
    auroc_24h_comparison,
    score_lab + ca_lab ~ h24_method,
    value.var = "auroc_fmt"
  )
  setnames(auroc_24h_sens_table, 
           c("primary", "sensitivity_bootstrap"),
           c("Primary (all obs)", "Sensitivity (bootstrap)"),
           skip_absent = TRUE)
} else {
  auroc_24h_sens_table = data.table()
}

## difference summary ----------------------------------------------------------

auroc_diff = unique(auroc_all[
  analysis == "main" & (is.na(h24_method) | h24_method == "primary") & !is.na(diff_auc), 
  .(metric, score_lab, diff_auc, p_delong)
])

auroc_diff[, `:=`(
  diff_fmt = sprintf("%.3f", diff_auc),
  p_fmt    = fifelse(!is.na(p_delong), sprintf("%.4f", p_delong), "—"),
  sig      = fcase(
    is.na(p_delong), "",
    p_delong < 0.001, "***",
    p_delong < 0.01,  "**",
    p_delong < 0.05,  "*",
    default = ""
  )
)]

## sensitivity analysis comparison (encounter-level) ---------------------------

sens_comparison = auroc_all[ca_01 == 1 & metric == "Encounter max", .(
  analysis_lab, score_lab, auroc
)]

if (nrow(sens_comparison) > 0) {
  sens_comparison = dcast(sens_comparison, score_lab ~ analysis_lab, value.var = "auroc")
}

# ==============================================================================
# EXPORTS
# ==============================================================================

message("\n== Discrimination analysis complete ==")
message("  Encounter-level estimates: ", nrow(auroc_encounter))
message("  24-hour horizon estimates: ", nrow(auroc_24h))
message("    - Primary (all obs): ", nrow(auroc_24h[h24_method == "primary"]))
message("    - Sensitivity (bootstrap): ", nrow(auroc_24h[h24_method == "sensitivity_bootstrap"]))
message("  Total AUROC estimates: ", nrow(auroc_all))

# Objects for figures and other scripts
auroc_results_final    = auroc_all
auroc_main_table       = auroc_main_wide
auroc_diff_table       = auroc_diff
sens_comparison_table  = sens_comparison
auroc_24h_sens_table   = auroc_24h_sens_table  # New: primary vs bootstrap comparison

rm(auroc_encounter, auroc_24h, delong_encounter)
gc()
