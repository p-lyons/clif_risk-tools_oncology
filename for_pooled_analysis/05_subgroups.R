# ==============================================================================
# 05_subgroups.R
# Subgroup analyses: site-level heterogeneity, hematologic vs solid tumors
# ==============================================================================

# setup ------------------------------------------------------------------------

library(pROC)
library(ggplot2)

source(here::here("code_pooled", "00_load.R"))

# SITE-LEVEL AUROCs ------------------------------------------------------------

message("\n== Site-level AUROCs ==")

## calculate AUROC for each site × score × cancer combination ------------------

df_main = maxscores_ca_raw[analysis == "main"]

# Expand aggregated data
df_expanded = df_main[rep(1:.N, times = n), .(site, score_name, ca_01, outcome, max_value)]

# Calculate AUROCs
site_aurocs = df_expanded[, {
  roc_obj = roc(outcome, max_value, levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)
  ci_vals = ci.auc(roc_obj)
  list(
    n        = .N,
    n_events = sum(outcome),
    auroc    = as.numeric(auc(roc_obj)),
    ci_lower = as.numeric(ci_vals[1]),
    ci_upper = as.numeric(ci_vals[3])
  )
}, by = .(site, score_name, ca_01)]

message("  Calculated AUROCs for ", uniqueN(site_aurocs$site), " sites")

## add labels ------------------------------------------------------------------

site_aurocs[, `:=`(
  ca_lab    = fifelse(ca_01 == 1, "Cancer", "Non-cancer"),
  score_lab = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
)]

## rank sites by performance ---------------------------------------------------

# Rank based on cancer patient AUROCs
site_rankings = site_aurocs[ca_01 == 1, .(
  site, score_name, auroc
)][, rank := rank(-auroc), by = score_name]

# Mean rank across all scores
site_order = site_rankings[, .(mean_rank = mean(rank)), by = site][order(mean_rank)]

site_aurocs[, site_ordered := factor(site, levels = site_order$site)]

## calculate heterogeneity metrics ---------------------------------------------

site_heterogeneity = site_aurocs[, .(
  n_sites    = .N,
  mean_auroc = mean(auroc),
  sd_auroc   = sd(auroc),
  min_auroc  = min(auroc),
  max_auroc  = max(auroc),
  range      = max(auroc) - min(auroc)
), by = .(score_name, ca_01)]

site_heterogeneity[, `:=`(
  ca_lab    = fifelse(ca_01 == 1, "Cancer", "Non-cancer"),
  score_lab = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
)]

message("  Heterogeneity summary:")
print(site_heterogeneity[, .(score_lab, ca_lab, mean_auroc = round(mean_auroc, 3), 
                              sd_auroc = round(sd_auroc, 3), range = round(range, 3))])

# HEMATOLOGIC VS SOLID TUMOR SUBGROUP ------------------------------------------

message("\n== Hematologic vs solid tumor subgroup ==")

if (nrow(maxscores_liquid_raw) > 0) {
  
  df_liquid = maxscores_liquid_raw
  
  # Expand and calculate AUROCs
  df_liquid_exp = df_liquid[rep(1:.N, times = n), .(score_name, liquid_01, outcome, max_value)]
  
  liquid_aurocs = df_liquid_exp[, {
    roc_obj = roc(outcome, max_value, levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)
    ci_vals = ci.auc(roc_obj)
    list(
      n        = .N,
      n_events = sum(outcome),
      auroc    = as.numeric(auc(roc_obj)),
      ci_lower = as.numeric(ci_vals[1]),
      ci_upper = as.numeric(ci_vals[3])
    )
  }, by = .(score_name, liquid_01)]
  
  liquid_aurocs[, `:=`(
    cancer_type = fifelse(liquid_01 == 1, "Hematologic", "Solid"),
    score_lab   = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
  )]
  
  message("  Hematologic vs solid AUROCs:")
  print(liquid_aurocs[, .(score_lab, cancer_type, auroc = round(auroc, 3), 
                          ci = paste0("(", round(ci_lower, 3), "-", round(ci_upper, 3), ")"))])
  
  ## DeLong comparison between liquid and solid --------------------------------
  
  liquid_delong = df_liquid_exp[, {
    if (sum(liquid_01 == 0) > 0 & sum(liquid_01 == 1) > 0 &
        length(unique(outcome[liquid_01 == 0])) == 2 &
        length(unique(outcome[liquid_01 == 1])) == 2) {
      
      roc_solid  = roc(outcome[liquid_01 == 0], max_value[liquid_01 == 0], quiet = TRUE)
      roc_liquid = roc(outcome[liquid_01 == 1], max_value[liquid_01 == 1], quiet = TRUE)
      tst        = roc.test(roc_solid, roc_liquid, method = "delong", paired = FALSE)
      
      list(
        auroc_solid   = as.numeric(auc(roc_solid)),
        auroc_liquid  = as.numeric(auc(roc_liquid)),
        diff_auc      = as.numeric(auc(roc_liquid) - auc(roc_solid)),
        p_delong      = unname(tst$p.value)
      )
    } else {
      list(auroc_solid = NA_real_, auroc_liquid = NA_real_, 
           diff_auc = NA_real_, p_delong = NA_real_)
    }
  }, by = score_name]
  
  liquid_delong[, score_lab := factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)]
  
  message("  DeLong comparison (liquid - solid):")
  print(liquid_delong[, .(score_lab, diff_auc = round(diff_auc, 3), p_delong = round(p_delong, 4))])
  
} else {
  message("  WARNING: No liquid/solid data available")
  liquid_aurocs = data.table()
  liquid_delong = data.table()
}

# 24-HOUR HORIZON BY CANCER TYPE -----------------------------------------------

message("\n== 24-hour AUROCs by cancer type ==")

if (nrow(counts_liquid_raw) > 0) {
  
  df_liquid_24h = counts_liquid_raw
  df_liquid_24h[, score_name := str_remove(score_name, "_total")]
  
  # Expand and calculate
  df_liq24_exp = df_liquid_24h[rep(1:.N, times = n), .(score_name, liquid_01, outcome, value)]
  
  liquid_aurocs_24h = df_liq24_exp[, {
    roc_obj = roc(outcome, value, levels = c(0, 1), direction = "<", ci = TRUE, quiet = TRUE)
    ci_vals = ci.auc(roc_obj)
    list(
      n        = .N,
      n_events = sum(outcome),
      auroc    = as.numeric(auc(roc_obj)),
      ci_lower = as.numeric(ci_vals[1]),
      ci_upper = as.numeric(ci_vals[3])
    )
  }, by = .(score_name, liquid_01)]
  
  liquid_aurocs_24h[, `:=`(
    cancer_type = fifelse(liquid_01 == 1, "Hematologic", "Solid"),
    score_lab   = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
  )]
  
  message("  24h horizon by cancer type calculated")
  
} else {
  liquid_aurocs_24h = data.table()
}

# SUMMARY TABLES ---------------------------------------------------------------

message("\n== Creating summary tables ==")

## Site heterogeneity table ----------------------------------------------------

site_het_table = dcast(
  site_heterogeneity[, .(score_lab, ca_lab, 
                         summary = sprintf("%.3f (%.3f-%.3f)", mean_auroc, min_auroc, max_auroc))],
  score_lab ~ ca_lab,
  value.var = "summary"
)

## Cancer type comparison table ------------------------------------------------

if (nrow(liquid_aurocs) > 0) {
  cancer_type_table = dcast(
    liquid_aurocs[, .(score_lab, cancer_type,
                      summary = sprintf("%.3f (%.3f-%.3f)", auroc, ci_lower, ci_upper))],
    score_lab ~ cancer_type,
    value.var = "summary"
  )
  
  # Add p-values
  cancer_type_table = merge(
    cancer_type_table,
    liquid_delong[, .(score_lab, p_delong = sprintf("%.4f", p_delong))],
    by = "score_lab"
  )
} else {
  cancer_type_table = data.table()
}

# EXPORTS ----------------------------------------------------------------------

message("\n== Subgroup analysis complete ==")

# Objects for figures
site_aurocs_final       = site_aurocs
site_heterogeneity_final = site_heterogeneity
site_rankings_final     = site_rankings
site_order_final        = site_order
liquid_aurocs_final     = liquid_aurocs
liquid_delong_final     = liquid_delong
liquid_aurocs_24h_final = liquid_aurocs_24h
site_het_table_final    = site_het_table
cancer_type_table_final = cancer_type_table
