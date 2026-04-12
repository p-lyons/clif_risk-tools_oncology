# ==============================================================================
# 05_subgroups.R
# Subgroup analyses: site-level heterogeneity, hematologic vs solid tumors
# ==============================================================================

# setup ------------------------------------------------------------------------

library(pROC)
library(ggplot2)

# SITE-LEVEL AUROCs ------------------------------------------------------------

message("\n== Site-level AUROCs ==")

## use pre-computed site-level AUROCs from auroc_enc_raw -----------------------
## these were calculated with proper DeLong CIs at each site in 03_analysis.R

site_aurocs = auroc_enc_raw[analysis == "main", .(
  site, score_name, ca_01,
  n        = n_obs,
  n_events,
  auroc,
  auroc_se,
  ci_lower,
  ci_upper
)]

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

if (nrow(auroc_liquid_enc) > 0) {
  
  ## use pre-computed site-level liquid/solid AUROCs ----------------------------
  ## meta-analyze across sites for each score × liquid_01 combination
  
  liquid_combos = unique(auroc_liquid_enc[, .(score_name, liquid_01)])
  
  liquid_meta_list = lapply(seq_len(nrow(liquid_combos)), function(i) {
    sub = auroc_liquid_enc[
      score_name == liquid_combos$score_name[i] & 
        liquid_01  == liquid_combos$liquid_01[i]
    ]
    ma = meta_analyze_aurocs(sub)
    ma[, `:=`(score_name = liquid_combos$score_name[i], liquid_01 = liquid_combos$liquid_01[i])]
  })
  
  liquid_aurocs = rbindlist(liquid_meta_list, use.names = TRUE, fill = TRUE)
  
  liquid_aurocs[, `:=`(
    cancer_type = fifelse(liquid_01 == 1, "Hematologic", "Solid"),
    score_lab   = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
  )]
  
  message("  Hematologic vs solid AUROCs:")
  print(liquid_aurocs[, .(score_lab, cancer_type, auroc = round(auroc, 3), 
                          ci = paste0("(", round(ci_lower, 3), "-", round(ci_upper, 3), ")"))])
  
  ## z-test comparison between liquid and solid --------------------------------
  
  liquid_delong = liquid_aurocs[, {
    auc_s  = auroc[liquid_01 == 0]
    auc_l  = auroc[liquid_01 == 1]
    se_s   = auroc_se[liquid_01 == 0]
    se_l   = auroc_se[liquid_01 == 1]
    
    if (length(auc_s) == 1 && length(auc_l) == 1 &&
        !is.na(auc_s) && !is.na(auc_l) && !is.na(se_s) && !is.na(se_l)) {
      diff   = auc_l - auc_s
      se_d   = sqrt(se_s^2 + se_l^2)
      z      = diff / se_d
      p      = 2 * pnorm(-abs(z))
      list(auroc_solid = auc_s, auroc_liquid = auc_l, diff_auc = diff, p_delong = p)
    } else {
      list(auroc_solid = NA_real_, auroc_liquid = NA_real_, 
           diff_auc = NA_real_, p_delong = NA_real_)
    }
  }, by = score_name]
  
  liquid_delong[, score_lab := factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)]
  
  message("  Liquid vs solid comparison (z-test on meta-analyzed AUROCs):")
  print(liquid_delong[, .(score_lab, diff_auc = round(diff_auc, 3), p_delong = round(p_delong, 4))])
  
} else if (nrow(maxscores_liquid_raw) > 0) {
  message("  WARNING: No site-level liquid AUROCs; skipping subgroup")
  liquid_aurocs = data.table()
  liquid_delong = data.table()
} else {
  message("  WARNING: No liquid/solid data available")
  liquid_aurocs = data.table()
  liquid_delong = data.table()
}

# 24-HOUR HORIZON BY CANCER TYPE -----------------------------------------------

message("\n== 24-hour AUROCs by cancer type ==")

if (nrow(auroc_liquid_h24) > 0) {
  
  ## meta-analyze site-level 24h liquid/solid AUROCs ---------------------------
  
  liq24_combos = unique(auroc_liquid_h24[, .(score_name, liquid_01)])
  
  liq24_meta_list = lapply(seq_len(nrow(liq24_combos)), function(i) {
    sub = auroc_liquid_h24[
      score_name == liq24_combos$score_name[i] & 
        liquid_01  == liq24_combos$liquid_01[i]
    ]
    ma = meta_analyze_aurocs(sub)
    ma[, `:=`(score_name = liq24_combos$score_name[i], liquid_01 = liq24_combos$liquid_01[i])]
  })
  
  liquid_aurocs_24h = rbindlist(liq24_meta_list, use.names = TRUE, fill = TRUE)
  
  liquid_aurocs_24h[, `:=`(
    cancer_type = fifelse(liquid_01 == 1, "Hematologic", "Solid"),
    score_lab   = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
  )]
  
  message("  24h horizon by cancer type calculated (meta-analyzed)")
  
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
