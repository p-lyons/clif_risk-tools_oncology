# ==============================================================================
# 05_subgroups.R
# Subgroup analyses: site-level heterogeneity, hematologic vs solid tumors,
# metastatic vs non-metastatic solid tumors, and subgroup threshold performance
#
# Round-two additions, all consuming artifacts that 00_load.R already read and
# that no script referenced:
#
#   mets stratum            auroc_mets_raw, auroc_mets_h24_raw, events_mets_raw
#   liquid thresholds       sesp_liquid_raw, ever_liquid_raw
#   outcome components      events_liquid_raw, events_mets_raw
#
# mets_01 is three-valued by construction: 1 for a solid tumour with a C77-C80
# secondary-site code, 0 for a solid tumour without one, and missing for every
# hematologic and every non-cancer encounter. A hematologic encounter can carry
# a secondary-site code without that meaning metastatic solid disease, so the
# stratum is restricted to solid tumours and the missing level is never
# collapsed to zero. R28 must describe it in exactly those terms.
#
# Discharges R28 (hematologic versus solid, metastatic versus non-metastatic,
# and outcome components by subgroup).
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

# METASTATIC VS NON-METASTATIC SOLID TUMOURS -----------------------------------

message("\n== Metastatic vs non-metastatic solid tumours ==")

#' Meta-analyze site-level AUROCs over a two-level stratum column.
#' Mirrors the liquid block above rather than introducing a second idiom.
meta_by_stratum = function(dt, stratum_col, label) {

  if (is.null(dt) || nrow(dt) == 0) {
    message("  ", label, ": no data; skipping")
    return(data.table())
  }

  d = as.data.table(dt)

  if (!(stratum_col %in% names(d))) {
    message("  ", label, ": no ", stratum_col, " column; skipping")
    return(data.table())
  }

  # mets_01 is three-valued; the missing level is hematologic and non-cancer
  # encounters, which are not part of this comparison.
  d = d[!is.na(get(stratum_col))]

  if ("estimable" %in% names(d)) d = d[estimable == TRUE | is.na(estimable)]
  if (nrow(d) == 0) {
    message("  ", label, ": no estimable rows; skipping")
    return(data.table())
  }

  combos = unique(d[, .(score_name, stratum = get(stratum_col))])
  out    = vector("list", nrow(combos))

  for (i in seq_len(nrow(combos))) {
    sub = d[score_name == combos$score_name[i] & get(stratum_col) == combos$stratum[i]]
    ma  = meta_analyze_aurocs(sub)
    ma[, `:=`(score_name = combos$score_name[i], stratum = combos$stratum[i])]
    out[[i]] = ma
  }

  rbindlist(out, use.names = TRUE, fill = TRUE)
}

if (exists("auroc_mets_raw") && nrow(auroc_mets_raw) > 0) {

  mets_src = as.data.table(auroc_mets_raw)
  if ("outcome_key" %in% names(mets_src)) mets_src = mets_src[outcome_key == "composite"]
  if ("analysis"    %in% names(mets_src)) mets_src = mets_src[analysis == "main"]

  mets_aurocs = meta_by_stratum(mets_src, "mets_01", "mets encounter AUROC")

  if (nrow(mets_aurocs) > 0) {

    mets_aurocs[, `:=`(
      mets_type = fifelse(stratum == 1, "Metastatic", "Non-metastatic"),
      score_lab = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
    )]

    message("  Metastatic vs non-metastatic AUROCs (solid tumours only):")
    print(mets_aurocs[, .(score_lab, mets_type, auroc = round(auroc, 3),
                          ci = paste0("(", round(ci_lower, 3), "-", round(ci_upper, 3), ")"))])

    ## z-test comparison -------------------------------------------------------

    mets_delong = mets_aurocs[, {
      auc_n = auroc[stratum == 0]
      auc_m = auroc[stratum == 1]
      se_n  = auroc_se[stratum == 0]
      se_m  = auroc_se[stratum == 1]

      if (length(auc_n) == 1 && length(auc_m) == 1 &&
          !is.na(auc_n) && !is.na(auc_m) && !is.na(se_n) && !is.na(se_m)) {
        diff = auc_m - auc_n
        se_d = sqrt(se_n^2 + se_m^2)
        p    = 2 * pnorm(-abs(diff / se_d))
        list(auroc_nonmets = auc_n, auroc_mets = auc_m, diff_auc = diff, p_delong = p)
      } else {
        list(auroc_nonmets = NA_real_, auroc_mets = NA_real_,
             diff_auc = NA_real_, p_delong = NA_real_)
      }
    }, by = score_name]

    mets_delong[, score_lab := factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)]

    message("  Metastatic vs non-metastatic comparison:")
    print(mets_delong[, .(score_lab, diff_auc = round(diff_auc, 3),
                          p_delong = round(p_delong, 4))])

  } else {
    mets_delong = data.table()
  }

} else {
  message("  auroc_mets_raw unavailable; skipping metastatic subgroup")
  mets_aurocs = data.table()
  mets_delong = data.table()
}

## 24-hour horizon by metastatic status ---------------------------------------

if (exists("auroc_mets_h24_raw") && nrow(auroc_mets_h24_raw) > 0) {

  mets24_src = as.data.table(auroc_mets_h24_raw)
  if ("outcome_key" %in% names(mets24_src)) mets24_src = mets24_src[outcome_key == "composite"]
  if ("analysis"    %in% names(mets24_src)) mets24_src = mets24_src[analysis == "main"]

  mets_aurocs_24h = meta_by_stratum(mets24_src, "mets_01", "mets 24h AUROC")

  if (nrow(mets_aurocs_24h) > 0) {
    mets_aurocs_24h[, `:=`(
      mets_type = fifelse(stratum == 1, "Metastatic", "Non-metastatic"),
      score_lab = factor(score_name, levels = names(SCORE_LABS), labels = SCORE_LABS)
    )]
    message("  24h horizon by metastatic status calculated (meta-analyzed)")
  }

} else {
  mets_aurocs_24h = data.table()
}

# SUBGROUP THRESHOLD PERFORMANCE -----------------------------------------------
# R28 asks for threshold-specific operating characteristics by cancer subtype,
# not only discrimination. sesp_liquid and ever_liquid carry them at the
# standard threshold. Counts pool by summation, as in 03_threshold.R.

message("\n== Subgroup threshold performance ==")

liquid_sesp = data.table()

if (exists("sesp_liquid_raw") && nrow(sesp_liquid_raw) > 0) {

  ls_src = as.data.table(sesp_liquid_raw)
  if ("outcome_key" %in% names(ls_src)) ls_src = ls_src[outcome_key %in% c("composite", "nohospice")]

  liquid_sesp = ls_src[, .(
    n_total      = sum(n_total),
    n_outcome    = sum(n_outcome),
    n_pos        = sum(n_pos),
    tp           = sum(tp),
    fp           = sum(fp),
    tn           = sum(tn),
    fn           = sum(fn)
  ), by = .(score_name, liquid_01, outcome_key)]

  liquid_sesp[, `:=`(
    sensitivity = tp / (tp + fn),
    specificity = tn / (tn + fp),
    ppv         = tp / (tp + fp),
    npv         = tn / (tn + fn),
    positivity  = n_pos / n_total,
    prevalence  = n_outcome / n_total,
    nne         = fifelse(tp + fp > 0 & tp > 0, (tp + fp) / tp, NA_real_)
  )]

  liquid_sesp[, `:=`(
    cancer_type = fifelse(liquid_01 == 1, "Hematologic", "Solid"),
    score_lab   = factor(score_name,  levels = names(SCORE_LABS),   labels = SCORE_LABS),
    outcome_lab = factor(outcome_key, levels = names(OUTCOME_LABS), labels = OUTCOME_LABS)
  )]

  setorder(liquid_sesp, outcome_lab, score_lab, -liquid_01)

  message("  Subgroup operating characteristics: ", nrow(liquid_sesp), " rows")

} else {
  message("  sesp_liquid_raw unavailable; skipping subgroup threshold performance")
}

# OUTCOME COMPONENTS BY SUBGROUP -----------------------------------------------
# R28 notes that hospice discharge, ICU transfer, and ward death may differ
# substantially across oncology populations, so each component is reported by
# subgroup rather than only the composite. The events artifacts carry n_enc and
# n_events per stratum and outcome.

message("\n== Outcome components by subgroup ==")

pool_events = function(nm, stratum_col, stratum_labs, label) {

  if (!exists(nm) || nrow(get(nm)) == 0) {
    message("  ", label, ": unavailable; skipping")
    return(data.table())
  }

  d = as.data.table(get(nm))

  if (!(stratum_col %in% names(d))) {
    message("  ", label, ": no ", stratum_col, " column; skipping")
    return(data.table())
  }

  d = d[!is.na(get(stratum_col))]

  agg = d[, .(n_enc = sum(n_enc), n_events = sum(n_events)),
          by = c(stratum_col, "outcome_key")]

  agg[, `:=`(
    pct         = n_events / n_enc,
    subgroup    = stratum_labs[as.character(get(stratum_col))],
    outcome_lab = factor(outcome_key, levels = names(OUTCOME_LABS), labels = OUTCOME_LABS)
  )]

  setorder(agg, outcome_lab, subgroup)
  agg
}

events_liquid_pooled = pool_events(
  "events_liquid_raw", "liquid_01",
  c("0" = "Solid", "1" = "Hematologic"),
  "outcome components, hematologic vs solid"
)

events_mets_pooled = pool_events(
  "events_mets_raw", "mets_01",
  c("0" = "Non-metastatic", "1" = "Metastatic"),
  "outcome components, metastatic vs non-metastatic"
)

subgroup_events = rbindlist(
  list(
    if (nrow(events_liquid_pooled) > 0) events_liquid_pooled[, .(subgroup, outcome_lab, n_enc, n_events, pct)] else NULL,
    if (nrow(events_mets_pooled)   > 0) events_mets_pooled[,   .(subgroup, outcome_lab, n_enc, n_events, pct)] else NULL
  ),
  use.names = TRUE
)

if (nrow(subgroup_events) > 0) {
  message("  Outcome components pooled for ", uniqueN(subgroup_events$subgroup), " subgroup(s)")
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

## Metastatic comparison table -------------------------------------------------

if (nrow(mets_aurocs) > 0) {

  mets_type_table = dcast(
    mets_aurocs[, .(score_lab, mets_type,
                    summary = sprintf("%.3f (%.3f-%.3f)", auroc, ci_lower, ci_upper))],
    score_lab ~ mets_type,
    value.var = "summary"
  )

  if (nrow(mets_delong) > 0) {
    mets_type_table = merge(
      mets_type_table,
      mets_delong[, .(score_lab, p_delong = sprintf("%.4f", p_delong))],
      by = "score_lab", all.x = TRUE
    )
  }

} else {
  mets_type_table = data.table()
}

## Subgroup threshold table ----------------------------------------------------

if (nrow(liquid_sesp) > 0) {

  subgroup_threshold_table = liquid_sesp[, .(
    Outcome        = outcome_lab,
    Score          = score_lab,
    Subgroup       = cancer_type,
    N              = format_n(n_total),
    `Outcome rate` = sprintf("%.1f%%", 100 * prevalence),
    `Alert rate`   = sprintf("%.1f%%", 100 * positivity),
    Sensitivity    = sprintf("%.1f%%", 100 * sensitivity),
    Specificity    = sprintf("%.1f%%", 100 * specificity),
    PPV            = sprintf("%.1f%%", 100 * ppv),
    NPV            = sprintf("%.1f%%", 100 * npv),
    `NNE (1/PPV)`  = fifelse(is.na(nne), "—", sprintf("%.1f", nne))
  )]

} else {
  subgroup_threshold_table = data.table()
}

## Outcome components by subgroup ----------------------------------------------

if (nrow(subgroup_events) > 0) {

  subgroup_events_table = dcast(
    subgroup_events[, .(subgroup, outcome_lab,
                        cell = paste0(format_n(n_events), " (",
                                      sprintf("%.1f%%", 100 * pct), ")"))],
    outcome_lab ~ subgroup,
    value.var = "cell",
    fill = "—"
  )

  setnames(subgroup_events_table, "outcome_lab", "Outcome")

} else {
  subgroup_events_table = data.table()
}

## write the new tables --------------------------------------------------------

if (requireNamespace("flextable", quietly = TRUE)) {

  sub_path = function(stem) {
    if (exists("today")) {
      here("output", "tables", paste0(stem, "_", today, ".docx"))
    } else {
      here("output", "tables", paste0(stem, ".docx"))
    }
  }

  write_subgroup = function(dt, stem, label) {
    if (is.null(dt) || nrow(dt) == 0) {
      message("  ", label, ": nothing to write")
      return(invisible(NULL))
    }
    ft = flextable::flextable(as.data.frame(dt)) |>
      flextable::autofit() |>
      flextable::bold(part = "header")
    p = sub_path(stem)
    flextable::save_as_docx(ft, path = p)
    message("  Saved: ", p, " (", nrow(dt), " rows)")
  }

  write_subgroup(mets_type_table,          "table_s_mets_auroc",         "Metastatic AUROC")
  write_subgroup(subgroup_threshold_table, "table_s_subgroup_threshold", "Subgroup thresholds")
  write_subgroup(subgroup_events_table,    "table_s_subgroup_events",    "Subgroup outcome components")
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

# round-two additions
mets_aurocs_final              = mets_aurocs
mets_delong_final              = mets_delong
mets_aurocs_24h_final          = mets_aurocs_24h
mets_type_table_final          = mets_type_table
liquid_sesp_final              = liquid_sesp
subgroup_threshold_table_final = subgroup_threshold_table
subgroup_events_final          = subgroup_events
subgroup_events_table_final    = subgroup_events_table
