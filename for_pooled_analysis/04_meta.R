# ==============================================================================
# 04_meta.R
# Threshold-equivalence driver.
#
# The analysis itself is unchanged and lives in 04_meta_body.R. This script runs
# it once per outcome rather than once for the composite.
#
# Round-one structure: a single linear script reading maxscores_ca_raw and
# sesp_raw, both of which legacy_view() had already filtered to the composite.
# Its closing note anticipated exactly this change — "if component outcomes are
# exported separately in a future run, we can rerun these analyses for the
# ICU-transfer-or-death subcomposite as a planned sensitivity analysis." The
# site pipeline now exports all five outcomes, so that sensitivity analysis is
# no longer hypothetical and R27 depends on it.
#
# The body is sourced with three objects set:
#   MX_IN       maxscores rows for one outcome
#   SESP_IN     sesp rows for the same outcome
#   OUT_SUFFIX  "" for the composite, "_<outcome>" otherwise
#
# Composite output is byte-identical to round one: it runs first, writes to the
# unsuffixed filenames, and its *_final objects are restored under their
# round-one names after the loop so 06_figures.R is untouched. Per-outcome
# results travel in thr_eq_by_outcome and its siblings.
#
# Discharges the threshold half of R27 (does the cancer gap in operating
# characteristics persist without hospice discharge) and of R28.
# ==============================================================================

# setup ------------------------------------------------------------------------

library(here)

## outcomes to run -------------------------------------------------------------
# Composite first, so its unsuffixed exports are written before any other
# outcome can overwrite a shared object name. The remaining four follow in the
# order they appear in OUTCOME_LABS.

META_OUTCOMES = c("composite", "nohospice", "wardicu", "warddeath", "hospicedc")

## guard -----------------------------------------------------------------------

if (!exists("maxscores_ca_all_raw") || nrow(maxscores_ca_all_raw) == 0) {
  stop("04_meta.R: maxscores_ca_all_raw is empty; run 00_load.R first.", call. = FALSE)
}
if (!exists("sesp_all_raw") || nrow(sesp_all_raw) == 0) {
  stop("04_meta.R: sesp_all_raw is empty; run 00_load.R first.", call. = FALSE)
}

mx_all   = as.data.table(maxscores_ca_all_raw)
sesp_all = as.data.table(sesp_all_raw)

available = intersect(META_OUTCOMES, unique(mx_all$outcome_key))
missing   = setdiff(META_OUTCOMES, available)

if (length(missing) > 0) {
  message("  04_meta: no maxscores rows for outcome(s) ",
          paste(missing, collapse = ", "), "; skipping those")
}
if (!("composite" %in% available)) {
  stop("04_meta.R: the composite outcome is absent; cannot proceed.", call. = FALSE)
}

## collectors ------------------------------------------------------------------
# The body assigns a fixed set of *_final names on every pass. They are harvested
# into per-outcome lists immediately after each source() call, before the next
# iteration overwrites them.

HARVEST = c(
  "rd_counts_final", "rd_wide_final", "rd_trend_final", "rd_table_final",
  "sweep_final", "thr_eq_final", "thr_eq_table_final", "eff_data_final",
  "site_ops_final", "news_fc_operating_final",
  "site_sweep_final", "het_diag_final"
)

meta_results = list()

# run --------------------------------------------------------------------------

for (ok in available) {

  message("\n══════════════════════════════════════════════════════════════════════")
  message("  04_meta: outcome = ", ok)
  message("══════════════════════════════════════════════════════════════════════")

  MX_IN      = mx_all[outcome_key == ok]
  SESP_IN    = sesp_all[outcome_key == ok]
  OUT_SUFFIX = if (ok == "composite") "" else paste0("_", ok)

  if (nrow(MX_IN) == 0) {
    message("  no maxscores rows; skipping")
    next
  }
  if (nrow(SESP_IN) == 0) {
    message("  no sesp rows; the per-site heterogeneity check will be skipped")
  }

  source(here::here("run_pooled_analysis", "04_meta_body.R"), local = FALSE)

  harvested = list()

  for (nm in HARVEST) {
    if (exists(nm)) {
      obj = get(nm)
      if (is.data.frame(obj)) {
        obj = as.data.table(copy(obj))
        obj[, outcome_key := ok]
      }
      harvested[[nm]] = obj
    }
  }

  meta_results[[ok]] = harvested

  message("  outcome '", ok, "' complete: ",
          nrow(harvested$thr_eq_final), " threshold-equivalence rows")
}

# ==============================================================================
# COMBINE AND RESTORE
# ==============================================================================

message("\n== 04_meta: combining across outcomes ==")

## stack each harvested object across outcomes ---------------------------------

stack_across = function(nm) {
  parts = lapply(meta_results, function(x) x[[nm]])
  parts = parts[!vapply(parts, is.null, logical(1))]
  parts = parts[vapply(parts, is.data.frame, logical(1))]
  if (length(parts) == 0) return(data.table())
  rbindlist(parts, use.names = TRUE, fill = TRUE)
}

thr_eq_by_outcome       = stack_across("thr_eq_final")
thr_eq_table_by_outcome = stack_across("thr_eq_table_final")
sweep_by_outcome        = stack_across("sweep_final")
eff_data_by_outcome     = stack_across("eff_data_final")
rd_table_by_outcome     = stack_across("rd_table_final")
rd_trend_by_outcome     = stack_across("rd_trend_final")
het_diag_by_outcome     = stack_across("het_diag_final")

message("  Threshold-equivalence rows across outcomes: ", nrow(thr_eq_by_outcome))
message("  Efficiency-curve points across outcomes: ",   nrow(eff_data_by_outcome))

## restore the composite objects under their round-one names -------------------
# 06_figures.R consumes sweep_final, eff_data_final, thr_eq_final, site_ops_final,
# news_fc_operating_final, site_sweep_final, and het_diag_final, and assumes one
# outcome per row. The loop leaves whichever outcome ran last in those names, so
# they are explicitly reset to the composite here. Any figure that wants another
# outcome should read the *_by_outcome objects and filter.

composite_objs = meta_results[["composite"]]

for (nm in HARVEST) {
  if (!is.null(composite_objs[[nm]])) {
    assign(nm, composite_objs[[nm]])
  }
}

## drop the loop-scoped inputs so a later script cannot read a stale slice ------

rm(MX_IN, SESP_IN, OUT_SUFFIX)

message("\n== 04_meta complete ==")
message("  Composite objects restored under their round-one names")
message("  Per-outcome results in thr_eq_by_outcome, sweep_by_outcome, ",
        "eff_data_by_outcome, rd_table_by_outcome, het_diag_by_outcome")

gc()
