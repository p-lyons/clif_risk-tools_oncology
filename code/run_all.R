# run_all.R — orchestrate the round-two site pipeline.
# Order: 00 setup, 01 cohort, 02 scores, 02b carry-forward, 02c monitoring,
#        03a artifacts, 03b lead time, 03c GEE models. Then a manifest of every
#        uploaded file for the coordinating center to reconcile.
#
# Orchestration note. Each stage is sourced into the global environment (so a
# stage can see the objects the previous stage left behind — cohort, scores,
# crosswalks, etc.). Several stages END with a whitelist cleanup of the form
#   rm(list = setdiff(ls(), keep)); gc()
# which, because the stage runs in the global environment, deletes every global
# object not on that stage's keep-list. Those keep-lists deliberately preserve
# `run_log` (and nothing else of this script's), so this orchestrator must keep
# ALL of its own state inside `run_log` and must not rely on any other top-level
# object or helper function surviving across a source() call. That is why the
# stage loop below stores its per-stage start time in run_log rather than in a
# local, and why PIPELINE_VERSION and the manifest helpers are defined AFTER the
# stages have run.

run_log = list()
run_log$t_start_num = as.numeric(Sys.time())   # wall-clock start (epoch seconds)

# stages, in order -------------------------------------------------------------
# Each entry is c(label, path). The loop is a single top-level for(): R captures
# the sequence internally, so a stage wiping the global `stages` binding mid-loop
# does not interrupt iteration.

stages = list(
  c("00",  "code/00_setup.R"),
  c("01",  "code/01_cohort.R"),
  c("02",  "code/02_scores.R"),
  c("02b", "code/02b_carryforward.R"),
  c("02c", "code/02c_monitoring.R"),
  c("03a", "code/03a_artifacts.R"),
  c("03b", "code/03b_leadtime.R"),
  c("03c", "code/03c_models.R")
)

message("Starting round-two pipeline...")

for (st in stages) {

  # Stash the label and start time in run_log BEFORE sourcing: the stage may wipe
  # every other global (including `st` itself), but run_log is always preserved.
  run_log[["._lab"]]    = st[1]
  run_log[["._t0_num"]] = as.numeric(Sys.time())

  gc(reset = TRUE)              # reset the max-used counters for this stage
  source(st[2])                 # evaluated in .GlobalEnv (source default local = FALSE)

  # Everything referenced from here to the next iteration must be either created
  # now (mins, peak) or preserved by the stage (run_log, and objects on the
  # stage's keep-list). Do NOT reference `st` — it may have been wiped above.
  lab  = run_log[["._lab"]]
  mins = (as.numeric(Sys.time()) - run_log[["._t0_num"]]) / 60
  peak = sum(gc()[, 6])         # column 6 = "max used (Mb)", summed over N/Vcells

  run_log[[paste0(lab, "_min")]]     = mins
  run_log[[paste0(lab, "_peak_mb")]] = peak
  message(sprintf("  [%-4s] %6.1f min | peak R mem %6.0f Mb", lab, mins, peak))

  # per-stage summary captures (cohort is preserved by 01; .n_score_rows by 02) -
  if (lab == "01") {
    run_log$n_cohort     = nrow(cohort)
    run_log$n_cancer     = sum(cohort$ca_01 == 1)
    run_log$n_ed_admit   = sum(cohort$ed_admit_01 == 1)
    # Composite outcome rate over ED-admit encounters, from the authoritative
    # o_composite_01 flag (ward-to-ICU transfer OR ward death OR hospice
    # discharge). The prior mean(dead_01 + hospice_01 > 0) omitted ward-to-ICU
    # transfers entirely, so it understated the composite rate.
    .ot      = as.data.table(read_parquet(here("proj_tables", "outcome_times.parquet")))
    .ed_jids = cohort$joined_hosp_id[cohort$ed_admit_01 == 1]
    run_log$outcome_rate = mean(.ot[joined_hosp_id %in% .ed_jids]$o_composite_01)
    rm(.ot, .ed_jids)
  }
  if (lab == "02") {
    run_log$n_score_rows = .n_score_rows
  }
}

run_log$total_minutes = (as.numeric(Sys.time()) - run_log$t_start_num) / 60

# drop the loop scratch keys so they do not linger in run_log ------------------
run_log[["._lab"]]    = NULL
run_log[["._t0_num"]] = NULL

# summary report ---------------------------------------------------------------

report = tidytable(
  site             = site_lowercase,
  n_encounters     = run_log$n_cohort,
  n_cancer         = run_log$n_cancer,
  n_ed_admit       = run_log$n_ed_admit,
  pct_cancer       = round(100 * run_log$n_cancer / run_log$n_cohort, 1),
  outcome_rate_pct = round(100 * run_log$outcome_rate, 1),
  n_score_rows     = run_log$n_score_rows,
  runtime_min      = round(run_log$total_minutes, 1),
  completed        = Sys.time()
)

print(report)
fwrite(report, here(BOX_DIR, paste0("run_report_", site_lowercase, ".csv")))

# upload manifest --------------------------------------------------------------
# Walk upload_to_box_v2/ and record every CSV with its dimensions and MD5, so the
# coordinating center can confirm a complete, uncorrupted upload before ingestion.
# Defined here, after the stages, so the stage cleanups cannot wipe it first.

PIPELINE_VERSION = "v2"   # bump when the artifact set changes; could be a git SHA

build_manifest = function(box, site, pipeline_version, run_date) {

  rels = list.files(box, pattern = "\\.csv$", recursive = TRUE, full.names = FALSE)
  rels = rels[!grepl("^MANIFEST-", basename(rels))]   # never list the manifest itself

  rbindlist(lapply(rels, function(rel) {
    fp = file.path(box, rel)
    dt = tryCatch(fread(fp, showProgress = FALSE), error = function(e) NULL)
    data.table::data.table(
      relative_path = rel,
      n_rows        = if (is.null(dt)) NA_integer_ else nrow(dt),
      n_cols        = if (is.null(dt)) NA_integer_ else ncol(dt),
      md5           = unname(tools::md5sum(fp))
    )
  }))[, `:=`(site = site, run_date = run_date, pipeline_version = pipeline_version)][]
}

# Representative expected artifacts — one per family/stage. A missing family means
# a stage failed silently. The coordinating center (P9) reconciles the full list.
EXPECTED_PATTERNS = c(
  "^admission_dx_chapter-ca-",
  "^admission_dx_stem-ca-",
  "^main/maxscores-ca-composite-",
  "^main/maxscores-mets-composite-",
  "^main/auroc-ca-composite-",
  "^main/firstscore-ca-",
  "^main/events-dxgroup-composite-",
  "^threshold/ever-ca-composite-",
  "^threshold/sesp-ca-composite-",
  "^threshold/cuminc-ca-composite-",
  "^threshold/first-ca-composite-",
  "^threshold/upset-ca-composite-",
  "^threshold/leadtime-ca-composite-",
  "^threshold/crossclass-ca-composite-",
  "^horizon/counts-ca-composite-h24-",
  "^sensitivity/maxscores-ca-composite-se_",
  "^sensitivity/maxscores-ca-cf6-",
  "^sensitivity/counts-ca-h24-cf6-",
  "^diagnostics/news_o2_resolution-",
  "^diagnostics/monitoring-ca-",
  "^diagnostics/monitoring_bins-ca-",
  "^diagnostics/missing_vlab-ca-",
  "^meta/coefficients-",
  "^meta/score_sds-",
  "^meta/gee_coefficients-"
)

manifest = build_manifest(here(BOX_DIR), site_lowercase, PIPELINE_VERSION, as.character(Sys.Date()))
fwrite(manifest, here(BOX_DIR, paste0("MANIFEST-", site_lowercase, ".csv")))

missing_families = EXPECTED_PATTERNS[
  !vapply(EXPECTED_PATTERNS, function(p) any(grepl(p, manifest$relative_path)), logical(1))
]
if (length(missing_families) > 0) {
  stop("MANIFEST QC failed — expected artifact families absent:\n  ",
       paste(missing_families, collapse = "\n  "), call. = FALSE)
}
message("✅ MANIFEST QC: all ", length(EXPECTED_PATTERNS), " expected artifact families present (",
        nrow(manifest), " files listed).")

message("Pipeline completed in ", round(run_log$total_minutes, 1), " minutes")
