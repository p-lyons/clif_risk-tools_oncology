# run_all.R

start_time = Sys.time()
log        = list()

message("Starting data import/validation...")
source("code/00_setup.R")
log$setup_time = Sys.time()

message("Cleaning data...")
source("code/01_cohort.R")
log$cohort_time = Sys.time()
log$n_cohort    = nrow(cohort)
log$n_cancer    = sum(cohort$ca_01 == 1)
log$n_ed_admit  = sum(cohort$ed_admit_01 == 1)
log$outcome_rate = mean(cohort$dead_01 + cohort$hospice_01 > 0)

message("Calculating scores...")
source("code/02_scores.R")
log$scores_time  = Sys.time()
log$n_score_rows = nrow(scores)

message("Performing analysis...")
source("code/03_analysis.R")
log$analysis_time = Sys.time()

# summary report
log$total_minutes = as.numeric(difftime(Sys.time(), start_time, units = "mins"))

report = tidytable(
  site             = site_lowercase,
  n_encounters     = log$n_cohort,
  n_cancer         = log$n_cancer,
  n_ed_admit       = log$n_ed_admit,
  pct_cancer       = round(100 * log$n_cancer / log$n_cohort, 1),
  outcome_rate_pct = round(100 * log$outcome_rate, 1),
  n_score_rows     = log$n_score_rows,
  runtime_min      = round(log$total_minutes, 1),
  completed        = Sys.time()
)

print(report)
fwrite(report, here("upload_to_box", paste0("run_report_", site_lowercase, ".csv")))

message("Pipeline completed in ", round(log$total_minutes, 1), " minutes")
