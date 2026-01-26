# ==============================================================================
# run_pooled.R
# Orchestration script for pooled analyses
# Run this after all sites have uploaded their aggregated data to Box
# ==============================================================================

# configuration ----------------------------------------------------------------

# Set to TRUE to run individual sections, FALSE to skip
RUN_TABLES         = TRUE
RUN_DISCRIMINATION = TRUE
RUN_THRESHOLD      = TRUE
RUN_META           = TRUE
RUN_SUBGROUPS      = TRUE
RUN_FIGURES        = TRUE

# Set working directory to project root (uncomment if needed)
# setwd(here::here())

# setup ------------------------------------------------------------------------

start_time = Sys.time()

message("\n")
message("╔══════════════════════════════════════════════════════════════════════╗")
message("║                    CLIF Oncology Risk Project                        ║")
message("║                       Pooled Analysis Pipeline                       ║")
message("╚══════════════════════════════════════════════════════════════════════╝")
message("\n")
message("Started: ", format(start_time, "%Y-%m-%d %H:%M:%S"))
message("")

# check dependencies -----------------------------------------------------------

message("== Checking dependencies ==\n")

required_packages = c(
  "data.table", "tidytable", "collapse", "stringr", "here",
  "pROC", "metafor", "glmmTMB",
  "ggplot2", "patchwork", "scales",
  "flextable", "officer",
  "DiagrammeR", "DiagrammeRsvg", "rsvg",
  "UpSetR"
)

missing = required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing) > 0) {
  stop(
    "Missing required packages: ", paste(missing, collapse = ", "), 
    "\nInstall with: install.packages(c('", paste(missing, collapse = "', '"), "'))",
    call. = FALSE
  )
}

message("  All ", length(required_packages), " required packages available\n")

# create output directories ----------------------------------------------------

message("== Setting up output directories ==\n")

dirs_needed = c(

  here::here("output"),
  here::here("output", "figures"),
  here::here("output", "tables")
)

for (d in dirs_needed) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    message("  Created: ", d)
  }
}

message("  Output directories ready\n")

# load data (00_load.R) --------------------------------------------------------

message("== Loading data ==\n")

source(here::here("code_pooled", "00_load.R"))

message("")

# run analyses -----------------------------------------------------------------

run_log = list()

## 01: tables ------------------------------------------------------------------

if (RUN_TABLES) {
  
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  message("  01_tables.R")
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  
  t1 = Sys.time()
  
  tryCatch({
    source(here::here("code_pooled", "01_tables.R"))
    run_log$tables = list(status = "success", time = Sys.time() - t1)
    message("\n  ✓ 01_tables.R complete (", round(difftime(Sys.time(), t1, units = "secs"), 1), " sec)\n")
  }, error = function(e) {
    run_log$tables <<- list(status = "failed", error = conditionMessage(e))
    message("\n  ✗ 01_tables.R FAILED: ", conditionMessage(e), "\n")
  })
  
} else {
  message("\n  ○ 01_tables.R skipped\n")
}

## 02: discrimination ----------------------------------------------------------

if (RUN_DISCRIMINATION) {
  
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  message("  02_discrimination.R")
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  
  t1 = Sys.time()
  
  tryCatch({
    source(here::here("code_pooled", "02_discrimination.R"))
    run_log$discrimination = list(status = "success", time = Sys.time() - t1)
    message("\n  ✓ 02_discrimination.R complete (", round(difftime(Sys.time(), t1, units = "secs"), 1), " sec)\n")
  }, error = function(e) {
    run_log$discrimination <<- list(status = "failed", error = conditionMessage(e))
    message("\n  ✗ 02_discrimination.R FAILED: ", conditionMessage(e), "\n")
  })
  
} else {
  message("\n  ○ 02_discrimination.R skipped\n")
}

## 03: threshold ---------------------------------------------------------------

if (RUN_THRESHOLD) {
  
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  message("  03_threshold.R")
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  
  t1 = Sys.time()
  
  tryCatch({
    source(here::here("code_pooled", "03_threshold.R"))
    run_log$threshold = list(status = "success", time = Sys.time() - t1)
    message("\n  ✓ 03_threshold.R complete (", round(difftime(Sys.time(), t1, units = "secs"), 1), " sec)\n")
  }, error = function(e) {
    run_log$threshold <<- list(status = "failed", error = conditionMessage(e))
    message("\n  ✗ 03_threshold.R FAILED: ", conditionMessage(e), "\n")
  })
  
} else {
  message("\n  ○ 03_threshold.R skipped\n")
}

## 04: meta-analysis -----------------------------------------------------------

if (RUN_META) {
  
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  message("  04_meta.R")
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  
  t1 = Sys.time()
  
  tryCatch({
    source(here::here("code_pooled", "04_meta.R"))
    run_log$meta = list(status = "success", time = Sys.time() - t1)
    message("\n  ✓ 04_meta.R complete (", round(difftime(Sys.time(), t1, units = "secs"), 1), " sec)\n")
  }, error = function(e) {
    run_log$meta <<- list(status = "failed", error = conditionMessage(e))
    message("\n  ✗ 04_meta.R FAILED: ", conditionMessage(e), "\n")
  })
  
} else {
  message("\n  ○ 04_meta.R skipped\n")
}

## 05: subgroups ---------------------------------------------------------------

if (RUN_SUBGROUPS) {
  
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  message("  05_subgroups.R")
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  
  t1 = Sys.time()
  
  tryCatch({
    source(here::here("code_pooled", "05_subgroups.R"))
    run_log$subgroups = list(status = "success", time = Sys.time() - t1)
    message("\n  ✓ 05_subgroups.R complete (", round(difftime(Sys.time(), t1, units = "secs"), 1), " sec)\n")
  }, error = function(e) {
    run_log$subgroups <<- list(status = "failed", error = conditionMessage(e))
    message("\n  ✗ 05_subgroups.R FAILED: ", conditionMessage(e), "\n")
  })
  
} else {
  message("\n  ○ 05_subgroups.R skipped\n")
}

## 06: figures -----------------------------------------------------------------

if (RUN_FIGURES) {
  
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  message("  06_figures.R")
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  
  t1 = Sys.time()
  
  tryCatch({
    source(here::here("code_pooled", "06_figures.R"))
    run_log$figures = list(status = "success", time = Sys.time() - t1)
    message("\n  ✓ 06_figures.R complete (", round(difftime(Sys.time(), t1, units = "secs"), 1), " sec)\n")
  }, error = function(e) {
    run_log$figures <<- list(status = "failed", error = conditionMessage(e))
    message("\n  ✗ 06_figures.R FAILED: ", conditionMessage(e), "\n")
  })
  
} else {
  message("\n  ○ 06_figures.R skipped\n")
}

# summary ----------------------------------------------------------------------

end_time   = Sys.time()
total_time = difftime(end_time, start_time, units = "mins")

message("\n")
message("╔══════════════════════════════════════════════════════════════════════╗")
message("║                           Run Summary                                ║")
message("╚══════════════════════════════════════════════════════════════════════╝")
message("")

# status table
status_df = data.frame(
  Script = c("01_tables", "02_discrimination", "03_threshold", 
             "04_meta", "05_subgroups", "06_figures"),
  Status = sapply(c("tables", "discrimination", "threshold", "meta", "subgroups", "figures"), function(x) {
    if (is.null(run_log[[x]])) return("skipped")
    run_log[[x]]$status
  })
)

n_success = sum(status_df$Status == "success")
n_failed  = sum(status_df$Status == "failed")
n_skipped = sum(status_df$Status == "skipped")

for (i in seq_len(nrow(status_df))) {
  icon = switch(status_df$Status[i],
                success = "✓",
                failed  = "✗",
                skipped = "○")
  message("  ", icon, " ", status_df$Script[i], ": ", status_df$Status[i])
}

message("")
message("  Total time: ", round(total_time, 1), " minutes")
message("  Completed:  ", format(end_time, "%Y-%m-%d %H:%M:%S"))
message("")

if (n_failed > 0) {
  message("  ⚠ ", n_failed, " script(s) failed - check messages above for details")
  message("")
}

# list outputs -----------------------------------------------------------------

message("== Output Files ==\n")

if (dir.exists(here::here("output", "figures"))) {
  fig_files = list.files(here::here("output", "figures"), pattern = "\\.(pdf|png)$")
  message("  Figures (", length(fig_files), " files):")
  message("    ", here::here("output", "figures"))
}

if (dir.exists(here::here("output", "tables"))) {
  tbl_files = list.files(here::here("output", "tables"), pattern = "\\.(docx|csv)$")
  if (length(tbl_files) > 0) {
    message("  Tables (", length(tbl_files), " files):")
    message("    ", here::here("output", "tables"))
  }
}

message("")
message("╔══════════════════════════════════════════════════════════════════════╗")
message("║                              Done                                    ║")
message("╚══════════════════════════════════════════════════════════════════════╝")
message("")
