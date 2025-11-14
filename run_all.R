
# run_all.R

message("Starting data import/validation...")
source("00_setup.R")

message("Cleaning data...")
source("01_cohort.R")

message("Calculating scores...")
source("02_scores.R")

message("Performing analysis...")
source("03_analysis.R")

message("Pipeline completed")
