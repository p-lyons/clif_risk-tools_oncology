
# run_all.R

message("Starting data import/validation...")
source("code/00_setup.R")

message("Cleaning data...")
source("code/01_cohort.R")

message("Calculating scores...")
source("code/02_scores.R")

message("Performing analysis...")
source("code/03_analysis.R")

message("Pipeline completed")