
Purpose

Validate common risk scores (SIRS, qSOFA, MEWS, NEWS, MEWS+SF) in adult oncology inpatients using CLIF federated dataset. Pipeline builds a cohort, derives time-varying scores, defines outcomes, and exports analysis-ready artifacts with built-in QC and privacy safeguards.

Repository layout
config/
  ├─ config_clif_oncrisk.yaml   # site_lowercase, file_type, clif_data_location, project_location
  └─ clif_sites.csv             # list of valid sites; site_name must match site_lowercase
proj_tables/                    # intermediate tables (created if absent)
proj_output/                    # analysis outputs (created if absent)
  ├─ main/ horizon/ threshold/ sensitivity/ meta/
code/
  ├─ 00_setup.R                 # environment, config, load + validate CLIF tables
  ├─ 01_cohort.R                # cohort build, outcomes, care processes, table
  ├─ 02_scores.R                # time-varying score assignment and persistence rules
  └─ 03_analyses.R              # main + sensitivity analyses, bootstraps, subgroups

Data requirements

Required CLIF tables in config$clif_data_location as clif_<name>.<file_type>:

patient, hospital_diagnosis, hospitalization, adt, vitals, labs,
medication_admin_continuous, respiratory_support, code_status, patient_assessments

Supported file_type: parquet, csv, or fst.

Study window: 2016-01-01 to 2024-12-31. Change in 00_setup.R if needed.

Configuration (config/config_clif_oncrisk.yaml)

Keys used:

site_lowercase: one of clif_sites.csv$site_name (case-insensitive check)

file_type: parquet | csv | fst

clif_data_location: directory containing CLIF extracts

project_location: output root for proj_tables and proj_output

Dependencies

R packages auto-installed and loaded on first run:

install: comorbidity, data.table, tidyverse, tidytable, collapse, janitor,
glmmTMB, readxl, arrow, rvest, readr, yaml, here, pROC, fst, zoo, ps

load: data.table, tidytable, collapse, stringr, arrow, here, pROC, ps

System notes:

Arrow: ensure OS libraries are available if building from source.

Multithreading uses available cores and RAM; environment variables are set for Arrow and data.table.

How to run

Artifacts are written under proj_tables/ and upload_to_box/ (see below).

From repo root:
1. Populate config/config_clif_oncrisk.yaml.
2. Confirm clif_sites.csv contains your site_lowercase.
3. Open RProject
4. Open run_all.R

Alternatives
source("code/00_setup.R")    # environment, config, load Arrow/CSV/FST, validate
source("code/01_cohort.R")   # build cohort and outcomes, save cohort.parquet
source("code/02_scores.R")   # assign SIRS/qSOFA/MEWS/NEWS/MEWS+SF, save scores_full.parquet
source("code/03_analyses.R") # main + sensitivity + bootstraps + subgroups; write artifacts
source("code/run_all.R")     # runs all scripts in sequence

Pipeline details
00_setup: environment, loading, validation

Detects cores and RAM; conservatively sets threads.

Reads config_clif_oncrisk.yaml and clif_sites.csv.

Verifies required CLIF tables exist; supports .parquet, .csv, .fst.

Validation: per-table required columns and value domains (case-insensitive resolver). Fails fast with a compact report.

01_cohort: cohort and outcomes

Adults (≥18), study window, contiguous stays linked if gap < 6h.

Inpatient ward requirement; excludes L&D/psych/rehab.

Requires at least one complete ward vital set and ≥6h ward time before outcomes.

Excludes ICU before wards; excludes very early outcomes.

Cancer attribution: ICD-10-CM case-finding list (excludes C44, excludes “in remission”), hierarchy to prioritize metastatic/hemato/etc., encounter-level flags (ca_01, liquid_01).

Outcomes: first of ward→ICU transfer, in-hospital death, or hospice discharge; also a no-hospice composite.

Care processes: vasopressors and respiratory support (IMV, NIPPV, HFNC).

Comorbidity: Elixhauser (Quan ICD-10), van Walraven weighted score vw.

QC:

Duplicate joined encounters and “you only die once” checks

Post-death admissions removed

Outcome prevalence sanity bands; fails if implausible

Outputs:

proj_tables/outcome_times.csv

proj_tables/careprocess.parquet

proj_tables/cohort.parquet

proj_tables/hid_jid_crosswalk.parquet

proj_output/table_02_cont_<site>.csv

proj_output/table_02_cat_<site>.csv

02_scores: time-varying scoring

Vitals: SIRS, qSOFA, MEWS, NEWS components; SF ratio for MEWS+SF.

FiO₂ imputation from L/min with room-air handling; carryforward within 6h for FiO₂.

LOCF windows: vitals 4h, labs 12h; integer-safe NA handling.

Scores restricted to ward intervals; observation-level end time is min(outcome, ward end).

Outputs:

proj_tables/scores_full.parquet

03_analyses: main, sensitivity, and subgroups

Thresholds: SIRS≥2, qSOFA≥2, MEWS≥5, NEWS≥5, MEWS+SF≥7.

Ever-positive: time to first threshold met; exported counts by cancer status.

Horizon counts: 12h and 24h event labeling via make_y(), with:

Standard counts over all rows

Clustered bootstrap: one observation per encounter, 100 iterations

Encounter-level maxima: per-score max and outcome; optional small-cell collapsing.

Sensitivity variants:

main, se_no_ed_req, se_fullcode_only, se_win0_120h, se_one_enc_per_pt

Subgroup: liquid vs solid malignancy for 24h horizon and max-score analyses.

Meta coefficients: per-score logistic GLMM with cancer interaction; random intercept by hospital if >1 hospital present.

Small-cell handling:

Hard cap: any n < 5 → recoded to 5 for pooled upset counts.

collapse_small() merges sparse adjacent cells; stops if >1% of observations would be collapsed.

Outputs via write_artifact() naming:

Root: proj_output/<analysis>/

Filenames: <artifact>[_<strata>][_h<hours>][-<variant>]-<site>.csv

Examples:

proj_output/threshold/ever_ca-<site>.csv

proj_output/horizon/counts_ca_h12-<site>.csv

proj_output/sensitivity/counts_ca_h24-se_one_enc_per_pt-<site>.csv

proj_output/main/maxscores_ca-<site>.csv

proj_output/main/maxscores_liquid-<site>.csv

proj_output/horizon/counts_liquid_h24-<site>.csv

proj_output/threshold/upset_ca-<site>.csv

proj_output/meta/coefficients-<site>.csv

Key assumptions and rules

Case-insensitive schema: column matching tolerates site-specific casing.

Admissions: ward-based med-surg cohort; ICU before ward excluded.

Outcomes: OR/procedure ICU admissions not counted as deterioration.

Suppression: cells with 0<n<5 are suppressed or collapsed.

Determinism: fixed seeds for sampling; filenames include site and variant.

Performance

Threads chosen by available cores and RAM with headroom. Arrow, data.table, and collapse threads are aligned. Defaults are conservative to avoid OOM.

Troubleshooting

Validation fails: check required columns and value domains in each table. Resolver is case-insensitive, but names must exist.

Arrow errors: upgrade arrow or install system libs; fall back to CSV/FST if needed.

Sanity check stop: outcome proportions out of band. Inspect cohort.parquet and source tables for local coding issues.

Small-cell stop: collapse_small() exceeded 1% cap. Review binning or thresholds.

Reproducibility

Seeds: analysis and bootstrap use fixed seeds.

Versions: record sessionInfo() when submitting pooled artifacts.

Pathing: all paths via here(); set project_location in config for nonstandard layouts.






