# Early Warning Score Validation in Oncology Inpatients

> A federated analysis pipeline for validating common early warning scores (SIRS, qSOFA, MEWS, MEWS-SF, NEWS2) for clinical deterioration in hospitalized adults with and without cancer, using the CLIF consortium data model.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Repository Structure](#repository-structure)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Configuration](#configuration)
- [Usage](#usage)
- [Pipeline Details](#pipeline-details)
  - [00_setup.R — Environment and Validation](#00_setupr--environment-and-validation)
  - [01_cohort.R — Cohort Construction](#01_cohortr--cohort-construction)
  - [02_scores.R — Score Derivation](#02_scoresr--score-derivation)
  - [02b_carryforward.R — Carry-Forward Sensitivity](#02b_carryforwardr--carry-forward-sensitivity)
  - [02c_monitoring.R — Monitoring Intensity](#02c_monitoringr--monitoring-intensity)
  - [03a_artifacts.R — Discrimination, Thresholds, Sensitivity](#03a_artifactsr--discrimination-thresholds-sensitivity)
  - [03b_leadtime.R — Lead Time to Deterioration](#03b_leadtimer--lead-time-to-deterioration)
  - [03c_models.R — Longitudinal GEE Models](#03c_modelsr--longitudinal-gee-models)
  - [run_all.R — Orchestration and Manifest](#run_allr--orchestration-and-manifest)
- [Output Artifacts](#output-artifacts)
- [Privacy and Small-Cell Handling](#privacy-and-small-cell-handling)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)
- [Citation](#citation)

---

## Overview

This repository implements a federated analysis pipeline for a multi-site validation study of early warning scores in hospitalized adults, comparing performance between patients with and without cancer. The pipeline:

1. **Builds a cohort** of adult ward inpatients from CLIF-formatted EHR data
2. **Derives time-varying scores** (SIRS, qSOFA, MEWS, MEWS-SF, NEWS2) from vitals and labs
3. **Defines clinical outcomes** over a competing-risk universe of ICU transfer, in-hospital death, and hospice discharge
4. **Exports analysis-ready aggregate artifacts** with built-in QC and privacy safeguards

Each participating site runs this pipeline locally on its own data. Only aggregate summary statistics are shared for pooled analysis — no patient-level data ever leaves the site. The design keeps site outputs poolable without disclosure: continuous summaries, for example, are exported as sums and sums-of-squares (`age_sum`, `age_sumsq`) rather than means and standard deviations, so the coordinating center recovers pooled moments exactly without any site releasing a patient.

This is the **round-two (revision)** pipeline. It supersedes the single `03_analysis.R` of round one with a generalized artifact engine (`03a`) and adds three reviewer-driven analyses: a carry-forward sensitivity analysis (`02b`), a monitoring-intensity export (`02c`), a lead-time analysis (`03b`), and a longitudinal GEE model (`03c`). Round-two artifacts are written to `upload_to_box_v2/`; the round-one artifacts remain frozen in `upload_to_box/`.

### Study Design

| Element | Description |
|---------|-------------|
| **Population** | Adults (≥18 years) admitted through the ED to inpatient wards |
| **Exposure** | Time-varying early warning scores during ward stay |
| **Primary Outcome** | Composite of ward-to-ICU transfer, ward death, or hospice discharge |
| **Secondary Outcomes** | The composite excluding hospice; and each component individually |
| **Comparison** | Cancer vs. non-cancer, with hematologic/solid and metastatic strata |
| **Study Window** | January 1, 2016 – December 31, 2024 |

---

## Features

- **Federated architecture** — Patient-level data never leaves the site; only aggregate count tables are shared
- **Automated validation** — Schema checks, value-domain verification, and case-insensitive column resolution
- **Flexible input formats** — Supports Parquet, CSV, and FST file formats
- **Robust score derivation** — LOCF imputation, FiO₂ inference, ward-restricted calculation, NEWS2 supplemental-oxygen item
- **Five competing-risk outcomes** — Composite and its decompositions, each with per-outcome event and censoring times
- **Comprehensive analyses** — Discrimination (AUROC), threshold performance, fixed-horizon prediction, lead time, carry-forward and other sensitivity analyses, longitudinal GEE
- **Privacy safeguards** — Estimability floors for AUROC and small-cell suppression for diagnosis tallies
- **Reproducible** — Fixed random seeds, a per-run upload manifest with checksums, and standardized output naming

---

## Repository Structure

```
├── config/
│   ├── config_clif_oncrisk.yaml        # Site-specific configuration (tracked; see git note below)
│   ├── clif_sites.csv                  # Valid site identifiers
│   └── icd10cm_casefinding_2023.xlsx   # Cancer case-finding code list
│
├── code/
│   ├── 00_setup.R                      # Environment, config loading, CLIF validation
│   ├── 01_cohort.R                     # Cohort construction, outcomes, Table 2, admission diagnoses
│   ├── 02_scores.R                     # Time-varying score calculation (incl. NEWS2)
│   ├── 02b_carryforward.R              # Carry-forward window sensitivity analysis
│   ├── 02c_monitoring.R               # Monitoring-intensity and missingness exports
│   ├── 03a_artifacts.R                 # Discrimination, thresholds, horizons, sensitivity variants
│   ├── 03b_leadtime.R                  # Lead time from first alert to deterioration
│   ├── 03c_models.R                    # Longitudinal (GEE) score × cancer interaction
│   └── run_all.R                       # Orchestration, run report, and upload manifest
│
├── .here                               # Project-root sentinel for the here package (see Installation)
│
├── proj_tables/                        # Intermediate tables (auto-created, not uploaded)
│   ├── cohort.parquet
│   ├── outcome_times.parquet
│   ├── scores_full.parquet
│   ├── scores_components.parquet
│   ├── news_o2_stream.parquet
│   ├── vital_lab_extract.parquet
│   ├── ward_times.parquet
│   ├── careprocess.parquet
│   └── hid_jid_crosswalk.parquet
│
├── upload_to_box/                      # Round-one outputs (frozen — do not overwrite)
├── upload_to_box_v2/                   # Round-two outputs for pooling (auto-created)
│   ├── main/                           # Discrimination and max-score artifacts
│   ├── threshold/                      # Threshold metrics, lead time, crossings
│   ├── horizon/                        # Fixed-horizon prediction counts
│   ├── sensitivity/                    # Sensitivity-analysis variants
│   ├── diagnostics/                    # Monitoring, missingness, NEWS2 O₂ resolution
│   └── meta/                           # Coefficients for meta-analysis
│
└── README.md
```

---

## Prerequisites

### Required CLIF Tables

The pipeline requires the following tables in your `clif_data_location` directory, named as `clif_<table>.<file_type>`:

| Table | Purpose |
|-------|---------|
| `clif_patient` | Demographics, death dates |
| `clif_hospitalization` | Admission/discharge times, ages, dispositions |
| `clif_adt` | ADT events for ward/ICU/ED identification |
| `clif_hospital_diagnosis` | ICD-10-CM codes for cancer identification and admission diagnoses |
| `clif_vitals` | Heart rate, respiratory rate, temperature, systolic BP, SpO₂ |
| `clif_labs` | WBC and pCO₂ (SIRS); broader lab set validated at load |
| `clif_respiratory_support` | Device category, FiO₂, flow rates (MEWS-SF and NEWS2 O₂ item) |
| `clif_medication_admin_continuous` | Vasopressor administration |
| `clif_code_status` | Admission code status (for the full-code sensitivity analysis) |
| `clif_patient_assessments` | GCS total (for qSOFA, MEWS, NEWS2) |

Setup validates that each table exists and that key categorical columns contain the expected values (e.g. `location_category` includes `ed`, `icu`, `ward`; `discharge_category` includes `Hospice` and `Expired`). Column matching is case-insensitive, but the columns must be present.

### System Requirements

- **R** ≥ 4.1.0
- **RAM** ≥ 20 GB recommended (peak memory can exceed 20 GB at large sites during score calculation and analysis; `run_all.R` reports peak heap use per stage)
- **Disk** — Sufficient space for intermediate Parquet files under `proj_tables/`

### R Package Dependencies

Packages are installed on first run if not already present:

```r
# Data manipulation
data.table, tidytable, collapse, tidyverse, janitor

# File I/O
arrow, fst, readr, readxl, yaml

# Analysis
pROC, glmmTMB, geepack, comorbidity

# Utilities
here, stringr, zoo, ps, rvest
```

> **Note:** The `arrow` package may require system libraries on Linux. See [Apache Arrow installation guide](https://arrow.apache.org/docs/r/articles/install.html) if you encounter build errors.

---

## Installation

1. **Clone the repository**

   ```bash
   git clone https://github.com/your-org/oncrisk-clif.git
   cd oncrisk-clif
   ```

   Cloning (rather than downloading a zip) is recommended, but either works — see the note on the `.here` sentinel below.

2. **Project root and the `here` package**

   This project does **not** use an `.Rproj` file. File paths resolve through the `here` package, which finds the project root by locating the `.here` sentinel file at the repository root. Always start R with your working directory at the repository root so `here()` roots correctly, and do not move scripts out of `code/`.

   > The `.here` sentinel is what makes downloads work: without it, `here()` would fall back to rooting off the `.git` directory, which is absent from a zip download. If `here()` cannot find the root, confirm `.here` is present at the top level of the repository.

3. **Configure your site settings**

   Edit `config/config_clif_oncrisk.yaml` with your site-specific settings (see [Configuration](#configuration)).

---

## Configuration

Edit `config/config_clif_oncrisk.yaml` with your site-specific settings:

```yaml
# Site identifier (must match an entry in clif_sites.csv)
site_lowercase: "your_site"

# Input file format: parquet, csv, or fst
file_type: "parquet"

# Path to directory containing CLIF tables
clif_data_location: "/path/to/clif/data"

# Output directory root (retained for backward compatibility; see note)
project_location: "/path/to/project/root"
```

### Configuration Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `site_lowercase` | Your site's identifier (case-insensitive; validated against `clif_sites.csv`) | `"emory"`, `"ohsu"` |
| `file_type` | Format of CLIF input files | `"parquet"`, `"csv"`, `"fst"` |
| `clif_data_location` | Directory containing `clif_*.{file_type}` files | `"/data/clif/extract_2024"` |
| `project_location` | Legacy output root | `"/projects/oncrisk"` |

> **Note on `project_location`:** This field is retained in the config but is no longer used to resolve paths. All directory creation and every write resolve through `here()` off the repository root, so `proj_tables/` and the upload directory are created there regardless of `project_location`.

### A note on git hygiene

`config/config_clif_oncrisk.yaml` is tracked but holds machine-specific absolute paths, so it will show as modified in `git status` permanently and by design. **Stage code explicitly — `git add code/` — and never `git add -A` or `git commit -a`**, or you will push your machine's paths to every site. If the permanent modified state is distracting, `git update-index --skip-worktree config/config_clif_oncrisk.yaml` hides it in your clone only. Revision work is on the `revision-jco` branch.

### Validating Your Site

Ensure your site identifier appears in `config/clif_sites.csv`. Setup will stop with the list of allowed sites if it does not.

---

## Usage

> **Important:** Start R with the working directory at the repository root before running any scripts, so the `here` package resolves paths from the `.here` sentinel.

### Standard Execution

Run the full pipeline in order with:

```r
source("code/run_all.R")
```

This sources every stage in sequence in the global environment, times each stage and records its peak memory, writes a per-site run report, and builds an upload manifest with a completeness check (see [run_all.R](#run_allr--orchestration-and-manifest)).

### Memory-Constrained Execution

If RAM is limited, run each script individually, allowing R to garbage-collect between stages, or restart R between stages for the lowest peak footprint:

```r
source("code/00_setup.R")   # environment + CLIF validation
source("code/01_cohort.R")  # cohort + outcomes + Table 2 + admission dx
source("code/02_scores.R")  # time-varying scores (writes scores_full + intermediates)
source("code/02b_carryforward.R")
source("code/02c_monitoring.R")
source("code/03a_artifacts.R")
source("code/03b_leadtime.R")
source("code/03c_models.R")
```

Each `0x` stage depends on the intermediates written by the previous one (in `proj_tables/`), so run them in order. `run_all.R` also builds the upload manifest at the end; if you run stages by hand, source the manifest section of `run_all.R` afterward (or run `run_all.R` once end-to-end) to produce `MANIFEST-<site>.csv`.

### Expected Runtime

Runtime scales with encounter and vitals/lab volume. For orientation, a site of roughly 145,000 encounters completes the full pipeline in about six minutes; larger sites take proportionally longer.

| Stage | Notes |
|-------|-------|
| 00_setup | Fast; dominated by file I/O and validation |
| 01_cohort | Scales with hospitalization count |
| 02_scores | Heaviest stage; scales with vitals/labs volume |
| 02b / 02c | Replays scoring at alternate windows / summarizes monitoring |
| 03a | Discrimination, thresholds, horizons, sensitivity variants |
| 03b / 03c | Lead time; GEE (encounter-sampled for tractability) |

---

## Pipeline Details

### 00_setup.R — Environment and Validation

**Purpose:** Initialize the computing environment, load configuration, and validate CLIF data.

**Key Operations:**

1. **Package management** — Installs missing packages and loads the working set
2. **Resource detection** — Identifies cores and available RAM; sets conservative thread limits for Arrow, data.table, and collapse
3. **Configuration loading** — Reads `config_clif_oncrisk.yaml` and validates `site_lowercase` against `clif_sites.csv`
4. **Artifact destination** — Sets the single `BOX_DIR` constant (`upload_to_box_v2`) that governs every write in the pipeline
5. **CLIF table loading** — Opens all required tables (Parquet/CSV/FST)
6. **Schema and value-domain validation** — Verifies required columns exist (case-insensitive) and expected categorical values are present

**Outputs:** Loaded CLIF tables in memory; `BOX_DIR` and site constants; validation report to console.

---

### 01_cohort.R — Cohort Construction

**Purpose:** Build the analytic cohort, derive outcomes, and export Table 2 and admission-diagnosis summaries.

#### Inclusion Criteria

- Age ≥ 18 years at admission
- At least one inpatient ward stay during the study window (stepdown is reclassified as ward)
- At least one complete set of the five required vital signs during the ward stay
- At least 6 hours of ward time before any outcome

#### Exclusion Criteria

- Labor & delivery, psychiatry, or rehabilitation units
- ICU stay before any ward time
- Encounters with outcomes within 6 hours of ward admission
- Post-death admissions and duplicate death records (data-quality exclusions)

Contiguous hospitalizations within 6 hours are linked into a single `joined_hosp_id`.

#### Cancer Identification

Cancer status is determined from ICD-10-CM diagnosis codes using the bundled case-finding list (`config/icd10cm_casefinding_2023.xlsx`). A priority hierarchy assigns each encounter a single cancer category:

| Category | ICD-10-CM Codes | Priority (`rank_enc`) |
|----------|-----------------|-----------------------|
| Metastatic | C77–C79, C80 | 1 (highest) |
| Hematologic | C81–C96 | 2 |
| High-risk solid | C22, C25, C34 (hepatobiliary, pancreatic, lung) | 3 |
| Other solid | C18, C50, C61 | 4 |
| Other | remaining malignant/hematopoietic codes | 5 |

Exclusions: non-melanoma skin cancer (C44), codes marked "in remission" or "personal history," Z85 personal-history codes, and a short list of non-specific hematologic codes.

Two derived flags accompany `rank_enc`: `liquid_01` (any hematologic code on the encounter, assigned independently of the hierarchy) and `mets_01`, which is deliberately three-valued (1 / 0 / missing) and restricted to solid tumors, so it reads as a solid-tumor stage contrast rather than a whole-cohort flag.

#### Outcome Definition

Outcomes are defined over the competing-risk universe of ICU transfer, in-hospital death, and hospice discharge. `outcome_times.parquet` spans **every** encounter in the final cohort (not only those with an event) and carries, for each of five outcomes, an indicator (`o_*_01`), an event time (`event_*_dttm`), and a censoring time (`censor_*_dttm`):

| Outcome key | Event set | Complement (censoring) |
|-------------|-----------|------------------------|
| `composite` | ICU, death, hospice | none |
| `nohospice` | ICU, death | hospice |
| `wardicu` | ICU | death, hospice |
| `warddeath` | death | ICU, hospice |
| `hospicedc` | hospice | ICU, death |

The ward-to-ICU transfer excludes moves that reach the ICU through a procedural location. Event and censoring times use strict inequalities, so an exact tie leaves both missing; ties are counted and reported. The composite reproduces the round-one `outcome_dttm` exactly (asserted in QC).

#### Comorbidity Scoring

Elixhauser comorbidities via the Quan ICD-10 algorithm, summarized by the van Walraven weighted index.

#### Admission Diagnoses

Primary admission diagnoses are exported at ICD-10-CM chapter and three-character stem level, restricted to ED-admitted encounters, suppressed at n ≤ 5. Stems are exported unmapped; the coordinating center maps them to CCSR centrally.

**Outputs:**

| File | Description |
|------|-------------|
| `proj_tables/cohort.parquet` | Final analytic cohort (one row per encounter) |
| `proj_tables/outcome_times.parquet` | Five-outcome event and censoring times |
| `proj_tables/careprocess.parquet` | Vasopressor and respiratory-support events |
| `proj_tables/hid_jid_crosswalk.parquet` | Hospitalization-ID mapping |
| `<BOX_DIR>/table_02_cont_<site>.csv` | Continuous variable summaries (sums/sums-of-squares) |
| `<BOX_DIR>/table_02_cat_<site>.csv` | Categorical variable summaries |
| `<BOX_DIR>/figure_s01_flow_<site>.csv` | Inclusion flow-diagram counts |
| `<BOX_DIR>/cancer_codes_primary_<site>.csv` | Primary cancer-code frequencies (n > 5) |
| `<BOX_DIR>/missing_demog_<site>.csv` | Demographic missingness |
| `<BOX_DIR>/admission_dx_chapter-ca-<site>.csv` | Admission diagnosis by chapter (n > 5) |
| `<BOX_DIR>/admission_dx_stem-ca-<site>.csv` | Admission diagnosis by 3-char stem (n > 5) |

---

### 02_scores.R — Score Derivation

**Purpose:** Calculate time-varying early warning scores from vitals and laboratory values, on the wards.

#### Scores Implemented

| Score | Components | Standard Threshold |
|-------|------------|--------------------|
| **SIRS** | Temperature, heart rate, respiratory rate, WBC, pCO₂ | ≥ 2 |
| **qSOFA** | Respiratory rate, systolic BP, GCS | ≥ 2 |
| **MEWS** | Heart rate, respiratory rate, systolic BP, temperature, GCS | ≥ 5 |
| **MEWS-SF** | MEWS components + SpO₂/FiO₂ ratio | ≥ 7 |
| **NEWS2** | Respiratory rate, SpO₂, supplemental O₂, temperature, systolic BP, heart rate, GCS | ≥ 5 (or any single parameter = 3) |

> **NEWS2, not NEWS.** Round two upgrades the round-one NEWS to NEWS2 by adding the 2-point supplemental-oxygen item (`news_o2`): 2 points when the patient is on any supplemental oxygen (FiO₂ > 0.21 or flow > 0), 0 on room air. The column is still named `news_total`; it now equals the round-one NEWS total plus the O₂ item. The escalation threshold (≥5) and the single-parameter "red score" rule (any component = 3) are unchanged between NEWS and NEWS2 — the O₂ item scores only 0 or 2 and cannot trigger the single-parameter rule. The O₂ item is carried onto the existing score grid by a 6-hour rolling as-of join, which leaves SIRS, qSOFA, MEWS, and MEWS-SF rows unchanged from round one. A per-site diagnostic (`news_o2_resolution`) records how each O₂ value was resolved.

#### Imputation Rules

| Parameter | Rule |
|-----------|------|
| **Vitals** | Last observation carried forward (LOCF) up to 6 hours |
| **Labs** | LOCF up to 12 hours |
| **FiO₂** | Imputed from oxygen flow rate when not directly documented; carried forward up to 6 hours |
| **SpO₂/FiO₂ ratio** | Calculated only when SpO₂ < 97% (to avoid a ceiling effect) |

#### FiO₂ Imputation from Flow Rate

```
If flow = 0 LPM:  FiO₂ = 0.21
If flow ≤ 6 LPM:  FiO₂ = 0.21 + (0.04 × flow)
If flow > 6 LPM:  FiO₂ = 0.21 + (0.04 × 6) + (0.03 × (flow − 6))
Room air:         FiO₂ = 0.21
```

#### Score Calculation Window

- Scores are calculated only during **ward intervals**
- Each observation row is a point-in-time score; the observation series is truncated at the composite outcome (the last score time the composite admits)
- Five per-outcome end times (`end_*_dttm`) are carried through so downstream scripts can truncate each outcome at its own horizon

**Outputs:**

| File | Description |
|------|-------------|
| `proj_tables/scores_full.parquet` | Complete time-varying score dataset with outcomes and strata |
| `proj_tables/scores_components.parquet` | Point-assigned components before carry-forward (input to 02b) |
| `proj_tables/news_o2_stream.parquet` | Supplemental-oxygen measurement stream before the 6h carry (input to 02b) |
| `proj_tables/vital_lab_extract.parquet` | One row per measurement, whole encounter (input to 02c) |
| `proj_tables/ward_times.parquet` | Ward-stay intervals (input to 02b/02c) |
| `<BOX_DIR>/diagnostics/news_o2_resolution-<site>.csv` | NEWS2 O₂-item resolution counts by cancer status |

---

### 02b_carryforward.R — Carry-Forward Sensitivity

**Purpose:** Test the sensitivity of results to the LOCF window. Replays the post-extraction scoring pipeline at 2-, 6-, and 12-hour carry-forward windows from `scores_components.parquet` and `news_o2_stream.parquet`, then re-applies the ward restriction. The 6-hour pass must reproduce the main `02_scores.R` output exactly (asserted in QC). Emits maximum-score and horizon-count artifacts under the `cf6` (and companion) variant tokens in `sensitivity/`.

---

### 02c_monitoring.R — Monitoring Intensity

**Purpose:** Characterize how intensively patients are monitored, and vital/lab missingness, stratified by cancer status. From `vital_lab_extract.parquet` and the ward hours implied by `scores_full`, it counts measurements per 24 ward-hours (`monitoring-ca`), bins those rates (`monitoring_bins-ca`), and reports vital/lab missingness (`missing_vlab-ca`). This replaces the round-one non-stratified `missing_vlab` export. Outputs land in `diagnostics/`.

---

### 03a_artifacts.R — Discrimination, Thresholds, Sensitivity

**Purpose:** The generalized artifact engine. Generates the primary discrimination, threshold, horizon, and sensitivity artifacts for pooled meta-analysis, across strata and outcomes, through a single `write_artifact()` path gated on an allow-list.

**Analysis families:**

1. **Discrimination (AUROC)** — Encounter-level maximum score vs. outcome, by cancer status and finer strata, with an estimability floor (see [Privacy](#privacy-and-small-cell-handling)).
1a. **First-ward score distribution** (`firstscore`) — each score's value at the first ward observation, as a count distribution by cancer status. Descriptive (no discrimination); poolable to medians/IQRs/means. Answers Reviewer 1's presenting-acuity comment.
2. **Threshold performance** — Sensitivity, specificity, PPV/NPV at standard thresholds (`sesp`); ever-positive (`ever`); time to first positivity (`first`); cumulative incidence of positivity (`cuminc`); and score co-positivity patterns (`upset`).
3. **Fixed-horizon prediction** — Score value at each observation linked to deterioration within the following 24 hours (`counts … h24`), with a clustered bootstrap (one observation per encounter, fixed seed).
4. **Sensitivity analyses** — see table below.
5. **Meta-analysis inputs** — Per-score logistic-regression coefficients and the score × cancer interaction (glmmTMB), plus score standard deviations, in `meta/`.

| Variant | Description |
|---------|-------------|
| `main` | Primary analysis (ED admissions only) |
| `se_fullcode_only` | Restricted to full-code status (full and presume-full) |
| `se_no_ed_req` | No ED-admission requirement |
| `se_win0_96h` | Restricted to the 0–96 hour window |
| `se_one_enc_per_pt` | One randomly selected encounter per patient |
| `cf6` (from 02b) | 6-hour carry-forward replay |

**Strata** span cancer status (`ca`), hematologic vs. solid (`liquid`), metastatic (`mets`), and cancer diagnosis group (`dxgroup`). **Outcomes** span the five keys defined in `01_cohort.R`.

---

### 03b_leadtime.R — Lead Time to Deterioration

**Purpose:** Report how much warning an alert affords. For each encounter, score, and outcome, computes the hours from the first threshold-positive score to the event, exported as binned counts plus sufficient statistics (never individual intervals), together with a crossed × event classification of every encounter (`crossclass`). Positivity uses the standard thresholds with the NEWS2 single-parameter rule. Because the score series is truncated before the event, every crossing precedes it and lead times are strictly positive (asserted). A QC cross-check ties the crossed-and-event counts back to `03a`'s `sesp` artifact. Outcomes: composite and nohospice. Outputs in `threshold/`.

---

### 03c_models.R — Longitudinal GEE Models

**Purpose:** Provide the reviewer-requested longitudinal method that accounts for repeated within-encounter measurements. For each score, fits a GEE with an exchangeable working correlation clustered on the encounter and exports the score × cancer interaction with robust standard errors, for central meta-analysis.

- **Model:** `outcome_24h ~ value * ca_01 (+ hospital_id where >1 hospital)`, `family = binomial`, `id = joined_hosp_id`, `corstr = "exchangeable"`
- **Outcome:** composite at the 24-hour horizon
- **Sampling:** encounters (not observations) are sampled to a cap of 100,000 with a fixed seed, keeping all observations of the chosen encounters
- **QC:** convergence is reported (not fatal); the interaction sign is cross-checked against `03a`'s glmmTMB coefficients

Output: `meta/gee_coefficients-<site>.csv`.

---

### run_all.R — Orchestration and Manifest

**Purpose:** Run every stage in order, then produce a run report and an upload manifest.

- Sources each stage in the global environment, timing it and recording peak R heap use
- Writes `run_report_<site>.csv` (encounter counts, cancer fraction, outcome rate, runtime)
- Walks `upload_to_box_v2/` and writes `MANIFEST-<site>.csv` — every uploaded CSV with its row/column counts and MD5 checksum, plus site, run date, and pipeline version — so the coordinating center can confirm a complete, uncorrupted upload
- Runs a completeness check that fails loudly if any expected artifact family is absent (a silent stage failure)

---

## Output Artifacts

All outputs follow a standardized naming convention:

```
upload_to_box_v2/<family>/<artifact>[-<strata>][-<outcome>][-h<hours>][-<variant>]-<site>.csv
```

Tokens are optional and appear only where they apply: `strata` (e.g. `ca`, `liquid`, `mets`, `dxgroup`), `outcome` (`composite`, `nohospice`, `wardicu`, `warddeath`, `hospicedc`), `hours` (e.g. `h24`), and `variant` (e.g. `se_fullcode_only`, `cf6`).

### Directory Structure

```
upload_to_box_v2/
├── main/
│   ├── maxscores-<strata>-<outcome>-<site>.csv     # Encounter maximum scores
│   ├── auroc-<strata>-<outcome>-<site>.csv         # Site-level AUROCs (estimability-gated)
│   ├── events-<strata>-<outcome>-<site>.csv        # Event tallies
│   └── firstscore-ca-<site>.csv                    # Score value at first ward obs, by cancer status
│
├── threshold/
│   ├── ever-<strata>-<outcome>-<site>.csv          # Ever-positive analysis
│   ├── sesp-<strata>-<outcome>-<site>.csv          # Sensitivity / specificity / PPV / NPV
│   ├── first-<strata>-<outcome>-<site>.csv         # Time to first positivity
│   ├── cuminc-<strata>-<outcome>-<site>.csv        # Cumulative incidence of positivity
│   ├── upset-<strata>-<outcome>-<site>.csv         # Score co-positivity patterns
│   ├── leadtime-ca-<outcome>-<site>.csv            # Lead time (binned + sufficient stats)
│   └── crossclass-ca-<outcome>-<site>.csv          # Crossed × event classification
│
├── horizon/
│   └── counts-<strata>-<outcome>-h24-<site>.csv    # 24h prediction counts (+ bootstrap)
│
├── sensitivity/
│   ├── maxscores-<strata>-<outcome>-<variant>-<site>.csv
│   ├── counts-<strata>-<outcome>-h24-<variant>-<site>.csv
│   └── ...                                          # se_* and cf6 variants
│
├── diagnostics/
│   ├── news_o2_resolution-<site>.csv               # NEWS2 O₂-item resolution
│   ├── monitoring-ca-<site>.csv                    # Measurements per 24 ward-hours
│   ├── monitoring_bins-ca-<site>.csv               # Binned monitoring rates
│   └── missing_vlab-ca-<site>.csv                  # Vital/lab missingness by cancer status
│
├── meta/
│   ├── coefficients-<site>.csv                     # glmmTMB regression coefficients
│   ├── score_sds-<site>.csv                        # Score standard deviations
│   └── gee_coefficients-<site>.csv                 # GEE score × cancer interaction
│
├── table_02_cont_<site>.csv                        # Continuous variable summaries
├── table_02_cat_<site>.csv                         # Categorical variable summaries
├── figure_s01_flow_<site>.csv                      # Flow-diagram data
├── missing_demog_<site>.csv                        # Demographic missingness
├── cancer_codes_primary_<site>.csv                 # Cancer-code frequencies
├── admission_dx_chapter-ca-<site>.csv              # Admission diagnosis by chapter
├── admission_dx_stem-ca-<site>.csv                 # Admission diagnosis by stem
├── run_report_<site>.csv                           # Pipeline execution summary
└── MANIFEST-<site>.csv                             # Upload manifest with checksums
```

### Files to Upload

After a successful run, upload the entire `upload_to_box_v2/` directory (including `MANIFEST-<site>.csv`) to the coordinating center. Do **not** overwrite `upload_to_box/`, which holds the frozen round-one submission.

---

## Privacy and Small-Cell Handling

The pipeline applies disclosure protections that differ by artifact class, because the disclosure risk differs. These three rules are distinct and are easy to conflate:

| Artifact class | Rule |
|----------------|------|
| Score values and threshold counts | No suppression |
| AUROC | Estimability floor: computed only when there are at least 10 events **and** at least 10 non-events; otherwise returned as `NA` with the counts and an `estimable` flag retained |
| Diagnosis-code tallies (cancer codes, admission diagnoses) | Cells with n ≤ 5 are dropped |

Small-cell handling also applies to the score co-positivity (`upset`) counts.

### What Is NOT Shared

- Patient-level data
- Dates of service
- Identifiers of any kind
- Raw CLIF tables

---

## Troubleshooting

### Common Issues

| Problem | Likely Cause | Solution |
|---------|--------------|----------|
| "Column not found" | Schema mismatch | Check column names in CLIF tables; the resolver is case-insensitive but names must exist |
| Arrow read errors | Missing system libraries | Install Arrow system dependencies or use `file_type: "csv"` |
| "Sanity check failed" | Outcome rates out of expected range | Review `cohort.parquet` for local coding issues |
| `here()` cannot find root | Working directory not at repo root, or `.here` missing | Start R at the repository root; confirm the `.here` file is present |
| "BOX_DIR not found" | A stage was run before `00_setup.R` | Run stages in order, or use `run_all.R` |
| Missing intermediate parquet | A `02*`/`03*` stage run out of order | Run the earlier stage first; each stage reads the previous stage's `proj_tables/` outputs |
| Memory errors | Insufficient RAM | Run stages individually with `gc()` between, or restart R between stages |
| MANIFEST QC failure | An expected artifact family is absent | A stage failed quietly; check console output for the failing stage |

### Validation Failures

Setup performs extensive validation. Common issues:

1. **Missing columns** — Ensure all required CLIF columns are present
2. **Invalid values** — Check categorical columns against expected domains (e.g. `location_category`, `discharge_category`, `device_category`)
3. **Date formats** — Timestamps should parse as datetimes

### Getting Help

1. Check the console output for the specific stage and error message
2. Review `proj_tables/` for intermediate outputs
3. Contact the coordinating center with the error message, `sessionInfo()` output, and approximate dataset size

---

## Contributing

This is a collaborative project across CLIF consortium sites. To suggest improvements:

1. Open an issue describing the proposed change
2. For code changes, submit a pull request with a clear description, testing performed at your site, and any impact on output artifacts

Stage code explicitly (`git add code/`); never `git add -A`, to avoid pushing machine-specific config paths.

---

## License

[Add appropriate license information]

---

## Citation

If you use this pipeline, please cite:

> [Citation information to be added upon publication]

---

## Acknowledgments

This work is supported by the CLIF consortium. We thank all participating sites for their contributions to this federated analysis.

---

## Contact

- **Coordinating Center:** [Contact information]
- **Technical Issues:** [GitHub issues or email]
