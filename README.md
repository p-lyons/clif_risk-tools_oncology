# Early Warning Score Validation in Oncology Inpatients

> A federated analysis pipeline for validating common early warning scores (SIRS, qSOFA, MEWS, NEWS, MEWS-SF) in hospitalized adults with cancer using the CLIF consortium data model.

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
  - [03_analysis.R — Statistical Analysis](#03_analysisr--statistical-analysis)
- [Output Artifacts](#output-artifacts)
- [Privacy and Small-Cell Handling](#privacy-and-small-cell-handling)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)
- [Citation](#citation)

---

## Overview

This repository implements a federated analysis pipeline for a multi-site validation study of early warning scores in oncology inpatients. The pipeline:

1. **Builds a cohort** of adult oncology inpatients from CLIF-formatted EHR data
2. **Derives time-varying scores** (SIRS, qSOFA, MEWS, NEWS, MEWS-SF) from vitals and labs
3. **Defines clinical outcomes** (ICU transfer, in-hospital death, hospice discharge)
4. **Exports analysis-ready artifacts** with built-in QC and privacy safeguards

Each participating site runs this pipeline locally on their data. Only aggregate summary statistics are shared for pooled analysis—no patient-level data leaves the site.

### Study Design

| Element | Description |
|---------|-------------|
| **Population** | Adults (≥18 years) admitted to inpatient wards with an oncologic diagnosis |
| **Exposure** | Time-varying early warning scores during ward stay |
| **Outcome** | Composite of ICU transfer, in-hospital death, or hospice discharge |
| **Study Window** | January 1, 2016 – December 31, 2024 |

---

## Features

- **Federated architecture** — Patient-level data never leaves the site
- **Automated validation** — Schema checks, value domain verification, and sanity bounds
- **Flexible input formats** — Supports Parquet, CSV, and FST file formats
- **Robust score derivation** — LOCF imputation, FiO₂ inference, ward-restricted calculations
- **Comprehensive analyses** — Discrimination (AUROC), calibration, threshold performance, sensitivity analyses
- **Privacy safeguards** — Small-cell suppression and collapsing (n < 5)
- **Reproducible** — Fixed random seeds, version tracking, standardized output naming

---

## Repository Structure

```
├── config/
│   ├── config_clif_oncrisk.yaml    # Site-specific configuration
│   └── clif_sites.csv              # Valid site identifiers
│
├── code/
│   ├── 00_setup.R                  # Environment, config loading, CLIF validation
│   ├── 01_cohort.R                 # Cohort construction and outcome derivation
│   ├── 02_scores.R                 # Time-varying score calculation
│   ├── 03_analysis.R               # Main, sensitivity, and subgroup analyses
│   └── run_all.R                   # Orchestration script
│
├── proj_tables/                    # Intermediate tables (auto-created)
│   ├── cohort.parquet
│   ├── scores_full.parquet
│   ├── careprocess.parquet
│   ├── outcome_times.csv
│   └── hid_jid_crosswalk.parquet
│
├── proj_output/                    # Final outputs for pooling (auto-created)
│   ├── main/                       # Primary analysis artifacts
│   ├── sensitivity/                # Sensitivity analysis variants
│   ├── horizon/                    # 12h and 24h prediction windows
│   ├── threshold/                  # Threshold-based metrics
│   ├── meta/                       # Coefficients for meta-analysis
│   └── diagnostics/                # QC and validation outputs
│
├── upload_to_box/                  # Files to share with coordinating center
│
├── oncrisk.Rproj                   # RStudio project file
└── README.md
```

---

## Prerequisites

### Required CLIF Tables

The pipeline requires the following tables in your `clif_data_location` directory, named as `clif_<table>.<file_type>`:

| Table | Purpose |
|-------|---------|
| `clif_patient` | Demographics, death dates |
| `clif_hospitalization` | Admission/discharge times and dispositions |
| `clif_adt` | ADT events for ward/ICU identification |
| `clif_hospital_diagnosis` | ICD-10-CM codes for cancer identification |
| `clif_vitals` | Heart rate, respiratory rate, temperature, blood pressure, SpO₂ |
| `clif_labs` | WBC, lactate (for SIRS criteria) |
| `clif_respiratory_support` | Device type, FiO₂, flow rates |
| `clif_medication_admin_continuous` | Vasopressor administration |
| `clif_code_status` | Code status for sensitivity analysis |
| `clif_patient_assessments` | GCS scores (for qSOFA, MEWS, NEWS) |

### System Requirements

- **R** ≥ 4.1.0
- **RAM** ≥ 16 GB recommended (scales with dataset size)
- **Disk** — Sufficient space for intermediate Parquet files

### R Package Dependencies

Packages are automatically installed on first run if not present:

```r
# Data manipulation
data.table, tidytable, collapse, tidyverse, janitor

# File I/O
arrow, fst, readr, readxl, yaml

# Analysis
pROC, glmmTMB, comorbidity

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

2. **Open the RStudio project**

   Double-click `oncrisk.Rproj` or open it from RStudio.

3. **Configure your site settings**

   Copy and edit the configuration file:

   ```bash
   cp config/config_clif_oncrisk_template.yaml config/config_clif_oncrisk.yaml
   ```

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

# Output directory (proj_tables/ and proj_output/ created here)
project_location: "/path/to/project/root"
```

### Configuration Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `site_lowercase` | Your site's identifier (case-insensitive) | `"emory"`, `"hopkins"` |
| `file_type` | Format of CLIF input files | `"parquet"`, `"csv"`, `"fst"` |
| `clif_data_location` | Directory containing `clif_*.{file_type}` files | `"/data/clif/extract_2024"` |
| `project_location` | Root directory for outputs | `"/projects/oncrisk"` |

### Validating Your Site

Ensure your site identifier appears in `config/clif_sites.csv`:

```csv
site_name
emory
hopkins
ohsu
rush
ucmc
umn
upenn
```

---

## Usage

### Quick Start

```r
# From within RStudio with the project open:
source("code/run_all.R")
```

This executes all pipeline stages in sequence.

### Running Individual Scripts

For debugging or partial runs:

```r
# Stage 1: Environment setup and CLIF validation
source("code/00_setup.R")

# Stage 2: Cohort construction and outcome derivation
source("code/01_cohort.R")

# Stage 3: Time-varying score calculation
source("code/02_scores.R")

# Stage 4: Statistical analyses and artifact export
source("code/03_analysis.R")
```

### Expected Runtime

| Stage | Typical Runtime | Notes |
|-------|-----------------|-------|
| 00_setup | 1–5 min | Depends on file format and I/O speed |
| 01_cohort | 5–15 min | Scales with hospitalization count |
| 02_scores | 10–30 min | Scales with vitals/labs volume |
| 03_analysis | 15–45 min | Bootstrap iterations are parallelized |

---

## Pipeline Details

### 00_setup.R — Environment and Validation

**Purpose:** Initialize the computing environment, load configuration, and validate CLIF data.

**Key Operations:**

1. **Resource detection** — Identifies available CPU cores and RAM; sets conservative thread limits for Arrow, data.table, and collapse
2. **Configuration loading** — Reads `config_clif_oncrisk.yaml` and validates against `clif_sites.csv`
3. **CLIF table loading** — Reads all required tables (Parquet/CSV/FST)
4. **Schema validation** — Verifies required columns exist (case-insensitive matching)
5. **Value domain checks** — Validates categorical values against expected domains

**Outputs:** Loaded CLIF tables in memory; validation report to console.

---

### 01_cohort.R — Cohort Construction

**Purpose:** Build the analytic cohort with inclusion/exclusion criteria and derive outcomes.

#### Inclusion Criteria

- Age ≥ 18 years at admission
- At least one inpatient ward admission during study window
- At least one complete vital sign set during ward stay
- Minimum 6 hours of ward time before any outcome

#### Exclusion Criteria

- Labor & delivery, psychiatry, or rehabilitation units
- ICU admission before any ward time
- Outcomes occurring within 6 hours of ward admission
- Post-death admissions (data quality exclusion)

#### Cancer Identification

Cancer status is determined from ICD-10-CM diagnosis codes:

| Category | ICD-10-CM Codes | Priority |
|----------|-----------------|----------|
| Metastatic | C77–C79, C80 | Highest |
| Hematologic | C81–C96 | High |
| High-mortality solid | C22, C25, C34 (hepatobiliary, pancreatic, lung) | Medium |
| Other solid tumors | C00–C76 (excluding C44 skin) | Lower |

Exclusions:
- Non-melanoma skin cancer (C44)
- Diagnoses documented as "in remission"

A hierarchy assigns each encounter a single cancer diagnosis when multiple are present.

#### Outcome Definition

The primary composite outcome is the **first occurrence** of:

1. **Ward → ICU transfer** (excluding scheduled OR/procedure transfers)
2. **In-hospital death**
3. **Hospice discharge**

A secondary outcome excludes hospice discharge.

#### Comorbidity Scoring

- **Elixhauser comorbidities** via the Quan ICD-10 algorithm
- **van Walraven weighted index** for summary scoring

**Outputs:**

| File | Description |
|------|-------------|
| `proj_tables/cohort.parquet` | Final analytic cohort |
| `proj_tables/outcome_times.csv` | Outcome timestamps |
| `proj_tables/careprocess.parquet` | Vasopressor and respiratory support |
| `proj_tables/hid_jid_crosswalk.parquet` | Hospitalization ID mapping |
| `proj_output/table_02_cont_<site>.csv` | Continuous variable summaries |
| `proj_output/table_02_cat_<site>.csv` | Categorical variable summaries |

---

### 02_scores.R — Score Derivation

**Purpose:** Calculate time-varying early warning scores from vitals and laboratory values.

#### Scores Implemented

| Score | Components | Threshold |
|-------|------------|-----------|
| **SIRS** | Temperature, heart rate, respiratory rate, WBC | ≥ 2 |
| **qSOFA** | Respiratory rate, systolic BP, GCS | ≥ 2 |
| **MEWS** | Heart rate, respiratory rate, systolic BP, temperature, AVPU | ≥ 5 |
| **NEWS** | Respiratory rate, SpO₂, supplemental O₂, temperature, systolic BP, heart rate, consciousness | ≥ 5 |
| **MEWS-SF** | MEWS components + SpO₂/FiO₂ ratio | ≥ 7 |

#### Imputation Rules

| Parameter | Rule |
|-----------|------|
| **Vitals** | Last observation carried forward (LOCF) up to 4 hours |
| **Labs** | LOCF up to 12 hours |
| **FiO₂** | Imputed from oxygen flow rate when not directly documented |
| **SpO₂/FiO₂ ratio** | Calculated only when SpO₂ < 97% (to avoid ceiling effects) |

#### FiO₂ Imputation from Flow Rate

```
If flow ≤ 6 LPM:  FiO₂ = 0.21 + (0.04 × flow)
If flow > 6 LPM:  FiO₂ = 0.21 + (0.04 × 6) + (0.02 × (flow - 6))
Room air:         FiO₂ = 0.21
```

FiO₂ values are carried forward up to 6 hours.

#### Score Calculation Window

- Scores are only calculated during **ward intervals**
- Observation end time is the minimum of: outcome time, ward end time, or discharge time
- Each observation row represents a point-in-time score

**Outputs:**

| File | Description |
|------|-------------|
| `proj_tables/scores_full.parquet` | Complete time-varying score dataset |

---

### 03_analysis.R — Statistical Analysis

**Purpose:** Generate all analysis artifacts for pooled meta-analysis.

#### Analysis Types

**1. Discrimination (AUROC)**
- Encounter-level maximum score vs. outcome
- Site-level AUROCs with DeLong confidence intervals
- Stratified by cancer status

**2. Threshold Performance**
- Sensitivity, specificity, PPV, NPV at standard thresholds
- Time to first positivity
- Cumulative incidence of score positivity

**3. Prediction Horizons**
- 12-hour and 24-hour outcome prediction windows
- Score value at each time point linked to future outcome

**4. Sensitivity Analyses**

| Variant | Description |
|---------|-------------|
| `main` | Primary analysis |
| `se_fullcode_only` | Restricted to full code status encounters |
| `se_no_ed_req` | No emergency department requirement |
| `se_win0_96h` | Expanded 96-hour prediction window |
| `se_one_enc_per_pt` | One randomly selected encounter per patient |

**5. Subgroup Analyses**
- Hematologic vs. solid malignancies
- Site-level heterogeneity

**6. Meta-Analysis Inputs**
- Per-score logistic regression coefficients
- Cancer × score interaction terms
- Random intercepts by hospital (if applicable)

**7. Bootstrap Confidence Intervals**
- Clustered bootstrap (one observation per encounter)
- 100 iterations with fixed seed

---

## Output Artifacts

All outputs follow a standardized naming convention:

```
proj_output/<analysis>/<artifact>[_<strata>][_h<hours>][-<variant>]-<site>.csv
```

### Directory Structure

```
proj_output/
├── main/
│   ├── auroc-ca-<site>.csv           # Site-level AUROCs by cancer status
│   ├── maxscores-ca-<site>.csv       # Encounter maximum scores
│   └── maxscores-liquid-<site>.csv   # Heme/solid subgroup
│
├── horizon/
│   ├── counts-ca-h24-<site>.csv      # 24h prediction counts
│   ├── counts-ca-h12-<site>.csv      # 12h prediction counts
│   ├── counts-ca-h24-boot-<site>.csv # Bootstrap samples
│   └── auroc-ca-h24-<site>.csv       # Horizon AUROCs
│
├── threshold/
│   ├── sesp-ca-<site>.csv            # Sensitivity/specificity
│   ├── ever-ca-<site>.csv            # Ever-positive analysis
│   ├── first-ca-<site>.csv           # Time to first positivity
│   ├── cuminc-ca-<site>.csv          # Cumulative incidence
│   └── upset-ca-<site>.csv           # Score co-positivity patterns
│
├── sensitivity/
│   ├── maxscores-ca-se_fullcode_only-<site>.csv
│   ├── counts-ca-h24-se_no_ed_req-<site>.csv
│   └── ...
│
├── meta/
│   ├── coefficients-<site>.csv       # Regression coefficients
│   └── score_sds-<site>.csv          # Score standard deviations
│
└── diagnostics/
    ├── overall-<site>.csv            # Cohort-level counts by variant
    ├── by_cancer-<site>.csv          # Counts by cancer status
    └── max_scores-<site>.csv         # Score distribution diagnostics
```

### Files to Upload

After successful completion, upload the entire `upload_to_box/` directory to the coordinating center. This includes:

- All files from `proj_output/`
- Flow diagram data (`figure_s01_flow_<site>.csv`)
- Table 2 summaries

---

## Privacy and Small-Cell Handling

The pipeline implements multiple safeguards to prevent identification of individuals:

### Small-Cell Suppression

Any cell with 0 < n < 5 is handled by:

1. **Suppression** — Value replaced with `NA` or redacted
2. **Collapsing** — Adjacent cells merged until n ≥ 5

### Collapsing Rules

The `collapse_small()` function:

- Merges sparse adjacent score values
- Stops if > 1% of total observations would be affected
- Logs all collapsing actions for audit

### Hard Caps

- Minimum cell size: n = 5
- Upset plot counts: any n < 5 recoded to 5

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
| "Column not found" | Schema mismatch | Check column names in CLIF tables; resolver is case-insensitive but names must exist |
| Arrow read errors | Missing system libraries | Install Arrow system dependencies or use `file_type: "csv"` |
| "Sanity check failed" | Outcome rates out of expected range | Review `cohort.parquet` for local coding issues |
| "Small-cell collapse exceeded 1%" | Too many sparse cells | Review score distributions; consider broader binning |
| Memory errors | Insufficient RAM | Reduce thread count in `00_setup.R` or process in chunks |

### Validation Failures

The pipeline performs extensive validation. Common issues:

1. **Missing columns** — Ensure all required CLIF columns are present
2. **Invalid values** — Check categorical columns against expected domains
3. **Date formats** — Timestamps should be in ISO 8601 format

### Getting Help

1. Check the console output for specific error messages
2. Review `proj_tables/` for intermediate outputs
3. Contact the coordinating center with:
   - Error message
   - `sessionInfo()` output
   - Approximate dataset size

---

## Contributing

This is a collaborative project across CLIF consortium sites. To suggest improvements:

1. Open an issue describing the proposed change
2. For code changes, submit a pull request with:
   - Clear description of the change
   - Testing performed at your site
   - Any impact on output artifacts

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
