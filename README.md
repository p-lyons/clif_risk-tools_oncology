# Early Warning Score Validation in Oncology Inpatients

> A federated analysis pipeline for validating common early warning scores (SIRS, qSOFA, MEWS, MEWS-SF, NEWS) for clinical deterioration in hospitalized adults with and without cancer, using the CLIF consortium data model.

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
  - [02_scores.R — Score Derivation and Descriptive Artifacts](#02_scoresr--score-derivation-and-descriptive-artifacts)
  - [02b_carryforward.R — Carry-Forward Sensitivity](#02b_carryforwardr--carry-forward-sensitivity)
  - [02c_monitoring.R — Monitoring Intensity](#02c_monitoringr--monitoring-intensity)
  - [03a_artifacts.R — Discrimination, Thresholds, Sensitivity](#03a_artifactsr--discrimination-thresholds-sensitivity)
  - [03b_leadtime.R — Lead Time to Deterioration](#03b_leadtimer--lead-time-to-deterioration)
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
2. **Derives time-varying scores** (SIRS, qSOFA, MEWS, MEWS-SF, NEWS) from vitals and labs
3. **Defines clinical outcomes** over a competing-risk universe of ICU transfer, in-hospital death, and hospice discharge
4. **Exports analysis-ready aggregate artifacts** with built-in QC and privacy safeguards

Each participating site runs this pipeline locally on its own data. Only aggregate summary statistics are shared for pooled analysis — no patient-level data ever leaves the site. The design keeps site outputs poolable without disclosure: continuous summaries, for example, are exported as sums and sums-of-squares (`age_sum`, `age_sumsq`) rather than means and standard deviations, so the coordinating center recovers pooled moments exactly without any site releasing a patient.

This is the **round-two (revision)** pipeline. It supersedes the single `03_analysis.R` of round one with a generalized artifact engine (`03a`) and adds three reviewer-driven analyses: a carry-forward sensitivity analysis (`02b`), a monitoring-intensity export (`02c`), and a lead-time analysis (`03b`). Round-two artifacts are written to `upload_to_box_v2/`; the round-one artifacts remain frozen in `upload_to_box/`.

Round two also replaces on-the-fly package installation with [uvr](https://github.com/nbafrank/uvr), which pins the R version and every package version through a committed lockfile. All eight sites therefore run identical package versions on whatever operating system they have, and nothing the pipeline installs touches the site's system R library. See [Prerequisites](#prerequisites) and [Installation](#installation).

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
- **Pinned, portable environment** — uvr installs the R version and the exact locked package set into a project-local library on macOS, Linux, and Windows
- **Automated validation** — Schema checks, value-domain verification, and case-insensitive column resolution
- **Flexible input formats** — Supports Parquet, CSV, and FST file formats
- **Score derivation with explicit imputation rules** — LOCF imputation, FiO₂ inference, ward-restricted calculation, the NEWS supplemental-oxygen item
- **Five competing-risk outcomes** — Composite and its decompositions, each with per-outcome event and censoring times
- **Comprehensive analyses** — Discrimination (AUROC), threshold performance, fixed-horizon prediction, lead time, carry-forward and other sensitivity analyses
- **Privacy safeguards** — Estimability floor for AUROC and small-cell suppression for diagnosis-code tallies
- **Reproducible** — Fixed random seeds, a committed lockfile, a per-run upload manifest with checksums, and standardized output naming

---

## Repository Structure

```
├── uvr.toml                            # R version constraint + direct dependencies (tracked)
├── uvr.lock                            # Exact resolved package set (tracked; do not hand-edit)
├── .r-version                          # Exact pinned R version, 4.6.1 (tracked)
├── .uvr/                               # Project-local R library (auto-created; git-ignored)
│
├── config/
│   ├── config_clif_oncrisk_template.yaml  # Template to copy (tracked)
│   ├── config_clif_oncrisk.yaml        # Site-specific configuration (tracked; see git note below)
│   ├── clif_sites.csv                  # Valid site identifiers and time zones
│   └── icd10cm_casefinding_2023.xlsx   # Cancer case-finding code list
│
├── code/
│   ├── 00_setup.R                      # Environment verification, config loading, CLIF validation
│   ├── 01_cohort.R                     # Cohort construction, outcomes, comorbidities, cancer codes
│   ├── 02_scores.R                     # Time-varying scores, then Table 2 and admission dx
│   ├── 02b_carryforward.R              # Carry-forward window sensitivity analysis
│   ├── 02c_monitoring.R                # Monitoring-intensity and missingness exports
│   ├── 03a_artifacts.R                 # Discrimination, thresholds, horizons, sensitivity variants
│   ├── 03b_leadtime.R                  # Lead time from first alert to deterioration
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
│   ├── hid_jid_crosswalk.parquet
│   └── site_lowercase.parquet
│
├── upload_to_box/                      # Round-one outputs (frozen — do not overwrite)
├── upload_to_box_v2/                   # Round-two outputs for pooling (auto-created)
│   ├── main/                           # Discrimination and max-score artifacts
│   ├── threshold/                      # Threshold metrics, lead time, crossings
│   ├── horizon/                        # Fixed-horizon prediction counts
│   ├── sensitivity/                    # Sensitivity-analysis variants
│   └── diagnostics/                    # Monitoring, missingness, NEWS O₂ resolution
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
| `clif_respiratory_support` | Device category, FiO₂, flow rates (MEWS-SF and the NEWS O₂ item) |
| `clif_medication_admin_continuous` | Vasopressor administration |
| `clif_code_status` | Admission code status (for the full-code sensitivity analysis) |
| `clif_patient_assessments` | GCS total (for qSOFA, MEWS, NEWS) |

Setup validates that each table exists and that key categorical columns contain the expected values (e.g. `location_category` includes `ed`, `icu`, `ward`; `discharge_category` includes `Hospice` and `Expired`). Column matching is case-insensitive, but the columns must be present.

> **`diagnosis_primary` is validated but tolerated when absent.** `00_setup.R` requires the column in `clif_hospital_diagnosis`. If your extract has the column but never populates it with `1`, `02_scores.R` writes header-only admission-diagnosis artifacts and prints a note rather than failing. Tell the coordinating center if you see that note.

### System Requirements

| Requirement | Value |
|-------------|-------|
| **R** | 4.6.1 exactly, pinned in `.r-version` and installed by uvr — you do not need a matching system R |
| **uvr** | A single ~8 MB executable; no administrator rights required |
| **RAM** | ≥ 20 GB recommended; peak use can exceed 20 GB at large sites during score calculation, and `run_all.R` reports peak heap use per stage |
| **Disk** | Space for intermediate Parquet files under `proj_tables/`, plus roughly 1 GB for the project library under `.uvr/` |

### R Package Dependencies

**You do not install any R package by hand.** Dependencies are declared in `uvr.toml` and resolved to exact versions in `uvr.lock`, both of which are committed. `uvr sync --frozen` installs precisely what the lockfile names, into a project-local library at `.uvr/library/`.

The direct dependencies are:

```
# Data manipulation
data.table, tidytable, collapse, dplyr, tibble, tidyr, janitor, rlang

# File I/O
arrow, fst, readxl, yaml, here

# Analysis
pROC, comorbidity

# Utilities
stringr, lubridate, ps
```

uvr resolves the full transitive closure of these into `uvr.lock`. Packages install as pre-built binaries from Posit Public Package Manager on macOS, Linux, and Windows, so nothing compiles in the normal case.

### R Version Pin

`uvr.toml` declares a floor of `>= 4.6.0`. `.r-version` pins the exact version, **4.6.1**, and that pin is what governs. Both files are committed.

The exact pin matters because pre-built package binaries are compiled against a specific R minor series. Two sites on 4.6.x and 4.7.x would not be running the same binaries even from an identical lockfile, which is how compiled dependencies fail at run time one layer down. The pin costs a site nothing: uvr installs R 4.6.1 into `~/.uvr/r-versions/` per project, with no administrator rights and no effect on the system R that machine already has.

> **Why this changed.** Round one called `install.packages()` at the top of `00_setup.R`, which meant each site ran whatever CRAN happened to serve on the day it ran, against whatever R that machine had. Version drift across sites is not always visible in results, and compiled packages in particular can fail at run time when built against mismatched versions. A committed lockfile removes that class of problem.
>
> **Windows only.** A package will occasionally fall back to a source install on Windows, which needs [Rtools](https://cran.r-project.org/bin/windows/Rtools/). `uvr doctor` reports whether Rtools is present. If `uvr sync --frozen` completes without error, you do not need it.

---

## Installation

### 1. Install the uvr binary (once per machine)

**macOS and Linux**

```bash
curl -fsSL https://raw.githubusercontent.com/nbafrank/uvr/main/install.sh | sh
```

**Windows**

Download `uvr-x86_64-pc-windows-msvc.zip` from the [latest release](https://github.com/nbafrank/uvr/releases/latest), extract `uvr.exe`, and add it to your `PATH`.

**From the R console, on any platform**

If you would rather stay in RStudio or Positron:

```r
# install.packages("pak")
pak::pak("nbafrank/uvr-r")
uvr::install_uvr()
```

Confirm the install with `uvr --version`. **The version must be 0.4.6 or newer.** If you already have uvr from another project, run `uvr self-update` before going further: releases in the 0.3.x series install R through a backend that their own version selector then rejects, which produces a "broken install" error that reinstalling does not clear.

### 2. Clone the repository

```bash
git clone https://github.com/p-lyons/clif_risk-tools_oncology.git
cd clif_risk-tools_oncology
```

Every command below runs from the repository root — the folder containing `uvr.toml`, `code/`, and `config/`. Cloning rather than downloading a zip is recommended, but either works; see the note on the `.here` sentinel below.

### 3. Install R and the locked package set

```bash
uvr sync --frozen
```

This reads `uvr.lock`, installs the pinned R version if your machine has no satisfying one, and installs every package at the exact version the lockfile names.

`--frozen` makes the command fail rather than silently re-resolve when `uvr.lock` and `uvr.toml` have drifted apart. That failure is informative: your checkout is inconsistent, and the fix is `git pull`, not a local re-lock.

If anything looks wrong, run `uvr doctor`. It reports your platform, whether binary packages are available for it, which R versions are installed, whether a compiler toolchain is present, and the state of the manifest and lockfile. Send its full output with any support request.

### 4. Project root and the `here` package

This project does **not** use an `.Rproj` file. File paths resolve through the `here` package, which finds the project root by locating the `.here` sentinel at the repository root. Launch the pipeline from the repository root and do not move scripts out of `code/`.

> The `.here` sentinel is what makes zip downloads work: without it, `here()` would fall back to rooting off the `.git` directory, which a zip does not contain. If `here()` cannot find the root, confirm `.here` is present at the top level.

### 5. Configure your site settings

Copy the template, then edit the copy — see [Configuration](#configuration).

```bash
cp config/config_clif_oncrisk_template.yaml config/config_clif_oncrisk.yaml
```

On Windows, use `copy`, or copy the file in Explorer.

---

## Configuration

Edit `config/config_clif_oncrisk.yaml`:

```yaml
# Path to directory containing CLIF tables
clif_data_location: "/path/to/clif/data"

# Output directory root (retained for backward compatibility; see note)
project_location: "/path/to/project/root"

# Input file format: parquet, csv, or fst
file_type: "parquet"

# Site identifier (must match an entry in clif_sites.csv)
site_lowercase: "your_site"

# Local time zone
time_zone: "America/Los_Angeles"

# Sparse-oxygen override — leave false unless instructed otherwise
allow_sparse_o2: false
```

### Configuration Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `clif_data_location` | Directory containing `clif_*.{file_type}` files | `"/data/clif/extract_2024"` |
| `project_location` | Legacy output root | `"/projects/oncrisk"` |
| `file_type` | Format of CLIF input files | `"parquet"`, `"csv"`, `"fst"` |
| `site_lowercase` | Your site's identifier; trimmed and lowercased, then validated against `clif_sites.csv` | `"emory"`, `"ohsu"` |
| `time_zone` | Your site's local time zone | `"America/Chicago"` |
| `allow_sparse_o2` | Downgrades the NEWS sparse-oxygen stop to a warning | `false` |

> **Note on `project_location`:** This field is retained in the config but is no longer used to resolve paths. All directory creation and every write resolve through `here()` off the repository root, so `proj_tables/` and the upload directory are created there regardless of `project_location`.

> **Note on `allow_sparse_o2`:** `02_scores.R` stops the run when more than 90% of NEWS score rows have no supplemental-oxygen measurement within 6 hours, which means `lpm_set` and `fio2_set` exist as columns in `clif_respiratory_support` but are effectively empty. Setting `allow_sparse_o2: true` downgrades that stop to a warning, and the run finishes with NEWS scored without its supplemental-oxygen item. **Set it to `true` only after the coordinating center has confirmed that your oxygen data are genuinely absent.** The key is optional — an older config that omits it reads as `false` — but the value must be unquoted `true` or `false`; a quoted `"false"` is rejected rather than coerced.

### A note on git hygiene

`config/config_clif_oncrisk.yaml` is tracked but holds machine-specific absolute paths, so it will show as modified in `git status` permanently and by design. **Stage code explicitly — `git add code/` — and never `git add -A` or `git commit -a`**, or you will push your machine's paths to every site. If the permanent modified state is distracting, `git update-index --skip-worktree config/config_clif_oncrisk.yaml` hides it in your clone only.

Add `.uvr/` to `.gitignore`; the project library is machine-specific and must not be committed. `uvr.toml`, `uvr.lock`, and `.r-version` must be committed.

### Validating Your Site

Ensure your site identifier appears in `config/clif_sites.csv`. Setup stops with the list of allowed sites if it does not.

---

## Usage

> **Important:** Launch through uvr, from the repository root. `00_setup.R` runs two startup guards: it confirms the session is using the project library at `.uvr/library/`, and it confirms the running R matches the version pinned in `.r-version`. Either mismatch stops the run in the first seconds rather than producing results against the wrong environment.
>
> The second guard is not redundant. uvr 0.4.x links the project library into a plain `Rscript` session, so packages resolve correctly even outside `uvr run` — but the R version in that case is whatever the shell finds, which at a site with an older system R would run 4.6.1 binaries under 4.5.x. `uvr run` is still the supported launch path; the version guard is what makes a mistake visible.

### Upgrading from an earlier round-two checkout

If this machine ever ran an earlier round-two version of the pipeline, delete the
output directory before the first run of this version:

```bash
rm -rf upload_to_box_v2
```

Earlier versions wrote artifact families this version no longer produces (a
`meta/` directory, `dxgroup` strata). The pipeline only deletes files it
regenerates, and the manifest lists **everything** in `upload_to_box_v2/`, so
stale files from an old run would otherwise ship to the coordinating center.
A clean directory costs nothing: `00_setup.R` recreates the structure.

### Standard Execution

```bash
uvr run code/run_all.R
```

Or from the R console:

```r
uvr::run("code/run_all.R")
```

This sources every stage in sequence in the global environment, times each stage and records its peak memory, writes a per-site run report, and builds an upload manifest with a completeness check (see [run_all.R](#run_allr--orchestration-and-manifest)).

### Memory-Constrained or Stage-by-Stage Execution

If RAM is limited, or you are debugging one stage, open an interactive R session with the project library active and source stages individually:

```bash
uvr run
```

Then, at the R prompt:

```r
source("code/00_setup.R")   # environment verification + CLIF validation
source("code/01_cohort.R")  # cohort + outcomes + comorbidities + cancer codes
source("code/02_scores.R")  # scores, then Table 2, missingness, admission dx
source("code/02b_carryforward.R")
source("code/02c_monitoring.R")
source("code/03a_artifacts.R")
source("code/03b_leadtime.R")
```

`uvr run` with no script argument starts R with `.uvr/library/` on the search path, which is what satisfies the guard in `00_setup.R`. Each `0x` stage depends on the intermediates the previous one wrote to `proj_tables/`, so run them in order. `run_all.R` also builds the upload manifest at the end; if you run stages by hand, source the manifest section of `run_all.R` afterward, or run `run_all.R` once end to end, to produce `MANIFEST-<site>.csv`.

### Expected Runtime

Runtime scales with encounter and vitals/lab volume. For orientation, a site of roughly 145,000 encounters completes the full pipeline in about 30 minutes (the bootstrap and the fixed-horizon families dominate); larger sites take proportionally longer. The first `uvr sync --frozen` adds a few minutes once, and subsequent syncs are near-instant.

| Stage | Notes |
|-------|-------|
| 00_setup | Fast; dominated by file I/O and validation |
| 01_cohort | Scales with hospitalization count |
| 02_scores | Heaviest stage; scales with vitals/labs volume |
| 02b / 02c | Replays scoring at alternate windows / summarizes monitoring |
| 03a | Discrimination, thresholds, horizons, sensitivity variants, bootstrap |
| 03b | Lead time |

---

## Pipeline Details

### 00_setup.R — Environment and Validation

**Purpose:** Verify the execution environment, load configuration, and validate CLIF data.

**Key Operations:**

1. **Environment verification** — Confirms the session is using the uvr project library, confirms the running R matches `.r-version`, then confirms every namespace the pipeline touches is installed. Installs nothing. A missing package stops the run in the first seconds with the command that repairs it, rather than deep into a stage that reaches for it. A checkout predating the pin has no `.r-version`, which reads as no constraint rather than as a failure.
2. **Resource detection** — Identifies cores and available RAM; sets conservative thread limits for Arrow, data.table, and collapse
3. **Configuration loading** — Reads `config_clif_oncrisk.yaml`, trims and lowercases `site_lowercase` and `file_type`, validates the site against `clif_sites.csv`, and validates `allow_sparse_o2` as an unquoted logical
4. **Artifact destination** — Sets the single `BOX_DIR` constant (`upload_to_box_v2`) that governs every write in the pipeline, and creates the six analysis subdirectories
5. **CLIF table loading** — Opens all required tables (Parquet/CSV/FST)
6. **Schema and value-domain validation** — Verifies required columns exist (case-insensitive) and expected categorical values are present

**Outputs:** Loaded CLIF tables in memory; `BOX_DIR`, `allow_sparse_o2`, and site constants; `<BOX_DIR>/env_manifest_<site>.csv`; validation report to console.

> **The environment manifest.** `00_setup.R` writes `env_manifest_<site>.csv`, a long-format table with one row per environment fact (R version, the pinned version from `.r-version`, platform, operating system, locale, time zone, library path) followed by one row per verified package and its installed version. It is written in `00` rather than assembled at the end of the run because each stage deletes every global not on its keep-list, so the objects would not survive to `run_all.R`. When one site's results diverge from the other seven, this file answers the first question asked, which is whether that site ran a different environment.

---

### 01_cohort.R — Cohort Construction

**Purpose:** Build the analytic cohort, derive outcomes and comorbidities, and export cancer-code frequencies.

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

Cancer status is determined from ICD-10-CM diagnosis codes using the bundled case-finding list (`config/icd10cm_casefinding_2023.xlsx`), restricted to its malignant-neoplasm and hematopoietic-neoplasm categories. A priority hierarchy assigns each encounter a single cancer category:

| Category | ICD-10-CM Codes | Priority (`rank_enc`) |
|----------|-----------------|-----------------------|
| Metastatic | C77–C79, C80 | 1 (highest) |
| Hematologic | C81–C96 | 2 |
| High-risk solid | C22, C25, C34 (hepatobiliary, pancreatic, lung) | 3 |
| Other solid | C18, C50, C61 | 4 |
| Other | remaining malignant/hematopoietic codes | 5 |

Exclusions: codes the case-finding list itself marks `drop` (which includes non-melanoma skin cancer, C44), codes whose definition contains "in remission" or "personal history," Z85 personal-history codes, and an `additional_drops` vector of non-specific myeloproliferative and thrombocytosis codes (D45, D47.02, D47.3, D47.9, D47.Z9, and the D72.11 family). The `additional_drops` vector is a deliberate study decision and is described in the eMethods.

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

> **`icu_01` and `wicu_01` are different variables.** `icu_01` marks any ICU stay during the encounter. `wicu_01` marks a ward-to-ICU transfer and is the outcome component. Do not substitute one for the other when reading artifacts.

#### Comorbidity Scoring

Elixhauser comorbidities via the Quan ICD-10 algorithm, summarized by the van Walraven weighted index.

**Outputs:**

| File | Description |
|------|-------------|
| `proj_tables/cohort.parquet` | Analytic cohort before score reconciliation (one row per encounter) |
| `proj_tables/outcome_times.parquet` | Five-outcome event and censoring times |
| `proj_tables/careprocess.parquet` | Vasopressor and respiratory-support events |
| `proj_tables/hid_jid_crosswalk.parquet` | Hospitalization-ID mapping |
| `<BOX_DIR>/cancer_codes_primary_<site>.csv` | Primary cancer-code frequencies (n > 5) |
| `<BOX_DIR>/hospital_types_<site>.csv` | Hospital counts by type; one type per hospital by majority vote over non-missing ADT rows |
| `<BOX_DIR>/figure_s01_flow_<site>.csv` | Inclusion flow-diagram counts (a final row is appended in `02`) |

---

### 02_scores.R — Score Derivation and Descriptive Artifacts

**Purpose:** Calculate time-varying early warning scores on the wards, reconcile the cohort against scorable encounters, and then export the encounter-level descriptive artifacts.

> **Why the descriptive artifacts live here.** Table 2, demographic missingness, and the admission-diagnosis tallies were exported from `01_cohort.R` in round one, before encounters with no calculable score were dropped. Their denominators therefore disagreed with every downstream artifact. They now run at the end of `02_scores.R`, after reconciliation, so a single reconciled cohort defines every denominator in the submission. `cohort.parquet` is rewritten here with the reconciled cohort.

#### Scores Implemented

| Score | Components | Standard Threshold |
|-------|------------|--------------------|
| **SIRS** | Temperature, heart rate, respiratory rate, WBC, pCO₂ | ≥ 2 |
| **qSOFA** | Respiratory rate, systolic BP, GCS | ≥ 2 |
| **MEWS** | Heart rate, respiratory rate, systolic BP, temperature, GCS | ≥ 5 |
| **MEWS-SF** | MEWS components + SpO₂/FiO₂ ratio | ≥ 7 |
| **NEWS** | Respiratory rate, SpO₂, supplemental O₂, temperature, systolic BP, heart rate, GCS | ≥ 5 (or any single parameter = 3) |

> **Completing NEWS.** Published NEWS (Smith 2013) includes a 2-point supplemental-oxygen item that the round-one build omitted; round two completes the score by adding it (`news_o2`): 2 points when the patient is on any supplemental oxygen (FiO₂ > 0.21 or flow > 0), 0 on room air. (This is NOT an upgrade to NEWS2 — NEWS2's distinguishing features, the alternative SpO₂ Scale 2 and the new-confusion level, are not implemented; see the manuscript's eTable 4 deviations.) The column is still named `news_total`; it now equals the round-one total plus the O₂ item. The escalation threshold (≥5) and the single-parameter "red score" rule (any component = 3) are unchanged — the O₂ item scores only 0 or 2 and cannot trigger the single-parameter rule. The O₂ item is carried onto the existing score grid by a 6-hour rolling as-of join, which leaves SIRS, qSOFA, MEWS, and MEWS-SF rows unchanged from round one. A per-site diagnostic (`news_o2_resolution`) records how each O₂ value was resolved, distinguishing values carried forward from a measurement (`carried_forward`) from values set to zero for want of one (`defaulted_zero`).

> **Sparse-oxygen guard.** If more than 90% of NEWS rows resolve to `defaulted_zero`, the run stops. That pattern means `lpm_set` and `fio2_set` are present as columns but effectively empty, and NEWS would silently degrade to the round-one build (no O₂ item) without anyone noticing. See `allow_sparse_o2` in [Configuration](#configuration).

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

#### Descriptive Artifacts (post-reconciliation)

- **Table 2** exports continuous variables as counts, sums, and sums-of-squares (`age_sum`, `age_sumsq`, `vw_sum`, `vw_sumsq`, `los_sum`, `los_sumsq`) with 2.5th and 97.5th percentiles, so the coordinating center recovers pooled means and standard deviations exactly. The categorical component carries banded age (18–39, 40–49, 50–59, 60–69, 70–79, 80+), banded van Walraven score, banded length of stay, and the five competing-risk outcome indicators. Age bands were computed but never bound into the table in round one, and are now included.
- **Admission diagnoses** are the first primary ICD-10-CM code of the earliest constituent hospitalization in each linked encounter, restricted to ED-admitted encounters, tallied at chapter and three-character stem level and suppressed at n ≤ 5. Stems are exported unmapped; the coordinating center maps them to CCSR centrally, which avoids distributing a crosswalk to eight sites.

**Outputs:**

| File | Description |
|------|-------------|
| `proj_tables/scores_full.parquet` | Complete time-varying score dataset with outcomes and strata |
| `proj_tables/scores_components.parquet` | Point-assigned components before carry-forward (input to 02b) |
| `proj_tables/news_o2_stream.parquet` | Supplemental-oxygen measurement stream before the 6h carry (input to 02b) |
| `proj_tables/vital_lab_extract.parquet` | One row per measurement, whole encounter (input to 02c) |
| `proj_tables/ward_times.parquet` | Ward-stay intervals (input to 02b/02c) |
| `proj_tables/cohort.parquet` | Rewritten with the reconciled cohort |
| `proj_tables/site_lowercase.parquet` | Site identifier, so a `03*` stage can run standalone |
| `<BOX_DIR>/table_02_cont_<site>.csv` | Continuous variable summaries (sums/sums-of-squares) |
| `<BOX_DIR>/table_02_cat_<site>.csv` | Categorical variable summaries, including age bands |
| `<BOX_DIR>/missing_demog_<site>.csv` | Demographic missingness |
| `<BOX_DIR>/figure_s01_flow_<site>.csv` | Reconciliation row appended to the flow diagram |
| `<BOX_DIR>/admission_dx_chapter-ca-<site>.csv` | Admission diagnosis by chapter (n > 5) |
| `<BOX_DIR>/admission_dx_stem-ca-<site>.csv` | Admission diagnosis by 3-char stem (n > 5) |
| `<BOX_DIR>/diagnostics/news_o2_resolution-<site>.csv` | NEWS O₂-item resolution counts by cancer status |

---

### 02b_carryforward.R — Carry-Forward Sensitivity

**Purpose:** Test the sensitivity of results to the LOCF window. Replays the post-extraction scoring pipeline at 2-, 6-, and 12-hour vital carry-forward windows from `scores_components.parquet` and `news_o2_stream.parquet` (labs stay at 12 hours throughout), then re-applies the ward restriction. The 6-hour pass must reproduce the main `02_scores.R` output exactly (asserted in QC). Emits maximum-score and 24-hour horizon-count artifacts under `cf2`, `cf6`, and `cf12` variant tokens in `sensitivity/`.

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
3. **Fixed-horizon prediction** — Score value at each observation linked to deterioration within the following 12 and 24 hours (`counts … h12`, `counts … h24`), with a cluster bootstrap at both horizons.
4. **Sensitivity analyses** — see table below.
5. **Liquid-stratum threshold block** — Standard-threshold operating characteristics (`sesp`, `ever`) for hematologic vs. solid malignancies within the cancer cohort, for the composite and nohospice outcomes.

> **The bootstrap resamples encounters, not observations.** 400 replicates draw `joined_hosp_id` values with replacement and carry every observation of a drawn encounter, weighted by its draw multiplicity, so the resampling unit matches the clustering unit and the bootstrap distribution belongs to the same point-level quantity the estimate does. Seed 2025.

| Variant | Description |
|---------|-------------|
| `main` | Primary analysis (ED admissions only) |
| `se_fullcode_only` | Restricted to full-code status (full and presume-full) |
| `se_no_ed_req` | No ED-admission requirement |
| `se_win0_96h` | Restricted to the 0–96 hour window |
| `se_one_enc_per_pt` | One randomly selected encounter per patient (seed 2025) |
| `cf2` / `cf6` / `cf12` (from 02b) | Carry-forward replays |

**Strata** span cancer status (`ca`), hematologic vs. solid (`liquid`, within cancer), and metastatic (`mets`, within solid-tumor cancer). **Outcomes** span the five keys defined in `01_cohort.R`.

---

### 03b_leadtime.R — Lead Time to Deterioration

**Purpose:** Report how much warning an alert affords. For each encounter, score, and outcome, computes the hours from the first threshold-positive score to the event, exported as binned counts plus sufficient statistics (never individual intervals), together with a crossed × event classification of every encounter (`crossclass`). Positivity uses the standard thresholds with the NEWS single-parameter rule. Because the score series is truncated before the event, every crossing precedes it and lead times are strictly positive (asserted). A QC cross-check ties the crossed-and-event counts back to `03a`'s `sesp` artifact. Outcomes: composite and nohospice. Outputs in `threshold/`.

---

### run_all.R — Orchestration and Manifest

**Purpose:** Run every stage in order, then produce a run report and an upload manifest.

- Sources each stage in the global environment, timing it and recording peak R heap use
- Writes `run_report_<site>.csv` (encounter counts, cancer fraction, outcome rate, runtime)
- Walks `upload_to_box_v2/` and writes `MANIFEST-<site>.csv` — every uploaded CSV with its row/column counts and MD5 checksum, plus site, run date, and pipeline version — so the coordinating center can confirm a complete, uncorrupted upload
- Runs a completeness check that fails loudly if any expected artifact family is absent (a silent stage failure)

> Each stage ends by deleting every global object not on its keep-list, which is why the orchestrator holds all of its own state inside `run_log`. If you add a stage or a cross-stage object, add it to the relevant keep-list, or it will vanish at the stage boundary.

---

## Output Artifacts

All outputs follow a standardized naming convention:

```
upload_to_box_v2/<family>/<artifact>[-<strata>][-<outcome>][-h<hours>][-<variant>]-<site>.csv
```

Tokens are optional and appear only where they apply: `strata` (e.g. `ca`, `liquid`, `mets`), `outcome` (`composite`, `nohospice`, `wardicu`, `warddeath`, `hospicedc`), `hours` (`h12`, `h24`), and `variant` (e.g. `se_fullcode_only`, `cf6`, `boot`).

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
│   └── counts-<strata>-<outcome>-h{12,24}-<site>.csv  # Horizon prediction counts (+ bootstrap, both horizons)
│
├── sensitivity/
│   ├── maxscores-<strata>-<outcome>-<variant>-<site>.csv
│   ├── counts-<strata>-<outcome>-h24-<variant>-<site>.csv
│   └── ...                                          # se_* and cf2/cf6/cf12 variants
│
├── diagnostics/
│   ├── news_o2_resolution-<site>.csv               # NEWS O₂-item resolution
│   ├── monitoring-ca-<site>.csv                    # Measurements per 24 ward-hours
│   ├── monitoring_bins-ca-<site>.csv               # Binned monitoring rates
│   └── missing_vlab-ca-<site>.csv                  # Vital/lab missingness by cancer status
│
├── table_02_cont_<site>.csv                        # Continuous variable summaries
├── table_02_cat_<site>.csv                         # Categorical variable summaries
├── figure_s01_flow_<site>.csv                      # Flow-diagram data
├── missing_demog_<site>.csv                        # Demographic missingness
├── cancer_codes_primary_<site>.csv                 # Cancer-code frequencies
├── hospital_types_<site>.csv                       # Hospital counts by type (eTable 1)
├── admission_dx_chapter-ca-<site>.csv              # Admission diagnosis by chapter
├── admission_dx_stem-ca-<site>.csv                 # Admission diagnosis by stem
├── env_manifest_<site>.csv                         # R version, platform, and package versions
├── run_report_<site>.csv                           # Pipeline execution summary
└── MANIFEST-<site>.csv                             # Upload manifest with checksums
```

### Files to Upload

After a successful run, upload the entire `upload_to_box_v2/` directory (including `MANIFEST-<site>.csv`) to the coordinating center. Do **not** overwrite `upload_to_box/`, which holds the frozen round-one submission.

---

## Privacy and Small-Cell Handling

The pipeline applies disclosure protections that differ by artifact class, because the disclosure risk differs. These rules are distinct and are easy to conflate:

| Artifact class | Rule |
|----------------|------|
| Score values, threshold counts, horizon counts, co-positivity (`upset`) counts | No suppression |
| AUROC | Estimability floor: computed only when there are at least 10 events **and** at least 10 non-events; otherwise returned as `NA` with `n_obs`, `n_events`, and an `estimable = FALSE` flag retained |
| Diagnosis-code tallies (cancer codes, admission diagnoses by chapter and stem) | Cells with n ≤ 5 are dropped |

The AUROC floor is an **estimability** rule, not a disclosure rule: an AUROC computed on fewer than ten events is not a usable estimate, and the underlying counts ship regardless so the coordinating center can see exactly which cells were withheld and why.

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
| "The active R library is not the uvr project library" | Pipeline launched outside uvr (RStudio Source button, or `Rscript`) | Use `uvr run code/run_all.R` or `uvr::run("code/run_all.R")` |
| "R version mismatch ... pinned to R 4.6.1" | Session running under a different R than the pin | Run `uvr r install 4.6.1`, then launch with `uvr run code/run_all.R` |
| "install ... is broken (no version response)" | R installed by uvr 0.3.x, which used a backend its own selector rejects | Update uvr with `uvr self-update`, then `uvr r uninstall 4.6.1 && uvr r install 4.6.1` |
| "Not found in the uvr project library: ..." | Library and lockfile disagree | Run `uvr sync --frozen`, then launch again |
| "R version 4.6.1 not found" or a long first sync | uvr is downloading and installing the pinned R | Expected on first run; let it finish. `uvr r list` confirms it afterward |
| `uvr sync --frozen` reports a stale lockfile | Inconsistent checkout | Run `git pull` and retry. Do **not** run plain `uvr sync` or `uvr lock` to work around this; that would re-resolve versions locally and your site would run a different package set from every other site. Contact the coordinating center |
| `uvr sync` reports missing system libraries (Linux) | Unmet system dependency | Run the `apt-get install` line uvr prints, then retry |
| Package compile failure (Windows) | Source fallback without a toolchain | Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/); confirm with `uvr doctor` |
| "Column not found" | Schema mismatch | Check column names in CLIF tables; the resolver is case-insensitive but names must exist |
| Arrow read errors | Corrupt or unreadable input files | Verify the files open outside R, or use `file_type: "csv"` |
| NEWS sparse-oxygen stop | `lpm_set` and `fio2_set` present but effectively empty | Confirm with the coordinating center, then set `allow_sparse_o2: true` |
| "Sanity check failed" | Outcome rates out of expected range | Review `cohort.parquet` for local coding issues |
| `here()` cannot find root | Launched outside the repository root, or `.here` missing | Launch from the repository root; confirm `.here` is present |
| "BOX_DIR not found" | A stage was run before `00_setup.R` | Run stages in order, or use `run_all.R` |
| Missing intermediate parquet | A `02*`/`03*` stage run out of order | Run the earlier stage first; each stage reads the previous stage's `proj_tables/` outputs |
| Memory errors | Insufficient RAM | Run stages individually inside `uvr run`, or restart R between stages |
| MANIFEST QC failure | An expected artifact family is absent | A stage failed quietly; check console output for the failing stage |

### Validation Failures

Setup performs extensive validation. Common issues:

1. **Missing columns** — Ensure all required CLIF columns are present
2. **Invalid values** — Check categorical columns against expected domains (e.g. `location_category`, `discharge_category`, `device_category`)
3. **Date formats** — Timestamps should parse as datetimes
4. **Quoted booleans** — `allow_sparse_o2` must be unquoted `true` or `false`

### Getting Help

1. Check the console output for the specific stage and error message
2. Run `uvr doctor` and capture its full output
3. Review `proj_tables/` for intermediate outputs
4. Contact the coordinating center with the error message, the `uvr doctor` output, and approximate dataset size

---

## Contributing

This is a collaborative project across CLIF consortium sites. To suggest improvements:

1. Open an issue describing the proposed change
2. For code changes, submit a pull request with a clear description, testing performed at your site, and any impact on output artifacts

Stage code explicitly (`git add code/`); never `git add -A`, to avoid pushing machine-specific config paths.

If a change adds an R package, add it with `uvr add <package>`, commit the resulting `uvr.toml` **and** `uvr.lock`, and say so in the pull request. Never edit `uvr.lock` by hand.

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
