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
