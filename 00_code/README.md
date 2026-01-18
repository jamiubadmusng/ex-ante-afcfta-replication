# Code Documentation

This folder contains all R scripts for replicating the analysis.

## Scripts Overview

| Script | Purpose | Runtime |
|--------|---------|---------|
| `code_to_clean_raw_170_industries.r` | Clean raw USITC data | ~5 min |
| `afcfta_descriptive.R` | Descriptive statistics & Figure 1-2 | ~10 min |
| `afcfta_main_code.R` | Main PPML gravity estimation (Table 2) | ~15 min |
| `afcfta_ge_data_preparation.R` | Prepare balanced GE datasets | ~20 min |
| `afcfta_ge_estimation_code_for_all_sectors_with_full_data.R` | GE welfare analysis, σ=7 (Table 3) | ~60 min |
| `afcfta_ge_estimation_code_for_all_sectors_with_full_data_higher_sigma.r` | GE welfare analysis, σ=10 (Table D.5) | ~60 min |
| `afcfta_ge_estimation_code_for_struc.r` | GE with Structural Gravity (Table D.7) | ~90 min |
| `170_industry_level_afcfta_brdr.R` | Industry-level estimates (Table D.2, D.4) | ~45 min |
| `afcfta_robustness_checks.R` | Robustness checks (Table D.6) | ~15 min |

## Execution Order

Run the scripts in the following order:

1. `code_to_clean_raw_170_industries.r` (only if starting from raw CSV files)
2. `afcfta_descriptive.R`
3. `afcfta_main_code.R`
4. `afcfta_ge_data_preparation.R`
5. `afcfta_ge_estimation_code_for_all_sectors_with_full_data.R`
6. `afcfta_ge_estimation_code_for_all_sectors_with_full_data_higher_sigma.r`
7. `afcfta_ge_estimation_code_for_struc.r`
8. `170_industry_level_afcfta_brdr.R`
9. `afcfta_robustness_checks.R`

## Important Notes

- **Working directory:** Each script sets `setwd()` at the beginning. Update this path to match your local setup.
- **Dependencies:** All scripts require the packages listed in the main README.
- **Memory:** GE estimation scripts may require 8+ GB RAM.

---

*Last updated: January 2026*
