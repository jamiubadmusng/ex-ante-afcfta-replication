
 Replication Package: An Ex-Ante Evaluation of the Economic Impact of the African Continental Free Trade Area (AfCFTA)

[![R](https://img.shields.io/badge/R-4.3.0+-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-10.xxxx/xxxxx-blue.svg)](https://doi.org/)

> **Replication materials for:** Badmus, J. O. (2026). "An Ex-Ante Evaluation of the Economic Impact of the African Continental Free Trade Area (AfCFTA)" *Erasmus Mundus Master's EGEI Dissertation*.

---

## 📋 Overview

The code in this replication package constructs the analysis files from two primary data sources: the International Trade and Production Database for Estimation (ITPD-E-R02) and the Dynamic Gravity Dataset (DGD, Release 2.1) from the U.S. International Trade Commission, using R. Nine R scripts run all of the code to generate the tables and figures in the paper. The replicator should expect the code to run for approximately 2-8 hours on a standard desktop machine.

---

## 📁 Repository Structure

```
├── README.md                 # This file
├── README.pdf                # Detailed replication guide (PDF)
├── README.tex                # LaTeX source for README.pdf
├── LICENSE                   # MIT License
├── CITATION.cff              # Citation metadata
│
├── 00_code/                  # R scripts for replication
│   ├── code_to_clean_raw_170_industries.R
│   ├── afcfta_descriptive.R
│   ├── afcfta_main_code.R
│   ├── afcfta_ge_data_preparation.R
│   ├── afcfta_ge_estimation_code_for_all_sectors_with_full_data.R
│   ├── afcfta_ge_estimation_code_for_all_sectors_with_full_data_higher_sigma.R
│   ├── afcfta_ge_estimation_code_for_struc.R
│   ├── 170_industry_level_afcfta_brdr.R
│   └── afcfta_robustness_checks.R
│
├── 01_input/                 # Input data files
│   ├── itpder2_*.rds         # ITPD-E industry-level trade data
│   ├── dgd_2_1.rds           # Dynamic Gravity Dataset
│   ├── ge_ppml_data_*.rds    # Prepared GE datasets
│   └── *.dta                 # Stata format data files
│
└── 02_output/                # Generated outputs
    ├── figures/              # PDF figures
    ├── tables/               # LaTeX tables
    └── csv/                  # CSV tables
```

---

## 💾 Data Availability

All data used in this study are **publicly available**:

**International Trade and Production Database for Estimation (ITPD-E-R02).** Data on bilateral trade flows and domestic trade across 170 industries were downloaded from the U.S. International Trade Commission. Data can be downloaded from https://www.usitc.gov/data/gravity/itpde.htm. The data are in the public domain.

Datafiles:
- `01_input/itpder2.rds` (aggregated trade data)
- `01_input/itpder2_1.rds` to `01_input/itpder2_170.rds` (industry-level trade data)

**Dynamic Gravity Dataset (DGD, Release 2.1).** Gravity covariates (distance, contiguity, colonial ties, common language, RTA membership, WTO membership, etc.) were downloaded from the U.S. International Trade Commission. Data can be downloaded from https://www.usitc.gov/data/gravity/dgd.htm.

Datafiles:
- `01_input/dgd_2_1.rds` (gravity covariates for main analysis)
- `01_input/dgd_2_1_rob.rds` (gravity covariates for robustness checks)

**Structural Gravity Database.** Manufacturing trade data from the World Trade Organization Structural Gravity Database were used for robustness checks. Data can be downloaded from https://www.wto.org/english/res_e/reser_e/structural_gravity_e.htm.

Datafile: `01_input/struc.rds`

---

## 🖥️ Computational Requirements

### Software

- **R** version 4.3.0 or later
- **RStudio** (recommended but not required)

### Required R Packages

```r
install.packages(c(
  "tidyverse",      # Data manipulation and visualization
  "fixest",         # PPML estimation with high-dimensional fixed effects
  "data.table",     # Efficient data manipulation
  "haven",          # Reading Stata files
  "readstata13",    # Reading Stata 13+ files
  "huxtable",       # Publication-quality tables
  "flextable",      # Flexible table formatting
  "msm",            # Delta method calculations
  "car",            # Companion to applied regression
  "broom",          # Tidying model outputs
  "purrr",          # Functional programming tools
  "kableExtra",     # Enhanced table formatting
  "ggplot2",        # Data visualization
  "scales",         # Scale functions for visualization
  "sf",             # Spatial data (for maps)
  "viridis"         # Color palettes
))
```

### Hardware

- **Minimum:** 8 GB RAM, quad-core processor
- **Recommended:** 16 GB RAM, 8-core processor
- **Runtime:** 2–8 hours on a standard desktop machine

---

## 📊 Description of Programs/Code

The folder `00_code/` contains nine R scripts that replicate all analyses in the paper:

| Script | Description |
|--------|-------------|
| `code_to_clean_raw_170_industries.R` | Loads 170 raw industry-level trade datasets downloaded from the USITC, filters them for the years 2000–2019, and saves the cleaned datasets as RDS files. This script should be run first if starting from raw CSV files. |
| `afcfta_descriptive.R` | Generates descriptive statistics and summary tables for intra-African trade, including trade volumes by sector, gravity covariates, and visualization of trade patterns (Figure 1 and Table 1). |
| `afcfta_main_code.R` | Estimates the main gravity model using PPML with exporter-time and importer-time fixed effects. Generates Table 2 (sectoral gravity estimates) with border effects for Total, Agriculture, Manufacturing, Mining and Energy, and Services sectors. |
| `afcfta_ge_data_preparation.R` | Prepares balanced squared datasets for General Equilibrium analysis. Creates balanced panels with imputed domestic trade for each of the four broad sectors (Agriculture, Manufacturing, Mining and Energy, Services) for the year 2019. |
| `afcfta_ge_estimation_code_for_all_sectors_with_full_data.R` | Implements the General Equilibrium PPML (GEPPML) procedure following Anderson et al. (2018) for all four sectors using ITPD-E data with σ = 7. Generates Table 3 (welfare effects by sector and country). |
| `afcfta_ge_estimation_code_for_all_sectors_with_full_data_higher_sigma.R` | Robustness check: Implements the GEPPML procedure with a higher elasticity of substitution (σ = 10). Generates Table D.5 in the Appendix. |
| `afcfta_ge_estimation_code_for_struc.R` | Robustness check: Implements the GEPPML procedure using manufacturing trade data from the Structural Gravity Database with multiple elasticity values (σ = 4, 5, 7, 10). Generates Table D.7 and Figure D.1 (welfare effects world map). |
| `170_industry_level_afcfta_brdr.R` | Estimates industry-level border effects for all 170 industries in the ITPD-E database. Generates Table D.2 (industry-level estimates with full sample) and Table D.4 (estimates excluding zero trade flows), along with Figure 3 (distribution of border effects by sector). |
| `afcfta_robustness_checks.R` | Robustness checks for the main gravity estimates using the Structural Gravity Database. Tests sensitivity to additional controls (WTO, PTA, EIA membership). Generates Table D.6. |


---

## 🚀 Instructions to Replicators

1. **Maintain the folder structure** of the replication package with three folders: `00_code/`, `01_input/`, and `02_output/`. The output folder should contain subfolders: `figures/`, `tables/`, `csv/`, and `rds/`.

2. **Download the required data files** from the sources listed above and place them in the `00_input/` folder in the appropriate format (RDS or DTA).

3. **Edit the `setwd()` command** at the beginning of each R script to point to the location of the replication package on your machine.

4. **Install all required R packages** (see Software Requirements above).

5. **Run the scripts in the following order:**
   1. `code_to_clean_raw_170_industries.R` (only if starting from raw CSV files)
   2. `afcfta_descriptive.R`
   3. `afcfta_main_code.R`
   4. `afcfta_ge_data_preparation.R`
   5. `afcfta_ge_estimation_code_for_all_sectors_with_full_data.R`
   6. `afcfta_ge_estimation_code_for_all_sectors_with_full_data_higher_sigma.R`
   7. `afcfta_ge_estimation_code_for_struc.R`
   8. `170_industry_level_afcfta_brdr.R`
   9. `afcfta_robustness_checks.R`

6. **Output files** (tables and figures) will be saved to the `02_ output/` folder.

---

## 📈 List of Tables and Figures

The provided code reproduces:

- [x] All numbers provided in text in the paper.
- [x] All tables and figures in the paper.

### Main Paper

| Output | Program | Output File |
|--------|---------|-------------|
| Figure 2 | `afcfta_descriptive.R` | `output/figures/figure_2.pdf` |
| Figure 3 | `170_industry_level_afcfta_brdr.R` | `output/figures/figure_3.pdf` |
| Table 1 | `afcfta_descriptive.R` | `output/tables/table_1a.tex` & `table_1b.tex` |
| Table 2 | `afcfta_main_code.R` | `output/tables/table_2.tex` |
| Table 3 | `afcfta_ge_estimation_code_for_all_sectors_with_full_data.R` | `output/tables/welfare_effects_table_full_data.tex` |

### Appendix

| Output | Program | Output File |
|--------|---------|-------------|
| Table D.2 | `170_industry_level_afcfta_brdr.R` | `output/csv/table_a_3.csv` |
| Table D.3 | `afcfta_main_code.R` | `output/tables/table_a_4.tex` |
| Table D.4 | `170_industry_level_afcfta_brdr.R` | `output/csv/table_a_5.csv` |
| Table D.5 | `afcfta_ge_estimation_code_for_all_sectors_with_full_data_higher_sigma.R` | `output/tables/welfare_effects_table_full_data_higher_sigma.tex` |
| Table D.6 | `afcfta_robustness_checks.R` | `output/tables/table_a_7.tex` |
| Table D.7 | `afcfta_ge_estimation_code_for_struc.R` | `output/tables/welfare_effects_table_struc_data.tex` |
| Figure D.1 | `afcfta_ge_estimation_code_for_struc.R` | `output/figures/welfare_effects_world_map_struc_viridis.pdf` |

---


## 📖 Citation

If you use this code or data, please cite:

```bibtex
@article{badmus2026afcfta,
  title={An Ex-Ante Evaluation of the Economic Impact of the African Continental Free Trade Area (AfCFTA)},
  author={Badmus, Jamiu Olamilekan},
  year={2026},
  volume={},
  pages={},
  doi={}
}
```

---

## 📧 Contact

**Jamiu Olamilekan Badmus**
[Personal Website](https://sites.google.com/view/jamiu-olamilekan-badmus/)
Email: [jamiubadmus001@gmail.com]

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

The data used in this project are subject to their original licenses from USITC and WTO.

---

## 🙏 Acknowledgments

- This replication package follows the [AEA Data Editor's guidelines](https://aeadataeditor.github.io/aea-de-guidance/) and the [Social Science Data Editors' template](https://social-science-data-editors.github.io/template_README/)

---

*Last updated: January 2026*
