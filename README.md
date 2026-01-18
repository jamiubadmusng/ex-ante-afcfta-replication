
 Replication Package: Ex-Ante Economic Impacts of the African Continental Free Trade Area (AfCFTA)

[![R](https://img.shields.io/badge/R-4.3.0+-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-10.xxxx/xxxxx-blue.svg)](https://doi.org/)

> **Replication materials for:** Badmus, J. O. (2026). "Ex-Ante Economic Impacts of the African Continental Free Trade Area: A Structural Gravity Analysis." *International Economics*.

---

## 📋 Overview

This repository contains the replication package for the empirical analysis of the ex-ante economic impacts of the African Continental Free Trade Area (AfCFTA) agreement. The analysis uses structural gravity models with Poisson Pseudo-Maximum Likelihood (PPML) estimation and General Equilibrium (GE) counterfactual simulations to quantify:

1. **Border effects** of the AfCFTA on intra-African trade across sectors
2. **Welfare effects** (% change in real GDP) from eliminating intra-African trade barriers
3. **Industry-level heterogeneity** in trade creation effects across 170 industries


### Key Findings

- Intra-African trade faces significant border frictions equivalent to **21–24% lower trade flows**
- AfCFTA implementation could increase African welfare by **27–108%** depending on the sector
- Substantial heterogeneity exists across industries, with services showing the largest gains

---

## 📁 Repository Structure

```
├── README.md                    # This file
├── README.pdf                   # Detailed replication guide (PDF)
├── LICENSE                      # MIT License
├── CITATION.cff                 # Citation file for this repository
├── .gitignore                   # Git ignore file
│
├── code/                        # R scripts for analysis
│   ├── 00_master.R              # Master script (run this to replicate all)
│   ├── 01_clean_industry_data.R # Data cleaning for 170 industries
│   ├── 02_descriptive.R         # Descriptive statistics and figures
│   ├── 03_main_estimation.R     # Main gravity estimation (Table 2)
│   ├── 04_ge_data_prep.R        # GE data preparation
│   ├── 05_ge_estimation.R       # GE welfare analysis (Table 3)
│   ├── 06_ge_robustness.R       # GE with alternative sigma
│   ├── 07_ge_struc_data.R       # GE with Structural Gravity data
│   ├── 08_industry_analysis.R   # Industry-level border effects
│   └── 09_robustness_checks.R   # Additional robustness checks
│
├── data/                        # Data files
│   ├── raw/                     # Original data (see data/README.md)
│   └── processed/               # Processed analysis datasets
│
└── output/                      # Results
    ├── figures/                 # Publication-ready figures
    ├── tables/                  # LaTeX and CSV tables
    └── logs/                    # Estimation logs
```

---

## 💾 Data Availability

All data used in this study are **publicly available**:

| Dataset | Source | URL |
|---------|--------|-----|
| ITPD-E-R02 | U.S. International Trade Commission | [usitc.gov/data/gravity/itpde.htm](https://www.usitc.gov/data/gravity/itpde.htm) |
| DGD Release 2.1 | U.S. International Trade Commission | [usitc.gov/data/gravity/dgd.htm](https://www.usitc.gov/data/gravity/dgd.htm) |
| Structural Gravity Database | World Trade Organization | [wto.org](https://www.wto.org/english/res_e/reser_e/structural_gravity_e.htm) |

See [`data/README.md`](data/README.md) for detailed data documentation.

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

## 🚀 Replication Instructions

### Quick Start

1. **Clone the repository:**
   ```bash
   git clone https://github.com/[username]/afcfta-replication.git
   cd afcfta-replication
   ```

2. **Download the data** from the sources listed above and place in `data/raw/`

3. **Open RStudio** and set the working directory to the repository root

4. **Run the master script:**
   ```r
   source("code/00_master.R")
   ```

### Step-by-Step Replication

If you prefer to run scripts individually:

| Order | Script | Description | Approx. Time |
|-------|--------|-------------|--------------|
| 1 | `01_clean_industry_data.R` | Clean raw industry data | 5 min |
| 2 | `02_descriptive.R` | Descriptive statistics | 10 min |
| 3 | `03_main_estimation.R` | Main gravity estimation | 15 min |
| 4 | `04_ge_data_prep.R` | Prepare GE datasets | 20 min |
| 5 | `05_ge_estimation.R` | GE welfare analysis | 60 min |
| 6 | `06_ge_robustness.R` | GE with σ = 10 | 60 min |
| 7 | `07_ge_struc_data.R` | GE with Structural Gravity | 90 min |
| 8 | `08_industry_analysis.R` | Industry-level analysis | 45 min |
| 9 | `09_robustness_checks.R` | Robustness checks | 15 min |

---

## 📊 Output Files

### Main Paper

| Output | Script | File |
|--------|--------|------|
| Figure 2 | `02_descriptive.R` | `output/figures/figure_2.pdf` |
| Figure 3 | `08_industry_analysis.R` | `output/figures/figure_3.pdf` |
| Table 1 | `02_descriptive.R` | `output/tables/table_1a.tex`, `table_1b.tex` |
| Table 2 | `03_main_estimation.R` | `output/tables/table_2.tex` |
| Table 3 | `05_ge_estimation.R` | `output/tables/welfare_effects_table_full_data.tex` |

### Appendix

| Output | Script | File |
|--------|--------|------|
| Table D.2 | `08_industry_analysis.R` | `output/csv/table_a_3.csv` |
| Table D.3 | `03_main_estimation.R` | `output/tables/table_a_4.tex` |
| Table D.4 | `08_industry_analysis.R` | `output/csv/table_a_5.csv` |
| Table D.5 | `06_ge_robustness.R` | `output/tables/welfare_effects_table_full_data_higher_sigma.tex` |
| Table D.6 | `09_robustness_checks.R` | `output/tables/table_a_7.tex` |
| Table D.7 | `07_ge_struc_data.R` | `output/tables/welfare_effects_table_struc_data.tex` |
| Figure D.1 | `07_ge_struc_data.R` | `output/figures/welfare_effects_world_map_struc_viridis.pdf` |

---

## 📖 Citation

If you use this code or data, please cite:

```bibtex
@article{badmus2026afcfta,
  title={Ex-Ante Economic Impacts of the African Continental Free Trade Area: 
         A Structural Gravity Analysis},
  author={Badmus, Jamiu Olamilekan},
  journal={International Economics},
  year={2026},
  volume={},
  pages={},
  doi={}
}
```

---

## 📧 Contact

**Jamiu Olamilekan Badmus**  
Email: [your.email@institution.edu]  
ORCID: [0000-0000-0000-0000]

---

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

The data used in this project are subject to their original licenses from USITC and WTO.

---

## 🙏 Acknowledgments

- This replication package follows the [AEA Data Editor's guidelines](https://aeadataeditor.github.io/aea-de-guidance/) and the [Social Science Data Editors' template](https://social-science-data-editors.github.io/template_README/)
- General Equilibrium PPML methodology based on [Anderson et al. (2018)](https://doi.org/10.1016/j.jinteco.2018.06.002)
- `fixest` package by [Bergé (2018)](https://lrberge.github.io/fixest/)

---

*Last updated: January 2026*
