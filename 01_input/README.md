# Data Documentation

This folder contains the data files used in the analysis.

## Data Sources

### 1. International Trade and Production Database for Estimation (ITPD-E-R02)

**Source:** U.S. International Trade Commission  
**URL:** https://www.usitc.gov/data/gravity/itpde.htm  
**License:** Public Domain  

The ITPD-E-R02 provides consistent bilateral trade and domestic trade data for 265 countries, 170 industries, and the years 1986–2019.

**Files:**
- `itpder2.rds` — Aggregated trade data (all sectors combined)
- `itpder2_1.rds` to `itpder2_170.rds` — Industry-level trade data (one file per industry)

### 2. Dynamic Gravity Dataset (DGD, Release 2.1)

**Source:** U.S. International Trade Commission  
**URL:** https://www.usitc.gov/data/gravity/dgd.htm  
**License:** Public Domain  

The DGD provides gravity covariates including bilateral distance, contiguity, colonial history, common language, and trade agreement membership.

**Files:**
- `dgd_2_1.rds` — Gravity covariates for main analysis (2000–2019)
- `dgd_2_1_rob.rds` — Extended gravity covariates for robustness checks (includes WTO, PTA, EIA)
- `release_2.1_2000_2019.zip` — Raw DGD data files

### 3. Structural Gravity Database

**Source:** World Trade Organization  
**URL:** https://www.wto.org/english/res_e/reser_e/structural_gravity_e.htm  
**License:** Public Domain  

Manufacturing trade data used for robustness checks.

**Files:**
- `struc.rds` — Processed Structural Gravity data
- `struc_2016_balanced.dta` — Balanced panel for GE analysis

---

*Last updated: January 2026*
