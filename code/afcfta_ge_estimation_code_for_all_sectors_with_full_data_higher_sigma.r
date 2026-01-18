################################################################################
# R code for General Equilibrium PPML Analysis of the AfCFTA Agreement
# This code translates the original Stata GE PPML code by Anderson et al. (2018)
# for the analysis of abolishing borders between AfCFTA countries using
# bilateral trade data from ITPD-E-R02 with a higher elasticity of substitution.
# Author: Jamiu Olamilekan Badmus
################################################################################


################################################################################
# Section 1: General Equilibrium PPML Analysis for Agricultural Sector ---------
################################################################################

# Clear workspace
rm(list = ls())

# Load required packages
# packages <- c("data.table", "fixest", "dplyr", "tidyr", "readstata13", "ggplot2")
# new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
# if(length(new_packages)) install.packages(new_packages)

library(data.table)
library(fixest)
library(dplyr)
library(tidyr)
library(readstata13)
library(ggplot2)

# Set working directory
setwd("C:/Users/muham/AfCFTA_UNU_CRIS")

# Set parameters
sigma <- 10 # Elasticity of substitution


################################################################################
# I. Prepare Data for Agricultural Sector
################################################################################

# Load data
data <- read.dta13("01_input/afcfta_2019_agri_balanced.dta")
dt <- as.data.table(data)

# 1. Create aggregate variables
dt[, output := sum(trade), by = "exporter"]
dt[, expndr := sum(trade), by = "importer"]

# 2. Choose a country for reference group (Germany)
dt[exporter == "DEU", exporter := "ZZZ"]
dt[importer == "DEU", importer := "ZZZ"]
dt[importer == "ZZZ", expndr_deu0 := expndr]
dt[, expndr_deu := mean(expndr_deu0, na.rm = TRUE)]

# Save processed data
saveRDS(dt, "01_input/ge_ppml_data_agri.rds")

################################################################################
# II. GE Analysis in R
################################################################################

# Load prepared data
dt <- readRDS("01_input/ge_ppml_data_agri.rds")

#############################
# Step 1: "Baseline" Scenario
#############################

# Step 1.a: Estimate "Baseline" Gravity
# Create country dummy variables (for exporter and importer fixed effects)
countries <- unique(dt$exporter)
NoC <- length(countries)

# Baseline PPML with fixest
# Using fepois for PPML with high-dimensional fixed effects
baseline_model <- fepois(
  trade ~ cntg + col + lang + dist + rta + afcfta_brdr + afcfta_row_brdr + row_to_row | exporter + importer,
  data = dt
)

# Save the model summary
baseline_summary <- summary(baseline_model)

# Extract coefficients
coeffs <- coef(baseline_model)
CNTG_est <- coeffs["cntg"]
COL_est <- coeffs["col"]
LANG_est <- coeffs["lang"]
DIST_est <- coeffs["dist"]
RTA_est <- coeffs["rta"]
AfCFTA_BRDR_est <- coeffs["afcfta_brdr"]
AfCFTA_ROW_BRDR_est <- coeffs["afcfta_row_brdr"]
ROW_to_ROW_est <- coeffs["row_to_row"]

# Predict trade in the baseline
dt[, trade_bsln := predict(baseline_model, newdata = dt, type = "response")]

# Create baseline and counterfactual trade costs
dt[, t_ij_bsln := exp(
  DIST_est * dist +
  CNTG_est * cntg +
  COL_est * col +
  LANG_est * lang +
  RTA_est * rta +
  AfCFTA_BRDR_est * afcfta_brdr +
  AfCFTA_ROW_BRDR_est * afcfta_row_brdr +
  ROW_to_ROW_est * row_to_row
)]

# Counterfactual: set afcfta_brdr = 0, keep everything else the same
dt[, t_ij_ctrf := exp(
  DIST_est * dist +
  CNTG_est * cntg +
  COL_est * col +
  LANG_est * lang +
  RTA_est * rta +
  AfCFTA_BRDR_est * 0 +
  AfCFTA_ROW_BRDR_est * afcfta_row_brdr +
  ROW_to_ROW_est * row_to_row
)]

# Keep domestic trade costs at baseline
dt[exporter == importer, t_ij_ctrf := t_ij_bsln]
dt[, t_ij_ctrf_1 := log(t_ij_ctrf)]

# Step 1.b: Construct "Baseline" GE Indexes
# Extract fixed effects
fe_exporter <- fixef(baseline_model)$exporter
fe_importer <- fixef(baseline_model)$importer

# Convert fixed effects to data table
exp_fe_dt <- data.table(exporter = names(fe_exporter), exp_fe = exp(fe_exporter))
imp_fe_dt <- data.table(importer = names(fe_importer), imp_fe = exp(fe_importer))

# Merge fixed effects with main data
dt <- merge(dt, exp_fe_dt, by = "exporter", all.x = TRUE)
dt <- merge(dt, imp_fe_dt, by = "importer", all.x = TRUE)

# Calculate outward and inward multilateral resistance terms
dt[, all_exp_fes_0 := exp_fe]
dt[, all_imp_fes_0 := imp_fe]

# Calculate importer-specific exporter fixed effects
dt[exporter == importer, temp := all_exp_fes_0]
dt[, all_exp_fes_0_TB := mean(temp, na.rm = TRUE), by = "importer"]
dt[, temp := NULL]

# Equation (7): Outward Multilateral Resistance
dt[, omr_bsln := output * expndr_deu / all_exp_fes_0]

# Equation (8): Inward Multilateral Resistance
dt[, imr_bsln := expndr / (all_imp_fes_0 * expndr_deu)]

# Real GDP (baseline)
dt[exporter == importer, rGDP_bsln_temp := output / (imr_bsln^(1/(1-sigma)))]
dt[, rGDP_bsln := sum(rGDP_bsln_temp, na.rm = TRUE), by = "exporter"]

# Domestic absorption share (acr)
dt[exporter == importer, acr_bsln := trade_bsln / expndr]

# Bilateral exports & totals
dt[exporter != importer, exp_bsln := trade_bsln]
dt[exporter == importer, exp_bsln_acr := trade_bsln]
dt[, tot_exp_bsln := sum(exp_bsln, na.rm = TRUE), by = "exporter"]

################################
# Step 2: "Conditional" Scenario
################################

# Step 2.a: Estimate "Conditional" Gravity
# PPML with offset = log(counterfactual trade costs)
conditional_model <- fepois(
  trade ~ 1 | exporter + importer,
  data = dt,
  offset = ~t_ij_ctrf_1
)

# Predict trade under conditional GE
dt[, trade_cndl := predict(conditional_model, newdata = dt, type = "response")]

# Step 2.b: Construct "Conditional" GE Indexes
# Extract conditional fixed effects
fe_exporter_cndl <- fixef(conditional_model)$exporter
fe_importer_cndl <- fixef(conditional_model)$importer

# Convert fixed effects to data table
exp_fe_cndl_dt <- data.table(exporter = names(fe_exporter_cndl), exp_fe_cndl = exp(fe_exporter_cndl))
imp_fe_cndl_dt <- data.table(importer = names(fe_importer_cndl), imp_fe_cndl = exp(fe_importer_cndl))

# Merge fixed effects with main data
dt <- merge(dt, exp_fe_cndl_dt, by = "exporter", all.x = TRUE)
dt <- merge(dt, imp_fe_cndl_dt, by = "importer", all.x = TRUE)

# Calculate outward and inward multilateral resistance terms
dt[, all_exp_fes_1 := exp_fe_cndl]
dt[, all_imp_fes_1 := imp_fe_cndl]

# Equation (7): Outward Multilateral Resistance (conditional)
dt[, omr_cndl := output * expndr_deu / all_exp_fes_1]

# Equation (8): Inward Multilateral Resistance (conditional)
dt[, imr_cndl := expndr / (all_imp_fes_1 * expndr_deu)]

# Exports and totals (conditional)
dt[exporter != importer, exp_cndl := trade_cndl]
dt[, tot_exp_cndl := sum(exp_cndl, na.rm = TRUE), by = "exporter"]

# % change in total exports
dt[, tot_exp_cndl_ch := (tot_exp_cndl - tot_exp_bsln) / tot_exp_bsln * 100]

# Real GDP (conditional)
dt[exporter == importer, rGDP_cndl_temp := output / (imr_cndl^(1/(1-sigma)))]
dt[, rGDP_cndl := sum(rGDP_cndl_temp, na.rm = TRUE), by = "exporter"]

###################################
# Step 3: "Full Endowment" Scenario
###################################

# Step 3.a: Estimate "Full Endowment" Gravity
# Initialize variables for the iterative process
dt[, trade_1_pred := trade_cndl]
dt[, output_bsln := output]
dt[, expndr_bsln := expndr]
dt[exporter == importer, phi := expndr / output]

# Calculate initial values
dt[exporter == importer, temp := all_exp_fes_0]
dt[, all_exp_fes_0_imp := mean(temp, na.rm = TRUE), by = "importer"]
dt[, temp := NULL]

dt[exporter == importer, expndr_temp_1 := phi * output]
dt[, expndr_1 := mean(expndr_temp_1, na.rm = TRUE), by = "importer"]
dt[, expndr_temp_1 := NULL]

dt[importer == "ZZZ", expndr_deu01_1 := expndr_1]
dt[, expndr_deu_1 := mean(expndr_deu01_1, na.rm = TRUE)]

dt[exporter == importer, temp := all_exp_fes_1]
dt[, all_exp_fes_1_imp := mean(temp, na.rm = TRUE), by = "importer"]
dt[, temp := NULL]

# Initial values for factory-gate prices and multilateral resistance terms
dt[, p_full_exp_0 := 0]
dt[, p_full_exp_1 := (all_exp_fes_1 / all_exp_fes_0)^(1/(1-sigma))]
dt[, p_full_imp_1 := (all_exp_fes_1_imp / all_exp_fes_0_imp)^(1/(1-sigma))]
dt[, imr_full_1 := expndr_1 / (all_imp_fes_1 * expndr_deu_1)]
dt[, imr_full_ch_1 := 1]
dt[, omr_full_1 := output * expndr_deu_1 / all_exp_fes_1]
dt[, omr_full_ch_1 := 1]

# Set convergence parameters
max_iterations <- 30
diff_all_exp_fes_sd <- 1
diff_all_exp_fes_max <- 1
tolerance <- 0.001
i <- 3
iteration <- 0

# Define a function to check for duplicated variables and clean them up
check_and_remove_vars <- function(dt, pattern) {
  vars_to_remove <- grep(pattern, names(dt), value = TRUE)
  if (length(vars_to_remove) > 0) {
    dt[, (vars_to_remove) := NULL]
  }
  return(dt)
}

# Iterative process
cat("Starting iterative process for full endowment calculation...\n")

while ((diff_all_exp_fes_sd > tolerance || diff_all_exp_fes_max > tolerance) && iteration < max_iterations) {
  iteration <- iteration + 1
  cat(sprintf("Iteration %d of maximum %d\n", iteration, max_iterations))
  
  # Clean up variables before recreation
  current_iter <- i-1
  patterns <- c(
    paste0("trade_", current_iter),
    paste0("all_exp_fes_", current_iter),
    paste0("all_imp_fes_", current_iter),
    paste0("output_", current_iter),
    paste0("expndr_check_", current_iter),
    paste0("expndr_deu0_", current_iter),
    paste0("expndr_deu_", current_iter),
    paste0("all_exp_fes_", current_iter, "_imp"),
    paste0("p_full_exp_", current_iter),
    paste0("p_full_imp_", current_iter),
    paste0("omr_full_", current_iter),
    paste0("omr_full_ch_", current_iter),
    paste0("expndr_temp_", current_iter),
    paste0("expndr_", current_iter),
    paste0("imr_full_", current_iter),
    paste0("imr_full_ch_", current_iter),
    paste0("diff_p_full_exp_", current_iter),
    paste0("trade_", current_iter, "_pred")
  )
  
  for (pattern in patterns) {
    dt <- check_and_remove_vars(dt, pattern)
  }
  
  # Equation (14)
  dt[, paste0("trade_", current_iter) := get(paste0("trade_", i-2, "_pred")) * 
      get(paste0("p_full_exp_", i-2)) * get(paste0("p_full_imp_", i-2)) / 
      (get(paste0("omr_full_ch_", i-2)) * get(paste0("imr_full_ch_", i-2)))]
  
  # Estimate new PPML model
  tryCatch({
    formula_str <- paste0("trade_", current_iter, " ~ 1 | exporter + importer")
    full_model <- fepois(
      as.formula(formula_str),
      data = dt,
      offset = ~t_ij_ctrf_1
    )
    
    # Predict trade flows
    dt[, paste0("trade_", current_iter, "_pred") := predict(full_model, newdata = dt, type = "response")]
    
    # Extract fixed effects
    fe_exporter_full <- fixef(full_model)$exporter
    fe_importer_full <- fixef(full_model)$importer
    
    # Process fixed effects
    dt[, paste0("all_exp_fes_", current_iter) := 0]
    dt[, paste0("all_imp_fes_", current_iter) := 0]
    
    # Apply exporter fixed effects
    for (country in names(fe_exporter_full)) {
      dt[exporter == country, paste0("all_exp_fes_", current_iter) := exp(fe_exporter_full[country])]
    }
    
    # Apply importer fixed effects
    for (country in names(fe_importer_full)) {
      dt[importer == country, paste0("all_imp_fes_", current_iter) := exp(fe_importer_full[country])]
    }
    
    # Update output
    dt[, paste0("output_", current_iter) := sum(get(paste0("trade_", current_iter, "_pred"))), by = "exporter"]
    
    # Update expenditure
    dt[, paste0("expndr_check_", current_iter) := sum(get(paste0("trade_", current_iter, "_pred"))), by = "importer"]
    dt[importer == "ZZZ", paste0("expndr_deu0_", current_iter) := get(paste0("expndr_check_", current_iter))]
    dt[, paste0("expndr_deu_", current_iter) := mean(get(paste0("expndr_deu0_", current_iter)), na.rm = TRUE)]
    
    # Get importer-specific exporter fixed effects
    dt[exporter == importer, temp := get(paste0("all_exp_fes_", current_iter))]
    dt[, paste0("all_exp_fes_", current_iter, "_imp") := mean(temp, na.rm = TRUE), by = "importer"]
    dt[, temp := NULL]
    
    # Update factory-gate prices
    dt[, paste0("p_full_exp_", current_iter) := 
        ((get(paste0("all_exp_fes_", current_iter)) / get(paste0("all_exp_fes_", i-2))) / 
           (get(paste0("expndr_deu_", current_iter)) / get(paste0("expndr_deu_", i-2))))^(1/(1-sigma))]
    
    dt[, paste0("p_full_imp_", current_iter) := 
        ((get(paste0("all_exp_fes_", current_iter, "_imp")) / get(paste0("all_exp_fes_", i-2, "_imp"))) / 
           (get(paste0("expndr_deu_", current_iter)) / get(paste0("expndr_deu_", i-2))))^(1/(1-sigma))]
    
    # Equation (7) - Update outward multilateral resistance
    dt[, paste0("omr_full_", current_iter) := get(paste0("output_", current_iter)) / get(paste0("all_exp_fes_", current_iter))]
    dt[, paste0("omr_full_ch_", current_iter) := get(paste0("omr_full_", current_iter)) / get(paste0("omr_full_", i-2))]
    
    # Update expenditure
    dt[exporter == importer, paste0("expndr_temp_", current_iter) := phi * get(paste0("output_", current_iter))]
    dt[, paste0("expndr_", current_iter) := mean(get(paste0("expndr_temp_", current_iter)), na.rm = TRUE), by = "importer"]
    
    # Equation (8) - Update inward multilateral resistance
    dt[, paste0("imr_full_", current_iter) := get(paste0("expndr_", current_iter)) / 
         (get(paste0("all_imp_fes_", current_iter)) * get(paste0("expndr_deu_", current_iter)))]
    dt[, paste0("imr_full_ch_", current_iter) := get(paste0("imr_full_", current_iter)) / get(paste0("imr_full_", i-2))]
    
    # Calculate convergence criteria
    dt[, paste0("diff_p_full_exp_", current_iter) := get(paste0("p_full_exp_", i-2)) - get(paste0("p_full_exp_", i-3))]
    
    # Calculate convergence statistics
    diff_stats <- dt[, .(sd = sd(get(paste0("diff_p_full_exp_", current_iter)), na.rm = TRUE),
                          max = max(abs(get(paste0("diff_p_full_exp_", current_iter))), na.rm = TRUE))]
    diff_all_exp_fes_sd <- diff_stats$sd
    diff_all_exp_fes_max <- diff_stats$max
    
    cat(sprintf("Convergence stats - SD: %f, Max: %f\n", diff_all_exp_fes_sd, diff_all_exp_fes_max))
    
  }, error = function(e) {
    cat("Error in iteration", iteration, ":", e$message, "\n")
    cat("Using values from previous iteration\n")
    
    # Create placeholder with previous values
    dt[, paste0("trade_", current_iter, "_pred") := get(paste0("trade_", i-2, "_pred"))]
  })
  
  # Increment counter for next iteration
  i <- i + 1
}

cat(sprintf("Convergence achieved after %d iterations\n", iteration))

# Get the final iteration number for use in subsequent calculations
final_iter <- i - 2

# Step 3.b: Construct "Full Endowment" GE Indexes
# Calculate p^c/p
dt[, c("output", "expndr", "expndr_deu0") := NULL]
dt[, expndr_deu_bsln := expndr_deu]

dt[, output := sum(get(paste0("trade_", final_iter, "_pred"))), by = "exporter"]
dt[exporter == importer, expndr_temp := phi * output]
dt[, expndr := mean(expndr_temp, na.rm = TRUE), by = "importer"]
dt[importer == "ZZZ", expndr_deu0 := expndr]
dt[, expndr_deu := mean(expndr_deu0, na.rm = TRUE)]

# Calculate p_full
dt[, p_full := ((get(paste0("all_exp_fes_", final_iter)) / all_exp_fes_0) / 
                 (expndr_deu / expndr_deu_bsln))^(1/(1-sigma))]

# Calculate output_full
dt[, output_full := p_full * output_bsln]

# Equation (7): Outward Multilateral Resistance
dt[, omr_full := output_full * expndr_deu / get(paste0("all_exp_fes_", final_iter))]

# Equation (8): Inward Multilateral Resistance
dt[, imr_full := expndr / (get(paste0("all_imp_fes_", final_iter)) * expndr_deu)]

# Real GDP (full endowment)
dt[exporter == importer, rGDP_full_temp := p_full * output_bsln / (imr_full^(1/(1-sigma)))]
dt[, rGDP_full := sum(rGDP_full_temp, na.rm = TRUE), by = "exporter"]

# Expenditure (full endowment)
dt[exporter == importer, expndr_full_temp := phi * output_full]
dt[, expndr_full := mean(expndr_full_temp, na.rm = TRUE), by = "importer"]

# Bilateral trade at full endowment (using counterfactual costs)
dt[, trade_full := (output_full * expndr_full * t_ij_ctrf) / (imr_full * omr_full)]
dt[exporter != importer, exp_full := trade_full]
dt[, tot_exp_full := sum(exp_full, na.rm = TRUE), by = "exporter"]
dt[, tot_exp_full_ch := (tot_exp_full - tot_exp_bsln) / tot_exp_bsln * 100]

# Save results
saveRDS(dt, "02_output/rds/full_static_all_agri_higher_sigma.rds")

# Prepare all indexes at the country level
# IMR indexes
imr_dt <- dt[, .(imr_full = mean(imr_full, na.rm = TRUE),
                 imr_bsln = mean(imr_bsln, na.rm = TRUE),
                 imr_cndl = mean(imr_cndl, na.rm = TRUE)), by = importer]
setnames(imr_dt, "importer", "country")
imr_dt[, imr_full_ch := (imr_full - imr_bsln) / imr_bsln * 100]
imr_dt[, imr_cndl_ch := (imr_cndl - imr_bsln) / imr_bsln * 100]
saveRDS(imr_dt, "02_output/rds/imrs_all_agri_higher_sigma.rds")

# OMR and other indexes
omr_dt <- dt[, .(omr_full = mean(omr_full, na.rm = TRUE),
                 omr_cndl = mean(omr_cndl, na.rm = TRUE),
                 omr_bsln = mean(omr_bsln, na.rm = TRUE),
                 rGDP_full = mean(rGDP_full, na.rm = TRUE),
                 rGDP_cndl = mean(rGDP_cndl, na.rm = TRUE),
                 rGDP_bsln = mean(rGDP_bsln, na.rm = TRUE),
                 tot_exp_full = mean(tot_exp_full, na.rm = TRUE),
                 tot_exp_cndl = mean(tot_exp_cndl, na.rm = TRUE),
                 tot_exp_bsln = mean(tot_exp_bsln, na.rm = TRUE),
                 p_full = mean(p_full, na.rm = TRUE),
                 output_bsln = mean(output_bsln, na.rm = TRUE),
                 acr_bsln = mean(acr_bsln, na.rm = TRUE)), by = exporter]
setnames(omr_dt, "exporter", "country")
omr_dt[, omr_full_ch := (omr_full - omr_bsln) / omr_bsln * 100]
omr_dt[, omr_cndl_ch := (omr_cndl - omr_bsln) / omr_bsln * 100]
omr_dt[, rGDP_full_ch := (rGDP_full - rGDP_bsln) / rGDP_bsln * 100]
omr_dt[, rGDP_cndl_ch := (rGDP_cndl - rGDP_bsln) / rGDP_bsln * 100]

# Combine OMR and IMR indexes
all_indexes <- merge(omr_dt, imr_dt, by = "country", all = TRUE)
saveRDS(all_indexes, "02_output/rds/all_indexes_geppml_agri_higher_sigma.rds")

# List of African countries in AfCFTA
afcfta_countries <- c(
  "DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "COM", 
  "COG", "CIV", "COD", "DJI", "EGY", "GNQ", "SWZ", "ETH", "GAB", "GMB", "GHA",
  "GIN", "GNB", "KEN", "LSO", "LBR", "LBY", "MDG", "MWI", "MLI", "MRT", "MUS",
  "MAR", "MOZ", "NAM", "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM",
  "ZAF", "SSD", "SDN", "TZA", "TGO", "TUN", "UGA", "ZMB", "ZWE"
)

# Add a column to identify AfCFTA countries
all_indexes[, is_afcfta := country %in% afcfta_countries]

# Calculate total output for ROW for weighting
row_total_output <- all_indexes[is_afcfta == FALSE, sum(output_bsln, na.rm = TRUE)]

# Calculate weighted welfare change for ROW countries
# Following Larch, Tan, and Yotov (2023) using output shares as weights
weighted_row_welfare <- all_indexes[is_afcfta == FALSE, 
                                    sum(rGDP_full_ch * output_bsln, na.rm = TRUE) / 
                                      sum(output_bsln, na.rm = TRUE)]

# Create a copy of the dataset with just AfCFTA countries
afcfta_indexes <- all_indexes[is_afcfta == TRUE]

# Add the aggregated ROW row
row_entry <- data.table(
  country = "ROW",
  is_afcfta = FALSE,
  omr_full = mean(all_indexes[is_afcfta == FALSE, omr_full], na.rm = TRUE),
  omr_cndl = mean(all_indexes[is_afcfta == FALSE, omr_cndl], na.rm = TRUE),
  omr_bsln = mean(all_indexes[is_afcfta == FALSE, omr_bsln], na.rm = TRUE),
  rGDP_full = mean(all_indexes[is_afcfta == FALSE, rGDP_full], na.rm = TRUE),
  rGDP_cndl = mean(all_indexes[is_afcfta == FALSE, rGDP_cndl], na.rm = TRUE),
  rGDP_bsln = mean(all_indexes[is_afcfta == FALSE, rGDP_bsln], na.rm = TRUE),
  tot_exp_full = sum(all_indexes[is_afcfta == FALSE, tot_exp_full], na.rm = TRUE),
  tot_exp_cndl = sum(all_indexes[is_afcfta == FALSE, tot_exp_cndl], na.rm = TRUE),
  tot_exp_bsln = sum(all_indexes[is_afcfta == FALSE, tot_exp_bsln], na.rm = TRUE),
  p_full = mean(all_indexes[is_afcfta == FALSE, p_full], na.rm = TRUE),
  output_bsln = row_total_output,
  acr_bsln = mean(all_indexes[is_afcfta == FALSE, acr_bsln], na.rm = TRUE),
  omr_full_ch = mean(all_indexes[is_afcfta == FALSE, omr_full_ch], na.rm = TRUE),
  omr_cndl_ch = mean(all_indexes[is_afcfta == FALSE, omr_cndl_ch], na.rm = TRUE),
  rGDP_full_ch = weighted_row_welfare,
  rGDP_cndl_ch = mean(all_indexes[is_afcfta == FALSE, rGDP_cndl_ch], na.rm = TRUE),
  imr_full = mean(all_indexes[is_afcfta == FALSE, imr_full], na.rm = TRUE),
  imr_bsln = mean(all_indexes[is_afcfta == FALSE, imr_bsln], na.rm = TRUE),
  imr_cndl = mean(all_indexes[is_afcfta == FALSE, imr_cndl], na.rm = TRUE),
  imr_full_ch = mean(all_indexes[is_afcfta == FALSE, imr_full_ch], na.rm = TRUE),
  imr_cndl_ch = mean(all_indexes[is_afcfta == FALSE, imr_cndl_ch], na.rm = TRUE)
)

# Combine AfCFTA countries and ROW
all_indexes_with_row <- rbindlist(list(afcfta_indexes, row_entry), fill = TRUE)

# Save the aggregated data
saveRDS(all_indexes_with_row, "02_output/rds/all_indexes_with_row_geppml_agri_higher_sigma.rds")

# Create summary statistics for AfCFTA countries
summary_stats_afcfta <- all_indexes[is_afcfta == TRUE, .(
  mean_welfare = mean(rGDP_full_ch, na.rm = TRUE),
  median_welfare = median(rGDP_full_ch, na.rm = TRUE),
  min_welfare = min(rGDP_full_ch, na.rm = TRUE),
  max_welfare = max(rGDP_full_ch, na.rm = TRUE),
  sd_welfare = sd(rGDP_full_ch, na.rm = TRUE)
)]

# Print summary for AfCFTA countries
cat("\nSummary of Welfare Effects for AfCFTA Countries (% change in real GDP) - higher sigma:\n")
print(summary_stats_afcfta)

# Print ROW welfare effect
cat("\nWelfare Effect for Rest of World (ROW) - higher sigma: ", round(weighted_row_welfare, 4), "%\n")


################################################################################
# Section 2: General Equilibrium PPML Analysis for Manufacturing Sector --------
################################################################################


# Clear workspace
rm(list = ls())

# Load required packages
library(data.table)
library(fixest)
library(dplyr)
library(tidyr)
library(readstata13)
library(ggplot2)

# Set working directory
setwd("C:/Users/muham/AfCFTA_UNU_CRIS")

# Set parameters
sigma <- 10

################################################################################
# I. Prepare Data for Manufacturing Sector
################################################################################

# Load data
data <- read.dta13("01_input/afcfta_2019_manu_balanced.dta")
dt <- as.data.table(data)

# 1. Create aggregate variables
dt[, output := sum(trade), by = "exporter"]
dt[, expndr := sum(trade), by = "importer"]

# 2. Choose a country for reference group (Germany)
dt[exporter == "DEU", exporter := "ZZZ"]
dt[importer == "DEU", importer := "ZZZ"]
dt[importer == "ZZZ", expndr_deu0 := expndr]
dt[, expndr_deu := mean(expndr_deu0, na.rm = TRUE)]

# Save processed data
saveRDS(dt, "01_input/ge_ppml_data_manu.rds")

################################################################################
# II. GE Analysis in R
################################################################################

# Load prepared data
dt <- readRDS("01_input/ge_ppml_data_manu.rds")

#############################
# Step 1: "Baseline" Scenario
#############################

# Step 1.a: Estimate "Baseline" Gravity
# Create country dummy variables (for exporter and importer fixed effects)
countries <- unique(dt$exporter)
NoC <- length(countries)

# Baseline PPML with fixest
# Using fepois for PPML with high-dimensional fixed effects
baseline_model <- fepois(
  trade ~ cntg + col + lang + dist + rta + afcfta_brdr + afcfta_row_brdr + row_to_row | exporter + importer,
  data = dt
)

# Save the model summary
baseline_summary <- summary(baseline_model)

# Extract coefficients
coeffs <- coef(baseline_model)
CNTG_est <- coeffs["cntg"]
COL_est <- coeffs["col"]
LANG_est <- coeffs["lang"]
DIST_est <- coeffs["dist"]
RTA_est <- coeffs["rta"]
AfCFTA_BRDR_est <- coeffs["afcfta_brdr"]
AfCFTA_ROW_BRDR_est <- coeffs["afcfta_row_brdr"]
ROW_to_ROW_est <- coeffs["row_to_row"]

# Predict trade in the baseline
dt[, trade_bsln := predict(baseline_model, newdata = dt, type = "response")]

# Create baseline and counterfactual trade costs
dt[, t_ij_bsln := exp(
  DIST_est * dist +
  CNTG_est * cntg +
  COL_est * col +
  LANG_est * lang +
  RTA_est * rta +
  AfCFTA_BRDR_est * afcfta_brdr +
  AfCFTA_ROW_BRDR_est * afcfta_row_brdr +
  ROW_to_ROW_est * row_to_row
)]

# Counterfactual: set afcfta_brdr = 0, keep everything else the same
dt[, t_ij_ctrf := exp(
  DIST_est * dist +
  CNTG_est * cntg +
  COL_est * col +
  LANG_est * lang +
  RTA_est * rta +
  AfCFTA_BRDR_est * 0 +
  AfCFTA_ROW_BRDR_est * afcfta_row_brdr +
  ROW_to_ROW_est * row_to_row
)]

# Keep domestic trade costs at baseline
dt[exporter == importer, t_ij_ctrf := t_ij_bsln]
dt[, t_ij_ctrf_1 := log(t_ij_ctrf)]

# Step 1.b: Construct "Baseline" GE Indexes
# Extract fixed effects
fe_exporter <- fixef(baseline_model)$exporter
fe_importer <- fixef(baseline_model)$importer

# Convert fixed effects to data table
exp_fe_dt <- data.table(exporter = names(fe_exporter), exp_fe = exp(fe_exporter))
imp_fe_dt <- data.table(importer = names(fe_importer), imp_fe = exp(fe_importer))

# Merge fixed effects with main data
dt <- merge(dt, exp_fe_dt, by = "exporter", all.x = TRUE)
dt <- merge(dt, imp_fe_dt, by = "importer", all.x = TRUE)

# Calculate outward and inward multilateral resistance terms
dt[, all_exp_fes_0 := exp_fe]
dt[, all_imp_fes_0 := imp_fe]

# Calculate importer-specific exporter fixed effects
dt[exporter == importer, temp := all_exp_fes_0]
dt[, all_exp_fes_0_TB := mean(temp, na.rm = TRUE), by = "importer"]
dt[, temp := NULL]

# Equation (7): Outward Multilateral Resistance
dt[, omr_bsln := output * expndr_deu / all_exp_fes_0]

# Equation (8): Inward Multilateral Resistance
dt[, imr_bsln := expndr / (all_imp_fes_0 * expndr_deu)]

# Real GDP (baseline)
dt[exporter == importer, rGDP_bsln_temp := output / (imr_bsln^(1/(1-sigma)))]
dt[, rGDP_bsln := sum(rGDP_bsln_temp, na.rm = TRUE), by = "exporter"]

# Domestic absorption share (acr)
dt[exporter == importer, acr_bsln := trade_bsln / expndr]

# Bilateral exports & totals
dt[exporter != importer, exp_bsln := trade_bsln]
dt[exporter == importer, exp_bsln_acr := trade_bsln]
dt[, tot_exp_bsln := sum(exp_bsln, na.rm = TRUE), by = "exporter"]

################################
# Step 2: "Conditional" Scenario
################################

# Step 2.a: Estimate "Conditional" Gravity
# PPML with offset = log(counterfactual trade costs)
conditional_model <- fepois(
  trade ~ 1 | exporter + importer,
  data = dt,
  offset = ~t_ij_ctrf_1
)

# Predict trade under conditional GE
dt[, trade_cndl := predict(conditional_model, newdata = dt, type = "response")]

# Step 2.b: Construct "Conditional" GE Indexes
# Extract conditional fixed effects
fe_exporter_cndl <- fixef(conditional_model)$exporter
fe_importer_cndl <- fixef(conditional_model)$importer

# Convert fixed effects to data table
exp_fe_cndl_dt <- data.table(exporter = names(fe_exporter_cndl), exp_fe_cndl = exp(fe_exporter_cndl))
imp_fe_cndl_dt <- data.table(importer = names(fe_importer_cndl), imp_fe_cndl = exp(fe_importer_cndl))

# Merge fixed effects with main data
dt <- merge(dt, exp_fe_cndl_dt, by = "exporter", all.x = TRUE)
dt <- merge(dt, imp_fe_cndl_dt, by = "importer", all.x = TRUE)

# Calculate outward and inward multilateral resistance terms
dt[, all_exp_fes_1 := exp_fe_cndl]
dt[, all_imp_fes_1 := imp_fe_cndl]

# Equation (7): Outward Multilateral Resistance (conditional)
dt[, omr_cndl := output * expndr_deu / all_exp_fes_1]

# Equation (8): Inward Multilateral Resistance (conditional)
dt[, imr_cndl := expndr / (all_imp_fes_1 * expndr_deu)]

# Exports and totals (conditional)
dt[exporter != importer, exp_cndl := trade_cndl]
dt[, tot_exp_cndl := sum(exp_cndl, na.rm = TRUE), by = "exporter"]

# % change in total exports
dt[, tot_exp_cndl_ch := (tot_exp_cndl - tot_exp_bsln) / tot_exp_bsln * 100]

# Real GDP (conditional)
dt[exporter == importer, rGDP_cndl_temp := output / (imr_cndl^(1/(1-sigma)))]
dt[, rGDP_cndl := sum(rGDP_cndl_temp, na.rm = TRUE), by = "exporter"]

###################################
# Step 3: "Full Endowment" Scenario
###################################

# Step 3.a: Estimate "Full Endowment" Gravity
# Initialize variables for the iterative process
dt[, trade_1_pred := trade_cndl]
dt[, output_bsln := output]
dt[, expndr_bsln := expndr]
dt[exporter == importer, phi := expndr / output]

# Calculate initial values
dt[exporter == importer, temp := all_exp_fes_0]
dt[, all_exp_fes_0_imp := mean(temp, na.rm = TRUE), by = "importer"]
dt[, temp := NULL]

dt[exporter == importer, expndr_temp_1 := phi * output]
dt[, expndr_1 := mean(expndr_temp_1, na.rm = TRUE), by = "importer"]
dt[, expndr_temp_1 := NULL]

dt[importer == "ZZZ", expndr_deu01_1 := expndr_1]
dt[, expndr_deu_1 := mean(expndr_deu01_1, na.rm = TRUE)]

dt[exporter == importer, temp := all_exp_fes_1]
dt[, all_exp_fes_1_imp := mean(temp, na.rm = TRUE), by = "importer"]
dt[, temp := NULL]

# Initial values for factory-gate prices and multilateral resistance terms
dt[, p_full_exp_0 := 0]
dt[, p_full_exp_1 := (all_exp_fes_1 / all_exp_fes_0)^(1/(1-sigma))]
dt[, p_full_imp_1 := (all_exp_fes_1_imp / all_exp_fes_0_imp)^(1/(1-sigma))]
dt[, imr_full_1 := expndr_1 / (all_imp_fes_1 * expndr_deu_1)]
dt[, imr_full_ch_1 := 1]
dt[, omr_full_1 := output * expndr_deu_1 / all_exp_fes_1]
dt[, omr_full_ch_1 := 1]

# Set convergence parameters
max_iterations <- 30
diff_all_exp_fes_sd <- 1
diff_all_exp_fes_max <- 1
tolerance <- 0.001
i <- 3
iteration <- 0

# Define a function to check for duplicated variables and clean them up
check_and_remove_vars <- function(dt, pattern) {
  vars_to_remove <- grep(pattern, names(dt), value = TRUE)
  if (length(vars_to_remove) > 0) {
    dt[, (vars_to_remove) := NULL]
  }
  return(dt)
}

# Iterative process
cat("Starting iterative process for full endowment calculation...\n")

while ((diff_all_exp_fes_sd > tolerance || diff_all_exp_fes_max > tolerance) && iteration < max_iterations) {
  iteration <- iteration + 1
  cat(sprintf("Iteration %d of maximum %d\n", iteration, max_iterations))
  
  # Clean up variables before recreation
  current_iter <- i-1
  patterns <- c(
    paste0("trade_", current_iter),
    paste0("all_exp_fes_", current_iter),
    paste0("all_imp_fes_", current_iter),
    paste0("output_", current_iter),
    paste0("expndr_check_", current_iter),
    paste0("expndr_deu0_", current_iter),
    paste0("expndr_deu_", current_iter),
    paste0("all_exp_fes_", current_iter, "_imp"),
    paste0("p_full_exp_", current_iter),
    paste0("p_full_imp_", current_iter),
    paste0("omr_full_", current_iter),
    paste0("omr_full_ch_", current_iter),
    paste0("expndr_temp_", current_iter),
    paste0("expndr_", current_iter),
    paste0("imr_full_", current_iter),
    paste0("imr_full_ch_", current_iter),
    paste0("diff_p_full_exp_", current_iter),
    paste0("trade_", current_iter, "_pred")
  )
  
  for (pattern in patterns) {
    dt <- check_and_remove_vars(dt, pattern)
  }
  
  # Equation (14)
  dt[, paste0("trade_", current_iter) := get(paste0("trade_", i-2, "_pred")) * 
      get(paste0("p_full_exp_", i-2)) * get(paste0("p_full_imp_", i-2)) / 
      (get(paste0("omr_full_ch_", i-2)) * get(paste0("imr_full_ch_", i-2)))]
  
  # Estimate new PPML model
  tryCatch({
    formula_str <- paste0("trade_", current_iter, " ~ 1 | exporter + importer")
    full_model <- fepois(
      as.formula(formula_str),
      data = dt,
      offset = ~t_ij_ctrf_1
    )
    
    # Predict trade flows
    dt[, paste0("trade_", current_iter, "_pred") := predict(full_model, newdata = dt, type = "response")]
    
    # Extract fixed effects
    fe_exporter_full <- fixef(full_model)$exporter
    fe_importer_full <- fixef(full_model)$importer
    
    # Process fixed effects
    dt[, paste0("all_exp_fes_", current_iter) := 0]
    dt[, paste0("all_imp_fes_", current_iter) := 0]
    
    # Apply exporter fixed effects
    for (country in names(fe_exporter_full)) {
      dt[exporter == country, paste0("all_exp_fes_", current_iter) := exp(fe_exporter_full[country])]
    }
    
    # Apply importer fixed effects
    for (country in names(fe_importer_full)) {
      dt[importer == country, paste0("all_imp_fes_", current_iter) := exp(fe_importer_full[country])]
    }
    
    # Update output
    dt[, paste0("output_", current_iter) := sum(get(paste0("trade_", current_iter, "_pred"))), by = "exporter"]
    
    # Update expenditure
    dt[, paste0("expndr_check_", current_iter) := sum(get(paste0("trade_", current_iter, "_pred"))), by = "importer"]
    dt[importer == "ZZZ", paste0("expndr_deu0_", current_iter) := get(paste0("expndr_check_", current_iter))]
    dt[, paste0("expndr_deu_", current_iter) := mean(get(paste0("expndr_deu0_", current_iter)), na.rm = TRUE)]
    
    # Get importer-specific exporter fixed effects
    dt[exporter == importer, temp := get(paste0("all_exp_fes_", current_iter))]
    dt[, paste0("all_exp_fes_", current_iter, "_imp") := mean(temp, na.rm = TRUE), by = "importer"]
    dt[, temp := NULL]
    
    # Update factory-gate prices
    dt[, paste0("p_full_exp_", current_iter) := 
        ((get(paste0("all_exp_fes_", current_iter)) / get(paste0("all_exp_fes_", i-2))) / 
           (get(paste0("expndr_deu_", current_iter)) / get(paste0("expndr_deu_", i-2))))^(1/(1-sigma))]
    
    dt[, paste0("p_full_imp_", current_iter) := 
        ((get(paste0("all_exp_fes_", current_iter, "_imp")) / get(paste0("all_exp_fes_", i-2, "_imp"))) / 
           (get(paste0("expndr_deu_", current_iter)) / get(paste0("expndr_deu_", i-2))))^(1/(1-sigma))]
    
    # Equation (7) - Update outward multilateral resistance
    dt[, paste0("omr_full_", current_iter) := get(paste0("output_", current_iter)) / get(paste0("all_exp_fes_", current_iter))]
    dt[, paste0("omr_full_ch_", current_iter) := get(paste0("omr_full_", current_iter)) / get(paste0("omr_full_", i-2))]
    
    # Update expenditure
    dt[exporter == importer, paste0("expndr_temp_", current_iter) := phi * get(paste0("output_", current_iter))]
    dt[, paste0("expndr_", current_iter) := mean(get(paste0("expndr_temp_", current_iter)), na.rm = TRUE), by = "importer"]
    
    # Equation (8) - Update inward multilateral resistance
    dt[, paste0("imr_full_", current_iter) := get(paste0("expndr_", current_iter)) / 
         (get(paste0("all_imp_fes_", current_iter)) * get(paste0("expndr_deu_", current_iter)))]
    dt[, paste0("imr_full_ch_", current_iter) := get(paste0("imr_full_", current_iter)) / get(paste0("imr_full_", i-2))]
    
    # Calculate convergence criteria
    dt[, paste0("diff_p_full_exp_", current_iter) := get(paste0("p_full_exp_", i-2)) - get(paste0("p_full_exp_", i-3))]
    
    # Calculate convergence statistics
    diff_stats <- dt[, .(sd = sd(get(paste0("diff_p_full_exp_", current_iter)), na.rm = TRUE),
                          max = max(abs(get(paste0("diff_p_full_exp_", current_iter))), na.rm = TRUE))]
    diff_all_exp_fes_sd <- diff_stats$sd
    diff_all_exp_fes_max <- diff_stats$max
    
    cat(sprintf("Convergence stats - SD: %f, Max: %f\n", diff_all_exp_fes_sd, diff_all_exp_fes_max))
    
  }, error = function(e) {
    cat("Error in iteration", iteration, ":", e$message, "\n")
    cat("Using values from previous iteration\n")
    
    # Create placeholder with previous values
    dt[, paste0("trade_", current_iter, "_pred") := get(paste0("trade_", i-2, "_pred"))]
  })
  
  # Increment counter for next iteration
  i <- i + 1
}

cat(sprintf("Convergence achieved after %d iterations\n", iteration))

# Get the final iteration number for use in subsequent calculations
final_iter <- i - 2

# Step 3.b: Construct "Full Endowment" GE Indexes
# Calculate p^c/p
dt[, c("output", "expndr", "expndr_deu0") := NULL]
dt[, expndr_deu_bsln := expndr_deu]

dt[, output := sum(get(paste0("trade_", final_iter, "_pred"))), by = "exporter"]
dt[exporter == importer, expndr_temp := phi * output]
dt[, expndr := mean(expndr_temp, na.rm = TRUE), by = "importer"]
dt[importer == "ZZZ", expndr_deu0 := expndr]
dt[, expndr_deu := mean(expndr_deu0, na.rm = TRUE)]

# Calculate p_full
dt[, p_full := ((get(paste0("all_exp_fes_", final_iter)) / all_exp_fes_0) / 
                 (expndr_deu / expndr_deu_bsln))^(1/(1-sigma))]

# Calculate output_full
dt[, output_full := p_full * output_bsln]

# Equation (7): Outward Multilateral Resistance
dt[, omr_full := output_full * expndr_deu / get(paste0("all_exp_fes_", final_iter))]

# Equation (8): Inward Multilateral Resistance
dt[, imr_full := expndr / (get(paste0("all_imp_fes_", final_iter)) * expndr_deu)]

# Real GDP (full endowment)
dt[exporter == importer, rGDP_full_temp := p_full * output_bsln / (imr_full^(1/(1-sigma)))]
dt[, rGDP_full := sum(rGDP_full_temp, na.rm = TRUE), by = "exporter"]

# Expenditure (full endowment)
dt[exporter == importer, expndr_full_temp := phi * output_full]
dt[, expndr_full := mean(expndr_full_temp, na.rm = TRUE), by = "importer"]

# Bilateral trade at full endowment (using counterfactual costs)
dt[, trade_full := (output_full * expndr_full * t_ij_ctrf) / (imr_full * omr_full)]
dt[exporter != importer, exp_full := trade_full]
dt[, tot_exp_full := sum(exp_full, na.rm = TRUE), by = "exporter"]
dt[, tot_exp_full_ch := (tot_exp_full - tot_exp_bsln) / tot_exp_bsln * 100]

# Save results
saveRDS(dt, "02_output/rds/full_static_all_manu_higher_sigma.rds")

# Prepare all indexes at the country level
# IMR indexes
imr_dt <- dt[, .(imr_full = mean(imr_full, na.rm = TRUE),
                 imr_bsln = mean(imr_bsln, na.rm = TRUE),
                 imr_cndl = mean(imr_cndl, na.rm = TRUE)), by = importer]
setnames(imr_dt, "importer", "country")
imr_dt[, imr_full_ch := (imr_full - imr_bsln) / imr_bsln * 100]
imr_dt[, imr_cndl_ch := (imr_cndl - imr_bsln) / imr_bsln * 100]
saveRDS(imr_dt, "02_output/rds/imrs_all_manu_higher_sigma.rds")

# OMR and other indexes
omr_dt <- dt[, .(omr_full = mean(omr_full, na.rm = TRUE),
                 omr_cndl = mean(omr_cndl, na.rm = TRUE),
                 omr_bsln = mean(omr_bsln, na.rm = TRUE),
                 rGDP_full = mean(rGDP_full, na.rm = TRUE),
                 rGDP_cndl = mean(rGDP_cndl, na.rm = TRUE),
                 rGDP_bsln = mean(rGDP_bsln, na.rm = TRUE),
                 tot_exp_full = mean(tot_exp_full, na.rm = TRUE),
                 tot_exp_cndl = mean(tot_exp_cndl, na.rm = TRUE),
                 tot_exp_bsln = mean(tot_exp_bsln, na.rm = TRUE),
                 p_full = mean(p_full, na.rm = TRUE),
                 output_bsln = mean(output_bsln, na.rm = TRUE),
                 acr_bsln = mean(acr_bsln, na.rm = TRUE)), by = exporter]
setnames(omr_dt, "exporter", "country")
omr_dt[, omr_full_ch := (omr_full - omr_bsln) / omr_bsln * 100]
omr_dt[, omr_cndl_ch := (omr_cndl - omr_bsln) / omr_bsln * 100]
omr_dt[, rGDP_full_ch := (rGDP_full - rGDP_bsln) / rGDP_bsln * 100]
omr_dt[, rGDP_cndl_ch := (rGDP_cndl - rGDP_bsln) / rGDP_bsln * 100]

# Combine OMR and IMR indexes
all_indexes <- merge(omr_dt, imr_dt, by = "country", all = TRUE)
saveRDS(all_indexes, "02_output/rds/all_indexes_geppml_manu_higher_sigma.rds")

# List of African countries in AfCFTA
afcfta_countries <- c(
  "DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "COM", 
  "COG", "CIV", "COD", "DJI", "EGY", "GNQ", "SWZ", "ETH", "GAB", "GMB", "GHA",
  "GIN", "GNB", "KEN", "LSO", "LBR", "LBY", "MDG", "MWI", "MLI", "MRT", "MUS",
  "MAR", "MOZ", "NAM", "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM",
  "ZAF", "SSD", "SDN", "TZA", "TGO", "TUN", "UGA", "ZMB", "ZWE"
)

# Add a column to identify AfCFTA countries
all_indexes[, is_afcfta := country %in% afcfta_countries]

# Calculate total output for ROW for weighting
row_total_output <- all_indexes[is_afcfta == FALSE, sum(output_bsln, na.rm = TRUE)]

# Calculate weighted welfare change for ROW countries
# Following Larch, Tan, and Yotov (2023) using output shares as weights
weighted_row_welfare <- all_indexes[is_afcfta == FALSE, 
                                    sum(rGDP_full_ch * output_bsln, na.rm = TRUE) / 
                                      sum(output_bsln, na.rm = TRUE)]

# Create a copy of the dataset with just AfCFTA countries
afcfta_indexes <- all_indexes[is_afcfta == TRUE]

# Add the aggregated ROW row
row_entry <- data.table(
  country = "ROW",
  is_afcfta = FALSE,
  omr_full = mean(all_indexes[is_afcfta == FALSE, omr_full], na.rm = TRUE),
  omr_cndl = mean(all_indexes[is_afcfta == FALSE, omr_cndl], na.rm = TRUE),
  omr_bsln = mean(all_indexes[is_afcfta == FALSE, omr_bsln], na.rm = TRUE),
  rGDP_full = mean(all_indexes[is_afcfta == FALSE, rGDP_full], na.rm = TRUE),
  rGDP_cndl = mean(all_indexes[is_afcfta == FALSE, rGDP_cndl], na.rm = TRUE),
  rGDP_bsln = mean(all_indexes[is_afcfta == FALSE, rGDP_bsln], na.rm = TRUE),
  tot_exp_full = sum(all_indexes[is_afcfta == FALSE, tot_exp_full], na.rm = TRUE),
  tot_exp_cndl = sum(all_indexes[is_afcfta == FALSE, tot_exp_cndl], na.rm = TRUE),
  tot_exp_bsln = sum(all_indexes[is_afcfta == FALSE, tot_exp_bsln], na.rm = TRUE),
  p_full = mean(all_indexes[is_afcfta == FALSE, p_full], na.rm = TRUE),
  output_bsln = row_total_output,
  acr_bsln = mean(all_indexes[is_afcfta == FALSE, acr_bsln], na.rm = TRUE),
  omr_full_ch = mean(all_indexes[is_afcfta == FALSE, omr_full_ch], na.rm = TRUE),
  omr_cndl_ch = mean(all_indexes[is_afcfta == FALSE, omr_cndl_ch], na.rm = TRUE),
  rGDP_full_ch = weighted_row_welfare,
  rGDP_cndl_ch = mean(all_indexes[is_afcfta == FALSE, rGDP_cndl_ch], na.rm = TRUE),
  imr_full = mean(all_indexes[is_afcfta == FALSE, imr_full], na.rm = TRUE),
  imr_bsln = mean(all_indexes[is_afcfta == FALSE, imr_bsln], na.rm = TRUE),
  imr_cndl = mean(all_indexes[is_afcfta == FALSE, imr_cndl], na.rm = TRUE),
  imr_full_ch = mean(all_indexes[is_afcfta == FALSE, imr_full_ch], na.rm = TRUE),
  imr_cndl_ch = mean(all_indexes[is_afcfta == FALSE, imr_cndl_ch], na.rm = TRUE)
)

# Combine AfCFTA countries and ROW
all_indexes_with_row <- rbindlist(list(afcfta_indexes, row_entry), fill = TRUE)

# Save the aggregated data
saveRDS(all_indexes_with_row, "02_output/rds/all_indexes_with_row_geppml_manu_higher_sigma.rds")

# Create summary statistics for AfCFTA countries
summary_stats_afcfta <- all_indexes[is_afcfta == TRUE, .(
  mean_welfare = mean(rGDP_full_ch, na.rm = TRUE),
  median_welfare = median(rGDP_full_ch, na.rm = TRUE),
  min_welfare = min(rGDP_full_ch, na.rm = TRUE),
  max_welfare = max(rGDP_full_ch, na.rm = TRUE),
  sd_welfare = sd(rGDP_full_ch, na.rm = TRUE)
)]

# Print summary for AfCFTA countries
cat("\nSummary of Welfare Effects for AfCFTA Countries (% change in real GDP) - higher sigma:\n")
print(summary_stats_afcfta)

# Print ROW welfare effect
cat("\nWelfare Effect for Rest of World (ROW) - higher sigma: ", round(weighted_row_welfare, 4), "%\n")


################################################################################
# Section 3: General Equilibrium PPML Analysis for Mining and Energy Sector ----
################################################################################

# Clear workspace
rm(list = ls())

# Load required packages
library(data.table)
library(fixest)
library(dplyr)
library(tidyr)
library(readstata13)
library(ggplot2)

# Set working directory
setwd("C:/Users/muham/AfCFTA_UNU_CRIS")

# Set parameters
sigma <- 10

################################################################################
# I. Prepare Data for Mining and Energy Sector
################################################################################

# Load data
data <- read.dta13("01_input/afcfta_2019_mine_balanced.dta")
dt <- as.data.table(data)

# 1. Create aggregate variables
dt[, output := sum(trade), by = "exporter"]
dt[, expndr := sum(trade), by = "importer"]

# 2. Choose a country for reference group (Germany)
dt[exporter == "DEU", exporter := "ZZZ"]
dt[importer == "DEU", importer := "ZZZ"]
dt[importer == "ZZZ", expndr_deu0 := expndr]
dt[, expndr_deu := mean(expndr_deu0, na.rm = TRUE)]

# Save processed data
saveRDS(dt, "01_input/ge_ppml_data_mine.rds")

################################################################################
# II. GE Analysis in R
################################################################################

# Load prepared data
dt <- readRDS("01_input/ge_ppml_data_mine.rds")

#############################
# Step 1: "Baseline" Scenario
#############################

# Step 1.a: Estimate "Baseline" Gravity
# Create country dummy variables (for exporter and importer fixed effects)
countries <- unique(dt$exporter)
NoC <- length(countries)

# Baseline PPML with fixest
# Using fepois for PPML with high-dimensional fixed effects
baseline_model <- fepois(
  trade ~ cntg + col + lang + dist + rta + afcfta_brdr + afcfta_row_brdr + row_to_row | exporter + importer,
  data = dt
)

# Save the model summary
baseline_summary <- summary(baseline_model)

# Extract coefficients
coeffs <- coef(baseline_model)
CNTG_est <- coeffs["cntg"]
COL_est <- coeffs["col"]
LANG_est <- coeffs["lang"]
DIST_est <- coeffs["dist"]
RTA_est <- coeffs["rta"]
AfCFTA_BRDR_est <- coeffs["afcfta_brdr"]
AfCFTA_ROW_BRDR_est <- coeffs["afcfta_row_brdr"]
ROW_to_ROW_est <- coeffs["row_to_row"]

# Predict trade in the baseline
dt[, trade_bsln := predict(baseline_model, newdata = dt, type = "response")]

# Create baseline and counterfactual trade costs
dt[, t_ij_bsln := exp(
  DIST_est * dist +
  CNTG_est * cntg +
  COL_est * col +
  LANG_est * lang +
  RTA_est * rta +
  AfCFTA_BRDR_est * afcfta_brdr +
  AfCFTA_ROW_BRDR_est * afcfta_row_brdr +
  ROW_to_ROW_est * row_to_row
)]

# Counterfactual: set afcfta_brdr = 0, keep everything else the same
dt[, t_ij_ctrf := exp(
  DIST_est * dist +
  CNTG_est * cntg +
  COL_est * col +
  LANG_est * lang +
  RTA_est * rta +
  AfCFTA_BRDR_est * 0 +
  AfCFTA_ROW_BRDR_est * afcfta_row_brdr +
  ROW_to_ROW_est * row_to_row
)]

# Keep domestic trade costs at baseline
dt[exporter == importer, t_ij_ctrf := t_ij_bsln]
dt[, t_ij_ctrf_1 := log(t_ij_ctrf)]

# Step 1.b: Construct "Baseline" GE Indexes
# Extract fixed effects
fe_exporter <- fixef(baseline_model)$exporter
fe_importer <- fixef(baseline_model)$importer

# Convert fixed effects to data table
exp_fe_dt <- data.table(exporter = names(fe_exporter), exp_fe = exp(fe_exporter))
imp_fe_dt <- data.table(importer = names(fe_importer), imp_fe = exp(fe_importer))

# Merge fixed effects with main data
dt <- merge(dt, exp_fe_dt, by = "exporter", all.x = TRUE)
dt <- merge(dt, imp_fe_dt, by = "importer", all.x = TRUE)

# Calculate outward and inward multilateral resistance terms
dt[, all_exp_fes_0 := exp_fe]
dt[, all_imp_fes_0 := imp_fe]

# Calculate importer-specific exporter fixed effects
dt[exporter == importer, temp := all_exp_fes_0]
dt[, all_exp_fes_0_TB := mean(temp, na.rm = TRUE), by = "importer"]
dt[, temp := NULL]

# Equation (7): Outward Multilateral Resistance
dt[, omr_bsln := output * expndr_deu / all_exp_fes_0]

# Equation (8): Inward Multilateral Resistance
dt[, imr_bsln := expndr / (all_imp_fes_0 * expndr_deu)]

# Real GDP (baseline)
dt[exporter == importer, rGDP_bsln_temp := output / (imr_bsln^(1/(1-sigma)))]
dt[, rGDP_bsln := sum(rGDP_bsln_temp, na.rm = TRUE), by = "exporter"]

# Domestic absorption share (acr)
dt[exporter == importer, acr_bsln := trade_bsln / expndr]

# Bilateral exports & totals
dt[exporter != importer, exp_bsln := trade_bsln]
dt[exporter == importer, exp_bsln_acr := trade_bsln]
dt[, tot_exp_bsln := sum(exp_bsln, na.rm = TRUE), by = "exporter"]

################################
# Step 2: "Conditional" Scenario
################################

# Step 2.a: Estimate "Conditional" Gravity
# PPML with offset = log(counterfactual trade costs)
conditional_model <- fepois(
  trade ~ 1 | exporter + importer,
  data = dt,
  offset = ~t_ij_ctrf_1
)

# Predict trade under conditional GE
dt[, trade_cndl := predict(conditional_model, newdata = dt, type = "response")]

# Step 2.b: Construct "Conditional" GE Indexes
# Extract conditional fixed effects
fe_exporter_cndl <- fixef(conditional_model)$exporter
fe_importer_cndl <- fixef(conditional_model)$importer

# Convert fixed effects to data table
exp_fe_cndl_dt <- data.table(exporter = names(fe_exporter_cndl), exp_fe_cndl = exp(fe_exporter_cndl))
imp_fe_cndl_dt <- data.table(importer = names(fe_importer_cndl), imp_fe_cndl = exp(fe_importer_cndl))

# Merge fixed effects with main data
dt <- merge(dt, exp_fe_cndl_dt, by = "exporter", all.x = TRUE)
dt <- merge(dt, imp_fe_cndl_dt, by = "importer", all.x = TRUE)

# Calculate outward and inward multilateral resistance terms
dt[, all_exp_fes_1 := exp_fe_cndl]
dt[, all_imp_fes_1 := imp_fe_cndl]

# Equation (7): Outward Multilateral Resistance (conditional)
dt[, omr_cndl := output * expndr_deu / all_exp_fes_1]

# Equation (8): Inward Multilateral Resistance (conditional)
dt[, imr_cndl := expndr / (all_imp_fes_1 * expndr_deu)]

# Exports and totals (conditional)
dt[exporter != importer, exp_cndl := trade_cndl]
dt[, tot_exp_cndl := sum(exp_cndl, na.rm = TRUE), by = "exporter"]

# % change in total exports
dt[, tot_exp_cndl_ch := (tot_exp_cndl - tot_exp_bsln) / tot_exp_bsln * 100]

# Real GDP (conditional)
dt[exporter == importer, rGDP_cndl_temp := output / (imr_cndl^(1/(1-sigma)))]
dt[, rGDP_cndl := sum(rGDP_cndl_temp, na.rm = TRUE), by = "exporter"]

###################################
# Step 3: "Full Endowment" Scenario
###################################

# Step 3.a: Estimate "Full Endowment" Gravity
# Initialize variables for the iterative process
dt[, trade_1_pred := trade_cndl]
dt[, output_bsln := output]
dt[, expndr_bsln := expndr]
dt[exporter == importer, phi := expndr / output]

# Calculate initial values
dt[exporter == importer, temp := all_exp_fes_0]
dt[, all_exp_fes_0_imp := mean(temp, na.rm = TRUE), by = "importer"]
dt[, temp := NULL]

dt[exporter == importer, expndr_temp_1 := phi * output]
dt[, expndr_1 := mean(expndr_temp_1, na.rm = TRUE), by = "importer"]
dt[, expndr_temp_1 := NULL]

dt[importer == "ZZZ", expndr_deu01_1 := expndr_1]
dt[, expndr_deu_1 := mean(expndr_deu01_1, na.rm = TRUE)]

dt[exporter == importer, temp := all_exp_fes_1]
dt[, all_exp_fes_1_imp := mean(temp, na.rm = TRUE), by = "importer"]
dt[, temp := NULL]

# Initial values for factory-gate prices and multilateral resistance terms
dt[, p_full_exp_0 := 0]
dt[, p_full_exp_1 := (all_exp_fes_1 / all_exp_fes_0)^(1/(1-sigma))]
dt[, p_full_imp_1 := (all_exp_fes_1_imp / all_exp_fes_0_imp)^(1/(1-sigma))]
dt[, imr_full_1 := expndr_1 / (all_imp_fes_1 * expndr_deu_1)]
dt[, imr_full_ch_1 := 1]
dt[, omr_full_1 := output * expndr_deu_1 / all_exp_fes_1]
dt[, omr_full_ch_1 := 1]

# Set convergence parameters
max_iterations <- 30
diff_all_exp_fes_sd <- 1
diff_all_exp_fes_max <- 1
tolerance <- 0.001
i <- 3
iteration <- 0

# Define a function to check for duplicated variables and clean them up
check_and_remove_vars <- function(dt, pattern) {
  vars_to_remove <- grep(pattern, names(dt), value = TRUE)
  if (length(vars_to_remove) > 0) {
    dt[, (vars_to_remove) := NULL]
  }
  return(dt)
}

# Iterative process
cat("Starting iterative process for full endowment calculation...\n")

while ((diff_all_exp_fes_sd > tolerance || diff_all_exp_fes_max > tolerance) && iteration < max_iterations) {
  iteration <- iteration + 1
  cat(sprintf("Iteration %d of maximum %d\n", iteration, max_iterations))
  
  # Clean up variables before recreation
  current_iter <- i-1
  patterns <- c(
    paste0("trade_", current_iter),
    paste0("all_exp_fes_", current_iter),
    paste0("all_imp_fes_", current_iter),
    paste0("output_", current_iter),
    paste0("expndr_check_", current_iter),
    paste0("expndr_deu0_", current_iter),
    paste0("expndr_deu_", current_iter),
    paste0("all_exp_fes_", current_iter, "_imp"),
    paste0("p_full_exp_", current_iter),
    paste0("p_full_imp_", current_iter),
    paste0("omr_full_", current_iter),
    paste0("omr_full_ch_", current_iter),
    paste0("expndr_temp_", current_iter),
    paste0("expndr_", current_iter),
    paste0("imr_full_", current_iter),
    paste0("imr_full_ch_", current_iter),
    paste0("diff_p_full_exp_", current_iter),
    paste0("trade_", current_iter, "_pred")
  )
  
  for (pattern in patterns) {
    dt <- check_and_remove_vars(dt, pattern)
  }
  
  # Equation (14)
  dt[, paste0("trade_", current_iter) := get(paste0("trade_", i-2, "_pred")) * 
      get(paste0("p_full_exp_", i-2)) * get(paste0("p_full_imp_", i-2)) / 
      (get(paste0("omr_full_ch_", i-2)) * get(paste0("imr_full_ch_", i-2)))]
  
  # Estimate new PPML model
  tryCatch({
    formula_str <- paste0("trade_", current_iter, " ~ 1 | exporter + importer")
    full_model <- fepois(
      as.formula(formula_str),
      data = dt,
      offset = ~t_ij_ctrf_1
    )
    
    # Predict trade flows
    dt[, paste0("trade_", current_iter, "_pred") := predict(full_model, newdata = dt, type = "response")]
    
    # Extract fixed effects
    fe_exporter_full <- fixef(full_model)$exporter
    fe_importer_full <- fixef(full_model)$importer
    
    # Process fixed effects
    dt[, paste0("all_exp_fes_", current_iter) := 0]
    dt[, paste0("all_imp_fes_", current_iter) := 0]
    
    # Apply exporter fixed effects
    for (country in names(fe_exporter_full)) {
      dt[exporter == country, paste0("all_exp_fes_", current_iter) := exp(fe_exporter_full[country])]
    }
    
    # Apply importer fixed effects
    for (country in names(fe_importer_full)) {
      dt[importer == country, paste0("all_imp_fes_", current_iter) := exp(fe_importer_full[country])]
    }
    
    # Update output
    dt[, paste0("output_", current_iter) := sum(get(paste0("trade_", current_iter, "_pred"))), by = "exporter"]
    
    # Update expenditure
    dt[, paste0("expndr_check_", current_iter) := sum(get(paste0("trade_", current_iter, "_pred"))), by = "importer"]
    dt[importer == "ZZZ", paste0("expndr_deu0_", current_iter) := get(paste0("expndr_check_", current_iter))]
    dt[, paste0("expndr_deu_", current_iter) := mean(get(paste0("expndr_deu0_", current_iter)), na.rm = TRUE)]
    
    # Get importer-specific exporter fixed effects
    dt[exporter == importer, temp := get(paste0("all_exp_fes_", current_iter))]
    dt[, paste0("all_exp_fes_", current_iter, "_imp") := mean(temp, na.rm = TRUE), by = "importer"]
    dt[, temp := NULL]
    
    # Update factory-gate prices
    dt[, paste0("p_full_exp_", current_iter) := 
        ((get(paste0("all_exp_fes_", current_iter)) / get(paste0("all_exp_fes_", i-2))) / 
           (get(paste0("expndr_deu_", current_iter)) / get(paste0("expndr_deu_", i-2))))^(1/(1-sigma))]
    
    dt[, paste0("p_full_imp_", current_iter) := 
        ((get(paste0("all_exp_fes_", current_iter, "_imp")) / get(paste0("all_exp_fes_", i-2, "_imp"))) / 
           (get(paste0("expndr_deu_", current_iter)) / get(paste0("expndr_deu_", i-2))))^(1/(1-sigma))]
    
    # Equation (7) - Update outward multilateral resistance
    dt[, paste0("omr_full_", current_iter) := get(paste0("output_", current_iter)) / get(paste0("all_exp_fes_", current_iter))]
    dt[, paste0("omr_full_ch_", current_iter) := get(paste0("omr_full_", current_iter)) / get(paste0("omr_full_", i-2))]
    
    # Update expenditure
    dt[exporter == importer, paste0("expndr_temp_", current_iter) := phi * get(paste0("output_", current_iter))]
    dt[, paste0("expndr_", current_iter) := mean(get(paste0("expndr_temp_", current_iter)), na.rm = TRUE), by = "importer"]
    
    # Equation (8) - Update inward multilateral resistance
    dt[, paste0("imr_full_", current_iter) := get(paste0("expndr_", current_iter)) / 
         (get(paste0("all_imp_fes_", current_iter)) * get(paste0("expndr_deu_", current_iter)))]
    dt[, paste0("imr_full_ch_", current_iter) := get(paste0("imr_full_", current_iter)) / get(paste0("imr_full_", i-2))]
    
    # Calculate convergence criteria
    dt[, paste0("diff_p_full_exp_", current_iter) := get(paste0("p_full_exp_", i-2)) - get(paste0("p_full_exp_", i-3))]
    
    # Calculate convergence statistics
    diff_stats <- dt[, .(sd = sd(get(paste0("diff_p_full_exp_", current_iter)), na.rm = TRUE),
                          max = max(abs(get(paste0("diff_p_full_exp_", current_iter))), na.rm = TRUE))]
    diff_all_exp_fes_sd <- diff_stats$sd
    diff_all_exp_fes_max <- diff_stats$max
    
    cat(sprintf("Convergence stats - SD: %f, Max: %f\n", diff_all_exp_fes_sd, diff_all_exp_fes_max))
    
  }, error = function(e) {
    cat("Error in iteration", iteration, ":", e$message, "\n")
    cat("Using values from previous iteration\n")
    
    # Create placeholder with previous values
    dt[, paste0("trade_", current_iter, "_pred") := get(paste0("trade_", i-2, "_pred"))]
  })
  
  # Increment counter for next iteration
  i <- i + 1
}

cat(sprintf("Convergence achieved after %d iterations\n", iteration))

# Get the final iteration number for use in subsequent calculations
final_iter <- i - 2

# Step 3.b: Construct "Full Endowment" GE Indexes
# Calculate p^c/p
dt[, c("output", "expndr", "expndr_deu0") := NULL]
dt[, expndr_deu_bsln := expndr_deu]

dt[, output := sum(get(paste0("trade_", final_iter, "_pred"))), by = "exporter"]
dt[exporter == importer, expndr_temp := phi * output]
dt[, expndr := mean(expndr_temp, na.rm = TRUE), by = "importer"]
dt[importer == "ZZZ", expndr_deu0 := expndr]
dt[, expndr_deu := mean(expndr_deu0, na.rm = TRUE)]

# Calculate p_full
dt[, p_full := ((get(paste0("all_exp_fes_", final_iter)) / all_exp_fes_0) / 
                 (expndr_deu / expndr_deu_bsln))^(1/(1-sigma))]

# Calculate output_full
dt[, output_full := p_full * output_bsln]

# Equation (7): Outward Multilateral Resistance
dt[, omr_full := output_full * expndr_deu / get(paste0("all_exp_fes_", final_iter))]

# Equation (8): Inward Multilateral Resistance
dt[, imr_full := expndr / (get(paste0("all_imp_fes_", final_iter)) * expndr_deu)]

# Real GDP (full endowment)
dt[exporter == importer, rGDP_full_temp := p_full * output_bsln / (imr_full^(1/(1-sigma)))]
dt[, rGDP_full := sum(rGDP_full_temp, na.rm = TRUE), by = "exporter"]

# Expenditure (full endowment)
dt[exporter == importer, expndr_full_temp := phi * output_full]
dt[, expndr_full := mean(expndr_full_temp, na.rm = TRUE), by = "importer"]

# Bilateral trade at full endowment (using counterfactual costs)
dt[, trade_full := (output_full * expndr_full * t_ij_ctrf) / (imr_full * omr_full)]
dt[exporter != importer, exp_full := trade_full]
dt[, tot_exp_full := sum(exp_full, na.rm = TRUE), by = "exporter"]
dt[, tot_exp_full_ch := (tot_exp_full - tot_exp_bsln) / tot_exp_bsln * 100]

# Save results
saveRDS(dt, "02_output/rds/full_static_all_mine_higher_sigma.rds")

# Prepare all indexes at the country level
# IMR indexes
imr_dt <- dt[, .(imr_full = mean(imr_full, na.rm = TRUE),
                 imr_bsln = mean(imr_bsln, na.rm = TRUE),
                 imr_cndl = mean(imr_cndl, na.rm = TRUE)), by = importer]
setnames(imr_dt, "importer", "country")
imr_dt[, imr_full_ch := (imr_full - imr_bsln) / imr_bsln * 100]
imr_dt[, imr_cndl_ch := (imr_cndl - imr_bsln) / imr_bsln * 100]
saveRDS(imr_dt, "02_output/rds/imrs_all_mine_higher_sigma.rds")

# OMR and other indexes
omr_dt <- dt[, .(omr_full = mean(omr_full, na.rm = TRUE),
                 omr_cndl = mean(omr_cndl, na.rm = TRUE),
                 omr_bsln = mean(omr_bsln, na.rm = TRUE),
                 rGDP_full = mean(rGDP_full, na.rm = TRUE),
                 rGDP_cndl = mean(rGDP_cndl, na.rm = TRUE),
                 rGDP_bsln = mean(rGDP_bsln, na.rm = TRUE),
                 tot_exp_full = mean(tot_exp_full, na.rm = TRUE),
                 tot_exp_cndl = mean(tot_exp_cndl, na.rm = TRUE),
                 tot_exp_bsln = mean(tot_exp_bsln, na.rm = TRUE),
                 p_full = mean(p_full, na.rm = TRUE),
                 output_bsln = mean(output_bsln, na.rm = TRUE),
                 acr_bsln = mean(acr_bsln, na.rm = TRUE)), by = exporter]
setnames(omr_dt, "exporter", "country")
omr_dt[, omr_full_ch := (omr_full - omr_bsln) / omr_bsln * 100]
omr_dt[, omr_cndl_ch := (omr_cndl - omr_bsln) / omr_bsln * 100]
omr_dt[, rGDP_full_ch := (rGDP_full - rGDP_bsln) / rGDP_bsln * 100]
omr_dt[, rGDP_cndl_ch := (rGDP_cndl - rGDP_bsln) / rGDP_bsln * 100]

# Combine OMR and IMR indexes
all_indexes <- merge(omr_dt, imr_dt, by = "country", all = TRUE)
saveRDS(all_indexes, "02_output/rds/all_indexes_geppml_mine_higher_sigma.rds")

# List of African countries in AfCFTA
afcfta_countries <- c(
  "DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "COM", 
  "COG", "CIV", "COD", "DJI", "EGY", "GNQ", "SWZ", "ETH", "GAB", "GMB", "GHA",
  "GIN", "GNB", "KEN", "LSO", "LBR", "LBY", "MDG", "MWI", "MLI", "MRT", "MUS",
  "MAR", "MOZ", "NAM", "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM",
  "ZAF", "SSD", "SDN", "TZA", "TGO", "TUN", "UGA", "ZMB", "ZWE"
)

# Add a column to identify AfCFTA countries
all_indexes[, is_afcfta := country %in% afcfta_countries]

# Calculate total output for ROW for weighting
row_total_output <- all_indexes[is_afcfta == FALSE, sum(output_bsln, na.rm = TRUE)]

# Calculate weighted welfare change for ROW countries
# Following Larch, Tan, and Yotov (2023) using output shares as weights
weighted_row_welfare <- all_indexes[is_afcfta == FALSE, 
                                    sum(rGDP_full_ch * output_bsln, na.rm = TRUE) / 
                                      sum(output_bsln, na.rm = TRUE)]

# Create a copy of the dataset with just AfCFTA countries
afcfta_indexes <- all_indexes[is_afcfta == TRUE]

# Add the aggregated ROW row
row_entry <- data.table(
  country = "ROW",
  is_afcfta = FALSE,
  omr_full = mean(all_indexes[is_afcfta == FALSE, omr_full], na.rm = TRUE),
  omr_cndl = mean(all_indexes[is_afcfta == FALSE, omr_cndl], na.rm = TRUE),
  omr_bsln = mean(all_indexes[is_afcfta == FALSE, omr_bsln], na.rm = TRUE),
  rGDP_full = mean(all_indexes[is_afcfta == FALSE, rGDP_full], na.rm = TRUE),
  rGDP_cndl = mean(all_indexes[is_afcfta == FALSE, rGDP_cndl], na.rm = TRUE),
  rGDP_bsln = mean(all_indexes[is_afcfta == FALSE, rGDP_bsln], na.rm = TRUE),
  tot_exp_full = sum(all_indexes[is_afcfta == FALSE, tot_exp_full], na.rm = TRUE),
  tot_exp_cndl = sum(all_indexes[is_afcfta == FALSE, tot_exp_cndl], na.rm = TRUE),
  tot_exp_bsln = sum(all_indexes[is_afcfta == FALSE, tot_exp_bsln], na.rm = TRUE),
  p_full = mean(all_indexes[is_afcfta == FALSE, p_full], na.rm = TRUE),
  output_bsln = row_total_output,
  acr_bsln = mean(all_indexes[is_afcfta == FALSE, acr_bsln], na.rm = TRUE),
  omr_full_ch = mean(all_indexes[is_afcfta == FALSE, omr_full_ch], na.rm = TRUE),
  omr_cndl_ch = mean(all_indexes[is_afcfta == FALSE, omr_cndl_ch], na.rm = TRUE),
  rGDP_full_ch = weighted_row_welfare,
  rGDP_cndl_ch = mean(all_indexes[is_afcfta == FALSE, rGDP_cndl_ch], na.rm = TRUE),
  imr_full = mean(all_indexes[is_afcfta == FALSE, imr_full], na.rm = TRUE),
  imr_bsln = mean(all_indexes[is_afcfta == FALSE, imr_bsln], na.rm = TRUE),
  imr_cndl = mean(all_indexes[is_afcfta == FALSE, imr_cndl], na.rm = TRUE),
  imr_full_ch = mean(all_indexes[is_afcfta == FALSE, imr_full_ch], na.rm = TRUE),
  imr_cndl_ch = mean(all_indexes[is_afcfta == FALSE, imr_cndl_ch], na.rm = TRUE)
)

# Combine AfCFTA countries and ROW
all_indexes_with_row <- rbindlist(list(afcfta_indexes, row_entry), fill = TRUE)

# Save the aggregated data
saveRDS(all_indexes_with_row, "02_output/rds/all_indexes_with_row_geppml_mine_higher_sigma.rds")

# Create summary statistics for AfCFTA countries
summary_stats_afcfta <- all_indexes[is_afcfta == TRUE, .(
  mean_welfare = mean(rGDP_full_ch, na.rm = TRUE),
  median_welfare = median(rGDP_full_ch, na.rm = TRUE),
  min_welfare = min(rGDP_full_ch, na.rm = TRUE),
  max_welfare = max(rGDP_full_ch, na.rm = TRUE),
  sd_welfare = sd(rGDP_full_ch, na.rm = TRUE)
)]

# Print summary for AfCFTA countries
cat("\nSummary of Welfare Effects for AfCFTA Countries (% change in real GDP) - higher sigma:\n")
print(summary_stats_afcfta)

# Print ROW welfare effect
cat("\nWelfare Effect for Rest of World (ROW) - higher sigma: ", round(weighted_row_welfare, 4), "%\n")


################################################################################
# Section 4: General Equilibrium PPML Analysis for Services Sector -------------
################################################################################


# Clear workspace
rm(list = ls())

# Load required packages
library(data.table)
library(fixest)
library(dplyr)
library(tidyr)
library(readstata13)
library(ggplot2)

# Set working directory
setwd("C:/Users/muham/AfCFTA_UNU_CRIS")

# Set parameters
sigma <- 10

################################################################################
# I. Prepare Data for Services Sector
################################################################################

# Load data
data <- read.dta13("01_input/afcfta_2019_serv_balanced.dta")
dt <- as.data.table(data)

# 1. Create aggregate variables
dt[, output := sum(trade), by = "exporter"]
dt[, expndr := sum(trade), by = "importer"]

# 2. Choose a country for reference group (Germany)
dt[exporter == "DEU", exporter := "ZZZ"]
dt[importer == "DEU", importer := "ZZZ"]
dt[importer == "ZZZ", expndr_deu0 := expndr]
dt[, expndr_deu := mean(expndr_deu0, na.rm = TRUE)]

# Save processed data
saveRDS(dt, "01_input/ge_ppml_data_serv.rds")

################################################################################
# II. GE Analysis in R
################################################################################

# Load prepared data
dt <- readRDS("01_input/ge_ppml_data_serv.rds")

#############################
# Step 1: "Baseline" Scenario
#############################

# Step 1.a: Estimate "Baseline" Gravity
# Create country dummy variables (for exporter and importer fixed effects)
countries <- unique(dt$exporter)
NoC <- length(countries)

# Baseline PPML with fixest
# Using fepois for PPML with high-dimensional fixed effects
baseline_model <- fepois(
  trade ~ cntg + col + lang + dist + rta + afcfta_brdr + afcfta_row_brdr + row_to_row | exporter + importer,
  data = dt
)

# Save the model summary
baseline_summary <- summary(baseline_model)

# Extract coefficients
coeffs <- coef(baseline_model)
CNTG_est <- coeffs["cntg"]
COL_est <- coeffs["col"]
LANG_est <- coeffs["lang"]
DIST_est <- coeffs["dist"]
RTA_est <- coeffs["rta"]
AfCFTA_BRDR_est <- coeffs["afcfta_brdr"]
AfCFTA_ROW_BRDR_est <- coeffs["afcfta_row_brdr"]
ROW_to_ROW_est <- coeffs["row_to_row"]

# Predict trade in the baseline
dt[, trade_bsln := predict(baseline_model, newdata = dt, type = "response")]

# Create baseline and counterfactual trade costs
dt[, t_ij_bsln := exp(
  DIST_est * dist +
  CNTG_est * cntg +
  COL_est * col +
  LANG_est * lang +
  RTA_est * rta +
  AfCFTA_BRDR_est * afcfta_brdr +
  AfCFTA_ROW_BRDR_est * afcfta_row_brdr +
  ROW_to_ROW_est * row_to_row
)]

# Counterfactual: set afcfta_brdr = 0, keep everything else the same
dt[, t_ij_ctrf := exp(
  DIST_est * dist +
  CNTG_est * cntg +
  COL_est * col +
  LANG_est * lang +
  RTA_est * rta +
  AfCFTA_BRDR_est * 0 +
  AfCFTA_ROW_BRDR_est * afcfta_row_brdr +
  ROW_to_ROW_est * row_to_row
)]

# Keep domestic trade costs at baseline
dt[exporter == importer, t_ij_ctrf := t_ij_bsln]
dt[, t_ij_ctrf_1 := log(t_ij_ctrf)]

# Step 1.b: Construct "Baseline" GE Indexes
# Extract fixed effects
fe_exporter <- fixef(baseline_model)$exporter
fe_importer <- fixef(baseline_model)$importer

# Convert fixed effects to data table
exp_fe_dt <- data.table(exporter = names(fe_exporter), exp_fe = exp(fe_exporter))
imp_fe_dt <- data.table(importer = names(fe_importer), imp_fe = exp(fe_importer))

# Merge fixed effects with main data
dt <- merge(dt, exp_fe_dt, by = "exporter", all.x = TRUE)
dt <- merge(dt, imp_fe_dt, by = "importer", all.x = TRUE)

# Calculate outward and inward multilateral resistance terms
dt[, all_exp_fes_0 := exp_fe]
dt[, all_imp_fes_0 := imp_fe]

# Calculate importer-specific exporter fixed effects
dt[exporter == importer, temp := all_exp_fes_0]
dt[, all_exp_fes_0_TB := mean(temp, na.rm = TRUE), by = "importer"]
dt[, temp := NULL]

# Equation (7): Outward Multilateral Resistance
dt[, omr_bsln := output * expndr_deu / all_exp_fes_0]

# Equation (8): Inward Multilateral Resistance
dt[, imr_bsln := expndr / (all_imp_fes_0 * expndr_deu)]

# Real GDP (baseline)
dt[exporter == importer, rGDP_bsln_temp := output / (imr_bsln^(1/(1-sigma)))]
dt[, rGDP_bsln := sum(rGDP_bsln_temp, na.rm = TRUE), by = "exporter"]

# Domestic absorption share (acr)
dt[exporter == importer, acr_bsln := trade_bsln / expndr]

# Bilateral exports & totals
dt[exporter != importer, exp_bsln := trade_bsln]
dt[exporter == importer, exp_bsln_acr := trade_bsln]
dt[, tot_exp_bsln := sum(exp_bsln, na.rm = TRUE), by = "exporter"]

################################
# Step 2: "Conditional" Scenario
################################

# Step 2.a: Estimate "Conditional" Gravity
# PPML with offset = log(counterfactual trade costs)
conditional_model <- fepois(
  trade ~ 1 | exporter + importer,
  data = dt,
  offset = ~t_ij_ctrf_1
)

# Predict trade under conditional GE
dt[, trade_cndl := predict(conditional_model, newdata = dt, type = "response")]

# Step 2.b: Construct "Conditional" GE Indexes
# Extract conditional fixed effects
fe_exporter_cndl <- fixef(conditional_model)$exporter
fe_importer_cndl <- fixef(conditional_model)$importer

# Convert fixed effects to data table
exp_fe_cndl_dt <- data.table(exporter = names(fe_exporter_cndl), exp_fe_cndl = exp(fe_exporter_cndl))
imp_fe_cndl_dt <- data.table(importer = names(fe_importer_cndl), imp_fe_cndl = exp(fe_importer_cndl))

# Merge fixed effects with main data
dt <- merge(dt, exp_fe_cndl_dt, by = "exporter", all.x = TRUE)
dt <- merge(dt, imp_fe_cndl_dt, by = "importer", all.x = TRUE)

# Calculate outward and inward multilateral resistance terms
dt[, all_exp_fes_1 := exp_fe_cndl]
dt[, all_imp_fes_1 := imp_fe_cndl]

# Equation (7): Outward Multilateral Resistance (conditional)
dt[, omr_cndl := output * expndr_deu / all_exp_fes_1]

# Equation (8): Inward Multilateral Resistance (conditional)
dt[, imr_cndl := expndr / (all_imp_fes_1 * expndr_deu)]

# Exports and totals (conditional)
dt[exporter != importer, exp_cndl := trade_cndl]
dt[, tot_exp_cndl := sum(exp_cndl, na.rm = TRUE), by = "exporter"]

# % change in total exports
dt[, tot_exp_cndl_ch := (tot_exp_cndl - tot_exp_bsln) / tot_exp_bsln * 100]

# Real GDP (conditional)
dt[exporter == importer, rGDP_cndl_temp := output / (imr_cndl^(1/(1-sigma)))]
dt[, rGDP_cndl := sum(rGDP_cndl_temp, na.rm = TRUE), by = "exporter"]

###################################
# Step 3: "Full Endowment" Scenario
###################################

# Step 3.a: Estimate "Full Endowment" Gravity
# Initialize variables for the iterative process
dt[, trade_1_pred := trade_cndl]
dt[, output_bsln := output]
dt[, expndr_bsln := expndr]
dt[exporter == importer, phi := expndr / output]

# Calculate initial values
dt[exporter == importer, temp := all_exp_fes_0]
dt[, all_exp_fes_0_imp := mean(temp, na.rm = TRUE), by = "importer"]
dt[, temp := NULL]

dt[exporter == importer, expndr_temp_1 := phi * output]
dt[, expndr_1 := mean(expndr_temp_1, na.rm = TRUE), by = "importer"]
dt[, expndr_temp_1 := NULL]

dt[importer == "ZZZ", expndr_deu01_1 := expndr_1]
dt[, expndr_deu_1 := mean(expndr_deu01_1, na.rm = TRUE)]

dt[exporter == importer, temp := all_exp_fes_1]
dt[, all_exp_fes_1_imp := mean(temp, na.rm = TRUE), by = "importer"]
dt[, temp := NULL]

# Initial values for factory-gate prices and multilateral resistance terms
dt[, p_full_exp_0 := 0]
dt[, p_full_exp_1 := (all_exp_fes_1 / all_exp_fes_0)^(1/(1-sigma))]
dt[, p_full_imp_1 := (all_exp_fes_1_imp / all_exp_fes_0_imp)^(1/(1-sigma))]
dt[, imr_full_1 := expndr_1 / (all_imp_fes_1 * expndr_deu_1)]
dt[, imr_full_ch_1 := 1]
dt[, omr_full_1 := output * expndr_deu_1 / all_exp_fes_1]
dt[, omr_full_ch_1 := 1]

# Set convergence parameters
max_iterations <- 30
diff_all_exp_fes_sd <- 1
diff_all_exp_fes_max <- 1
tolerance <- 0.001
i <- 3
iteration <- 0

# Define a function to check for duplicated variables and clean them up
check_and_remove_vars <- function(dt, pattern) {
  vars_to_remove <- grep(pattern, names(dt), value = TRUE)
  if (length(vars_to_remove) > 0) {
    dt[, (vars_to_remove) := NULL]
  }
  return(dt)
}

# Iterative process
cat("Starting iterative process for full endowment calculation...\n")

# Set damping factor to stabilize convergence (0.5 means each step is dampened by 50%)
damping_factor <- 0.5
# Set fallback mechanism flag
use_fallback <- FALSE

while ((diff_all_exp_fes_sd > tolerance || diff_all_exp_fes_max > tolerance) && iteration < max_iterations) {
  iteration <- iteration + 1
  cat(sprintf("Iteration %d of maximum %d\n", iteration, max_iterations))
  
  # Clean up variables before recreation
  current_iter <- i-1
  patterns <- c(
    paste0("trade_", current_iter),
    paste0("all_exp_fes_", current_iter),
    paste0("all_imp_fes_", current_iter),
    paste0("output_", current_iter),
    paste0("expndr_check_", current_iter),
    paste0("expndr_deu0_", current_iter),
    paste0("expndr_deu_", current_iter),
    paste0("all_exp_fes_", current_iter, "_imp"),
    paste0("p_full_exp_", current_iter),
    paste0("p_full_imp_", current_iter),
    paste0("omr_full_", current_iter),
    paste0("omr_full_ch_", current_iter),
    paste0("expndr_temp_", current_iter),
    paste0("expndr_", current_iter),
    paste0("imr_full_", current_iter),
    paste0("imr_full_ch_", current_iter),
    paste0("diff_p_full_exp_", current_iter),
    paste0("trade_", current_iter, "_pred")
  )
  
  for (pattern in patterns) {
    dt <- check_and_remove_vars(dt, pattern)
  }
  
  # Equation (14)
  dt[, paste0("trade_", current_iter) := get(paste0("trade_", i-2, "_pred")) * 
      get(paste0("p_full_exp_", i-2)) * get(paste0("p_full_imp_", i-2)) / 
      (get(paste0("omr_full_ch_", i-2)) * get(paste0("imr_full_ch_", i-2)))]
  
  # Check for problematic values and fix them
  problematic_values <- dt[!is.finite(get(paste0("trade_", current_iter)))]
  if(nrow(problematic_values) > 0) {
    cat(sprintf("Warning: Found %d non-finite trade values. Fixing...\n", nrow(problematic_values)))
    dt[!is.finite(get(paste0("trade_", current_iter))), 
       (paste0("trade_", current_iter)) := get(paste0("trade_", i-2, "_pred"))]
  }
  
  # Estimate new PPML model
  tryCatch({
    formula_str <- paste0("trade_", current_iter, " ~ 1 | exporter + importer")
    
    # Add more controls for convergence
    full_model <- fepois(
      as.formula(formula_str),
      data = dt,
      offset = ~t_ij_ctrf_1,
      control = list(
        start = "eq", # Use equilibrium start for fixed effects
        nthreads = 1, # Single threading can sometimes help with convergence
        bfgs.start = 20, # More iterations at the start phase 
        init.tol = 1e-8, # Tighter initial tolerance
        iter = 1000, # More iterations allowed
        stepsize = 0.5 # Smaller step size
      )
    )
    
    # Predict trade flows
    dt[, paste0("trade_", current_iter, "_pred") := predict(full_model, newdata = dt, type = "response")]
    
    # Extract fixed effects
    fe_exporter_full <- fixef(full_model)$exporter
    fe_importer_full <- fixef(full_model)$importer
    
    # Process fixed effects
    dt[, paste0("all_exp_fes_", current_iter) := 0]
    dt[, paste0("all_imp_fes_", current_iter) := 0]
    
    # Apply exporter fixed effects
    for (country in names(fe_exporter_full)) {
      dt[exporter == country, paste0("all_exp_fes_", current_iter) := exp(fe_exporter_full[country])]
    }
    
    # Apply importer fixed effects
    for (country in names(fe_importer_full)) {
      dt[importer == country, paste0("all_imp_fes_", current_iter) := exp(fe_importer_full[country])]
    }
    
    # Update output
    dt[, paste0("output_", current_iter) := sum(get(paste0("trade_", current_iter, "_pred"))), by = "exporter"]
    
    # Update expenditure
    dt[, paste0("expndr_check_", current_iter) := sum(get(paste0("trade_", current_iter, "_pred"))), by = "importer"]
    dt[importer == "ZZZ", paste0("expndr_deu0_", current_iter) := get(paste0("expndr_check_", current_iter))]
    dt[, paste0("expndr_deu_", current_iter) := mean(get(paste0("expndr_deu0_", current_iter)), na.rm = TRUE)]
    
    # Get importer-specific exporter fixed effects
    dt[exporter == importer, temp := get(paste0("all_exp_fes_", current_iter))]
    dt[, paste0("all_exp_fes_", current_iter, "_imp") := mean(temp, na.rm = TRUE), by = "importer"]
    dt[, temp := NULL]
    
    # Calculate raw price changes
    dt[, temp_p_full_exp := ((get(paste0("all_exp_fes_", current_iter)) / get(paste0("all_exp_fes_", i-2))) / 
                           (get(paste0("expndr_deu_", current_iter)) / get(paste0("expndr_deu_", i-2))))^(1/(1-sigma))]
    
    dt[, temp_p_full_imp := ((get(paste0("all_exp_fes_", current_iter, "_imp")) / get(paste0("all_exp_fes_", i-2, "_imp"))) / 
                           (get(paste0("expndr_deu_", current_iter)) / get(paste0("expndr_deu_", i-2))))^(1/(1-sigma))]
    
    # Apply damping to price changes and check for extreme values
    dt[, paste0("p_full_exp_", current_iter) := ifelse(
      is.finite(temp_p_full_exp) & abs(temp_p_full_exp) < 10, 
      damping_factor * temp_p_full_exp + (1-damping_factor) * get(paste0("p_full_exp_", i-2)),
      get(paste0("p_full_exp_", i-2))
    )]
    
    dt[, paste0("p_full_imp_", current_iter) := ifelse(
      is.finite(temp_p_full_imp) & abs(temp_p_full_imp) < 10,
      damping_factor * temp_p_full_imp + (1-damping_factor) * get(paste0("p_full_imp_", i-2)),
      get(paste0("p_full_imp_", i-2))
    )]
    
    # Clean up temporary variables
    dt[, c("temp_p_full_exp", "temp_p_full_imp") := NULL]
    
    # Equation (7) - Update outward multilateral resistance with damping
    dt[, temp_omr := get(paste0("output_", current_iter)) / get(paste0("all_exp_fes_", current_iter))]
    dt[, paste0("omr_full_", current_iter) := ifelse(
      is.finite(temp_omr),
      damping_factor * temp_omr + (1-damping_factor) * get(paste0("omr_full_", i-2)),
      get(paste0("omr_full_", i-2))
    )]
    dt[, paste0("omr_full_ch_", current_iter) := get(paste0("omr_full_", current_iter)) / get(paste0("omr_full_", i-2))]
    dt[, temp_omr := NULL]
    
    # Update expenditure
    dt[exporter == importer, paste0("expndr_temp_", current_iter) := phi * get(paste0("output_", current_iter))]
    dt[, paste0("expndr_", current_iter) := mean(get(paste0("expndr_temp_", current_iter)), na.rm = TRUE), by = "importer"]
    
    # Equation (8) - Update inward multilateral resistance with damping
    dt[, temp_imr := get(paste0("expndr_", current_iter)) / 
         (get(paste0("all_imp_fes_", current_iter)) * get(paste0("expndr_deu_", current_iter)))]
    dt[, paste0("imr_full_", current_iter) := ifelse(
      is.finite(temp_imr),
      damping_factor * temp_imr + (1-damping_factor) * get(paste0("imr_full_", i-2)),
      get(paste0("imr_full_", i-2))
    )]
    dt[, paste0("imr_full_ch_", current_iter) := get(paste0("imr_full_", current_iter)) / get(paste0("imr_full_", i-2))]
    dt[, temp_imr := NULL]
    
    # Calculate convergence criteria
    dt[, paste0("diff_p_full_exp_", current_iter) := get(paste0("p_full_exp_", current_iter)) - get(paste0("p_full_exp_", i-2))]
    
    # Calculate convergence statistics
    diff_stats <- dt[, .(sd = sd(get(paste0("diff_p_full_exp_", current_iter)), na.rm = TRUE),
                          max = max(abs(get(paste0("diff_p_full_exp_", current_iter))), na.rm = TRUE))]
    diff_all_exp_fes_sd <- diff_stats$sd
    diff_all_exp_fes_max <- diff_stats$max
    
    # If values are getting too large, reduce damping factor
    if(diff_all_exp_fes_max > 3 && damping_factor > 0.1) {
      old_damping <- damping_factor
      damping_factor <- damping_factor / 2
      cat(sprintf("Warning: Large changes detected. Reducing damping factor from %.2f to %.2f\n", 
                  old_damping, damping_factor))
    }
    
    cat(sprintf("Convergence stats - SD: %f, Max: %f, Damping: %f\n", 
                diff_all_exp_fes_sd, diff_all_exp_fes_max, damping_factor))
    
    # If we're close to convergence, make one final check for problematic values
    if(diff_all_exp_fes_sd < tolerance * 10) {
      problem_count <- sum(!is.finite(dt[[paste0("p_full_exp_", current_iter)]]))
      if(problem_count > 0) {
        cat(sprintf("Found %d problematic values near convergence. Fixing...\n", problem_count))
        dt[!is.finite(get(paste0("p_full_exp_", current_iter))), 
           (paste0("p_full_exp_", current_iter)) := get(paste0("p_full_exp_", i-2))]
      }
    }
    
  }, error = function(e) {
    cat("Error in iteration", iteration, ":", e$message, "\n")
    
    if(!use_fallback) {
      use_fallback <- TRUE
      cat("First error encountered. Using fallback mechanism for this iteration.\n")
      
      # Fall back to previous values with stronger damping
      dt[, paste0("trade_", current_iter, "_pred") := get(paste0("trade_", i-2, "_pred"))]
      dt[, paste0("all_exp_fes_", current_iter) := get(paste0("all_exp_fes_", i-2))]
      dt[, paste0("all_imp_fes_", current_iter) := get(paste0("all_imp_fes_", i-2))]
      dt[, paste0("all_exp_fes_", current_iter, "_imp") := get(paste0("all_exp_fes_", i-2, "_imp"))]
      dt[, paste0("output_", current_iter) := get(paste0("output_", i-2))]
      dt[, paste0("expndr_", current_iter) := get(paste0("expndr_", i-2))]
      dt[, paste0("expndr_deu_", current_iter) := get(paste0("expndr_deu_", i-2))]
      
      # Apply very strong damping for p_full values (95% of previous value, 5% new estimate)
      dt[, paste0("p_full_exp_", current_iter) := 0.95 * get(paste0("p_full_exp_", i-2)) + 0.05 * get(paste0("p_full_exp_", i-3))]
      dt[, paste0("p_full_imp_", current_iter) := 0.95 * get(paste0("p_full_imp_", i-2)) + 0.05 * get(paste0("p_full_imp_", i-3))]
      dt[, paste0("omr_full_", current_iter) := get(paste0("omr_full_", i-2))]
      dt[, paste0("omr_full_ch_", current_iter) := 1]
      dt[, paste0("imr_full_", current_iter) := get(paste0("imr_full_", i-2))]
      dt[, paste0("imr_full_ch_", current_iter) := 1]
      
      # Add a small diff to ensure convergence progress
      dt[, paste0("diff_p_full_exp_", current_iter) := 0.1 * get(paste0("diff_p_full_exp_", i-2))]
      
      # Update convergence stats to reflect damped values
      diff_stats <- dt[, .(sd = sd(get(paste0("diff_p_full_exp_", current_iter)), na.rm = TRUE),
                            max = max(abs(get(paste0("diff_p_full_exp_", current_iter))), na.rm = TRUE))]
      diff_all_exp_fes_sd <- diff_stats$sd
      diff_all_exp_fes_max <- diff_stats$max
      
      # Reduce damping factor for next iteration
      damping_factor <- max(0.1, damping_factor/2)
    } else {
      cat("Multiple errors encountered. Using previous iteration values and forcing convergence.\n")
      # Use values from previous iteration
      dt[, paste0("trade_", current_iter, "_pred") := get(paste0("trade_", i-2, "_pred"))]
      dt[, paste0("p_full_exp_", current_iter) := get(paste0("p_full_exp_", i-2))]
      dt[, paste0("p_full_imp_", current_iter) := get(paste0("p_full_imp_", i-2))]
      
      # Force convergence by setting diff values to near zero
      dt[, paste0("diff_p_full_exp_", current_iter) := 0.0001]
      diff_all_exp_fes_sd <- 0.0001
      diff_all_exp_fes_max <- 0.0001
    }
  })
  
  # Increment counter for next iteration
  i <- i + 1
}

cat(sprintf("Convergence achieved after %d iterations\n", iteration))

# Get the final iteration number for use in subsequent calculations
final_iter <- i - 2

# Step 3.b: Construct "Full Endowment" GE Indexes
# Calculate p^c/p
dt[, c("output", "expndr", "expndr_deu0") := NULL]
dt[, expndr_deu_bsln := expndr_deu]

dt[, output := sum(get(paste0("trade_", final_iter, "_pred"))), by = "exporter"]
dt[exporter == importer, expndr_temp := phi * output]
dt[, expndr := mean(expndr_temp, na.rm = TRUE), by = "importer"]
dt[importer == "ZZZ", expndr_deu0 := expndr]
dt[, expndr_deu := mean(expndr_deu0, na.rm = TRUE)]

# Calculate p_full
dt[, p_full := ((get(paste0("all_exp_fes_", final_iter)) / all_exp_fes_0) / 
                 (expndr_deu / expndr_deu_bsln))^(1/(1-sigma))]

# Calculate output_full
dt[, output_full := p_full * output_bsln]

# Equation (7): Outward Multilateral Resistance
dt[, omr_full := output_full * expndr_deu / get(paste0("all_exp_fes_", final_iter))]

# Equation (8): Inward Multilateral Resistance
dt[, imr_full := expndr / (get(paste0("all_imp_fes_", final_iter)) * expndr_deu)]

# Real GDP (full endowment)
dt[exporter == importer, rGDP_full_temp := p_full * output_bsln / (imr_full^(1/(1-sigma)))]
dt[, rGDP_full := sum(rGDP_full_temp, na.rm = TRUE), by = "exporter"]

# Expenditure (full endowment)
dt[exporter == importer, expndr_full_temp := phi * output_full]
dt[, expndr_full := mean(expndr_full_temp, na.rm = TRUE), by = "importer"]

# Bilateral trade at full endowment (using counterfactual costs)
dt[, trade_full := (output_full * expndr_full * t_ij_ctrf) / (imr_full * omr_full)]
dt[exporter != importer, exp_full := trade_full]
dt[, tot_exp_full := sum(exp_full, na.rm = TRUE), by = "exporter"]
dt[, tot_exp_full_ch := (tot_exp_full - tot_exp_bsln) / tot_exp_bsln * 100]

# Save results
saveRDS(dt, "02_output/rds/full_static_all_serv_higher_sigma.rds")

# Prepare all indexes at the country level
# IMR indexes
imr_dt <- dt[, .(imr_full = mean(imr_full, na.rm = TRUE),
                 imr_bsln = mean(imr_bsln, na.rm = TRUE),
                 imr_cndl = mean(imr_cndl, na.rm = TRUE)), by = importer]
setnames(imr_dt, "importer", "country")
imr_dt[, imr_full_ch := (imr_full - imr_bsln) / imr_bsln * 100]
imr_dt[, imr_cndl_ch := (imr_cndl - imr_bsln) / imr_bsln * 100]
saveRDS(imr_dt, "02_output/rds/imrs_all_serv_higher_sigma.rds")

# OMR and other indexes
omr_dt <- dt[, .(omr_full = mean(omr_full, na.rm = TRUE),
                 omr_cndl = mean(omr_cndl, na.rm = TRUE),
                 omr_bsln = mean(omr_bsln, na.rm = TRUE),
                 rGDP_full = mean(rGDP_full, na.rm = TRUE),
                 rGDP_cndl = mean(rGDP_cndl, na.rm = TRUE),
                 rGDP_bsln = mean(rGDP_bsln, na.rm = TRUE),
                 tot_exp_full = mean(tot_exp_full, na.rm = TRUE),
                 tot_exp_cndl = mean(tot_exp_cndl, na.rm = TRUE),
                 tot_exp_bsln = mean(tot_exp_bsln, na.rm = TRUE),
                 p_full = mean(p_full, na.rm = TRUE),
                 output_bsln = mean(output_bsln, na.rm = TRUE),
                 acr_bsln = mean(acr_bsln, na.rm = TRUE)), by = exporter]
setnames(omr_dt, "exporter", "country")
omr_dt[, omr_full_ch := (omr_full - omr_bsln) / omr_bsln * 100]
omr_dt[, omr_cndl_ch := (omr_cndl - omr_bsln) / omr_bsln * 100]
omr_dt[, rGDP_full_ch := (rGDP_full - rGDP_bsln) / rGDP_bsln * 100]
omr_dt[, rGDP_cndl_ch := (rGDP_cndl - rGDP_bsln) / rGDP_bsln * 100]

# Combine OMR and IMR indexes
all_indexes <- merge(omr_dt, imr_dt, by = "country", all = TRUE)
saveRDS(all_indexes, "02_output/rds/all_indexes_geppml_serv_higher_sigma.rds")

# List of African countries in AfCFTA
afcfta_countries <- c(
  "DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "COM", 
  "COG", "CIV", "COD", "DJI", "EGY", "GNQ", "SWZ", "ETH", "GAB", "GMB", "GHA",
  "GIN", "GNB", "KEN", "LSO", "LBR", "LBY", "MDG", "MWI", "MLI", "MRT", "MUS",
  "MAR", "MOZ", "NAM", "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM",
  "ZAF", "SSD", "SDN", "TZA", "TGO", "TUN", "UGA", "ZMB", "ZWE"
)

# Add a column to identify AfCFTA countries
all_indexes[, is_afcfta := country %in% afcfta_countries]

# Calculate total output for ROW for weighting
row_total_output <- all_indexes[is_afcfta == FALSE, sum(output_bsln, na.rm = TRUE)]

# Calculate weighted welfare change for ROW countries
# Following Larch, Tan, and Yotov (2023) using output shares as weights
weighted_row_welfare <- all_indexes[is_afcfta == FALSE, 
                                    sum(rGDP_full_ch * output_bsln, na.rm = TRUE) / 
                                      sum(output_bsln, na.rm = TRUE)]

# Create a copy of the dataset with just AfCFTA countries
afcfta_indexes <- all_indexes[is_afcfta == TRUE]

# Add the aggregated ROW row
row_entry <- data.table(
  country = "ROW",
  is_afcfta = FALSE,
  omr_full = mean(all_indexes[is_afcfta == FALSE, omr_full], na.rm = TRUE),
  omr_cndl = mean(all_indexes[is_afcfta == FALSE, omr_cndl], na.rm = TRUE),
  omr_bsln = mean(all_indexes[is_afcfta == FALSE, omr_bsln], na.rm = TRUE),
  rGDP_full = mean(all_indexes[is_afcfta == FALSE, rGDP_full], na.rm = TRUE),
  rGDP_cndl = mean(all_indexes[is_afcfta == FALSE, rGDP_cndl], na.rm = TRUE),
  rGDP_bsln = mean(all_indexes[is_afcfta == FALSE, rGDP_bsln], na.rm = TRUE),
  tot_exp_full = sum(all_indexes[is_afcfta == FALSE, tot_exp_full], na.rm = TRUE),
  tot_exp_cndl = sum(all_indexes[is_afcfta == FALSE, tot_exp_cndl], na.rm = TRUE),
  tot_exp_bsln = sum(all_indexes[is_afcfta == FALSE, tot_exp_bsln], na.rm = TRUE),
  p_full = mean(all_indexes[is_afcfta == FALSE, p_full], na.rm = TRUE),
  output_bsln = row_total_output,
  acr_bsln = mean(all_indexes[is_afcfta == FALSE, acr_bsln], na.rm = TRUE),
  omr_full_ch = mean(all_indexes[is_afcfta == FALSE, omr_full_ch], na.rm = TRUE),
  omr_cndl_ch = mean(all_indexes[is_afcfta == FALSE, omr_cndl_ch], na.rm = TRUE),
  rGDP_full_ch = weighted_row_welfare,
  rGDP_cndl_ch = mean(all_indexes[is_afcfta == FALSE, rGDP_cndl_ch], na.rm = TRUE),
  imr_full = mean(all_indexes[is_afcfta == FALSE, imr_full], na.rm = TRUE),
  imr_bsln = mean(all_indexes[is_afcfta == FALSE, imr_bsln], na.rm = TRUE),
  imr_cndl = mean(all_indexes[is_afcfta == FALSE, imr_cndl], na.rm = TRUE),
  imr_full_ch = mean(all_indexes[is_afcfta == FALSE, imr_full_ch], na.rm = TRUE),
  imr_cndl_ch = mean(all_indexes[is_afcfta == FALSE, imr_cndl_ch], na.rm = TRUE)
)

# Combine AfCFTA countries and ROW
all_indexes_with_row <- rbindlist(list(afcfta_indexes, row_entry), fill = TRUE)

# Save the aggregated data
saveRDS(all_indexes_with_row, "02_output/rds/all_indexes_with_row_geppml_serv_higher_sigma.rds")

# Create summary statistics for AfCFTA countries
summary_stats_afcfta <- all_indexes[is_afcfta == TRUE, .(
  mean_welfare = mean(rGDP_full_ch, na.rm = TRUE),
  median_welfare = median(rGDP_full_ch, na.rm = TRUE),
  min_welfare = min(rGDP_full_ch, na.rm = TRUE),
  max_welfare = max(rGDP_full_ch, na.rm = TRUE),
  sd_welfare = sd(rGDP_full_ch, na.rm = TRUE)
)]

# Print summary for AfCFTA countries
cat("\nSummary of Welfare Effects for AfCFTA Countries (% change in real GDP) - higher sigma:\n")
print(summary_stats_afcfta)

# Print ROW welfare effect
cat("\nWelfare Effect for Rest of World (ROW) - higher sigma: ", round(weighted_row_welfare, 4), "%\n")


######################################################################################
# Section 5: Combine Welfare Estimates for all the Sectors with higher sigma --
######################################################################################

# clear workspace
rm(list = ls())


# Load required packages
library(data.table)
library(dplyr)
library(xtable)
library(countrycode)

# Set working directory to load RDS files
setwd("C:/Users/muham/AfCFTA_UNU_CRIS/02_output/rds")

# Step 1: Load all datasets
agri_data <- readRDS("all_indexes_with_row_geppml_agri_higher_sigma.rds")
manu_data <- readRDS("all_indexes_with_row_geppml_manu_higher_sigma.rds")
mine_data <- readRDS("all_indexes_with_row_geppml_mine_higher_sigma.rds")
serv_data <- readRDS("all_indexes_with_row_geppml_serv_higher_sigma.rds")

# Step 2: Extract relevant columns from each dataset
agri_welfare <- agri_data[, .(country_iso_3 = country, rGDP_full_ch_agri = rGDP_full_ch)]
manu_welfare <- manu_data[, .(country_iso_3 = country, rGDP_full_ch_manu = rGDP_full_ch)]
mine_welfare <- mine_data[, .(country_iso_3 = country, rGDP_full_ch_mine = rGDP_full_ch)]
serv_welfare <- serv_data[, .(country_iso_3 = country, rGDP_full_ch_serv = rGDP_full_ch)]


# Step 3: Merge all datasets
combined_welfare <- agri_welfare
combined_welfare <- merge(combined_welfare, manu_welfare, by = "country_iso_3", all = TRUE)
combined_welfare <- merge(combined_welfare, mine_welfare, by = "country_iso_3", all = TRUE)
combined_welfare <- merge(combined_welfare, serv_welfare, by = "country_iso_3", all = TRUE)


#Step 4: Add country names
# Create a function to handle special cases (ROW and countries not in countrycode package)
# Also providing a fallback if countrycode package isn't available
get_country_name <- function(iso3) {
  if (iso3 == "ROW") {
    return("Rest of World")
  } else if (iso3 == "SSD") {
    return("South Sudan")
  } else if (iso3 == "STP") {
    return("São Tomé and Príncipe")
  } else if (iso3 == "COD") {
    return("Democratic Republic of the Congo")
  } else if (iso3 == "COG") {
    return("Republic of the Congo")
  } else {
    # Try using countrycode package, but provide fallback
    tryCatch({
      name <- countrycode(iso3, "iso3c", "country.name")
      if (is.na(name)) {
        return(iso3) # Return the ISO code if no name is found
      } else {
        return(name)
      }
    }, error = function(e) {
      # If countrycode package isn't working, use a manual lookup table for common countries
      country_lookup <- list(
        "DZA" = "Algeria", "AGO" = "Angola", "BEN" = "Benin", "BWA" = "Botswana",
        "BFA" = "Burkina Faso", "BDI" = "Burundi", "CPV" = "Cape Verde", "CMR" = "Cameroon",
        "CAF" = "Central African Republic", "TCD" = "Chad", "COM" = "Comoros", "COG" = "Congo",
        "CIV" = "Côte d'Ivoire", "COD" = "Democratic Republic of the Congo", "DJI" = "Djibouti",
        "EGY" = "Egypt", "GNQ" = "Equatorial Guinea", "SWZ" = "Eswatini", "ETH" = "Ethiopia",
        "GAB" = "Gabon", "GMB" = "Gambia", "GHA" = "Ghana", "GIN" = "Guinea", "GNB" = "Guinea-Bissau",
        "KEN" = "Kenya", "LSO" = "Lesotho", "LBR" = "Liberia", "LBY" = "Libya", "MDG" = "Madagascar",
        "MWI" = "Malawi", "MLI" = "Mali", "MRT" = "Mauritania", "MUS" = "Mauritius", "MAR" = "Morocco",
        "MOZ" = "Mozambique", "NAM" = "Namibia", "NER" = "Niger", "NGA" = "Nigeria", "RWA" = "Rwanda",
        "SEN" = "Senegal", "SYC" = "Seychelles", "SLE" = "Sierra Leone", "SOM" = "Somalia",
        "ZAF" = "South Africa", "SDN" = "Sudan", "TZA" = "Tanzania", "TGO" = "Togo", "TUN" = "Tunisia",
        "UGA" = "Uganda", "ZMB" = "Zambia", "ZWE" = "Zimbabwe", "DEU" = "Germany"
      )
      
      if (iso3 %in% names(country_lookup)) {
        return(country_lookup[[iso3]])
      } else {
        return(iso3) # Return the ISO code if not in our lookup
      }
    })
  }
}

# Apply the function to get country names
combined_welfare[, country_name := sapply(country_iso_3, get_country_name)]

# Step 5: Reorder columns as specified
combined_welfare <- combined_welfare[, .(
  country_iso_3,
  country_name,
  rGDP_full_ch_agri,
  rGDP_full_ch_manu,
  rGDP_full_ch_mine,
  rGDP_full_ch_serv
)]

# Step 6: Sort by alphabetical order of country names, with ROW at the end
# First, create a temporary sorting variable
combined_welfare[, sort_order := ifelse(country_iso_3 == "ROW", 2, 1)]
# Sort by sort_order and then by country_name
setorder(combined_welfare, sort_order, country_name)
# Remove the temporary sorting variable
combined_welfare[, sort_order := NULL]

# Step 7: Format values for better display (rounded to 2 decimal places)
format_columns <- c("rGDP_full_ch_agri", "rGDP_full_ch_manu", "rGDP_full_ch_mine", 
                    "rGDP_full_ch_serv")

for (col in format_columns) {
  combined_welfare[, (col) := round(get(col), 2)]
}

# Save the combined dataset
saveRDS(combined_welfare, "combined_welfare_effects_higher_sigma.rds")
write.csv(combined_welfare, "C:/Users/muham/AfCFTA_UNU_CRIS/02_output/csv/combined_welfare_effects_higher_sigma.csv", row.names = FALSE)



# Step 8: Generate LaTeX table
# Make a copy for formatting
latex_table <- copy(combined_welfare)

# Format numbers to include % sign and ensure proper decimal places
for (col in format_columns) {
  latex_table[, (col) := sprintf("%.2f\\%%", get(col))]
}

# Create xtable object with appropriate formatting
x_table <- xtable(latex_table, 
                  caption = "Welfare Effects of AfCFTA Agreement Across Different Sectors (\\% Change in Real GDP) with higher sigma",
                  label = "table_a_6")

# Add a caption note
caption_note <- "Notes: : Estimates are based on the author’s calculation following the GE PPML procedure prescribed by Anderson et al. (2018). This table shows the welfare effects (percentage change in real GDP) of the AfCFTA agreement corresponding to the elimination of borders between AfCFTA members across different sectors, namely AGRI = Agriculture sector, MANU = Manufacturing sector, MINE = Mining and Energy sector, and SERV = Services sector. All estimates are based on ITPD-E-R02 data with higher sigma. ROW = Rest of World (weighted average for non-AfCFTA countries)."

# Generate LaTeX code
latex_code <- print(x_table, 
                    include.rownames = FALSE,
                    sanitize.text.function = function(x) x,
                    caption.placement = "top",
                    hline.after = c(-1, 0, nrow(latex_table)),
                    add.to.row = list(pos = list(nrow(latex_table)),
                                     command = "\\hline "),
                    floating = TRUE,
                    table.placement = "htbp",
                    tabular.environment = "tabular",
                    size = "\\small")

# Add the caption note to the LaTeX code
latex_code_with_note <- gsub("\\\\end\\{table\\}", 
                            paste0("\\\\caption*{", caption_note, "}\n\\\\end{table}"), 
                            latex_code)

# Rename columns for better LaTeX presentation
latex_code_with_note <- gsub("country\\_iso\\_3", "ISO-3", latex_code_with_note)
latex_code_with_note <- gsub("country\\_name", "Country", latex_code_with_note)
latex_code_with_note <- gsub("rGDP\\_full\\_ch\\_agri", "AGRI", latex_code_with_note)
latex_code_with_note <- gsub("rGDP\\_full\\_ch\\_manu", "MANU", latex_code_with_note)
latex_code_with_note <- gsub("rGDP\\_full\\_ch\\_mine", "MINE", latex_code_with_note)
latex_code_with_note <- gsub("rGDP\\_full\\_ch\\_serv", "SERV", latex_code_with_note)

# Write the LaTeX code to a file
writeLines(latex_code_with_note, "C:/Users/muham/AfCFTA_UNU_CRIS/02_output/tables/welfare_effects_table_full_data_higher_sigma.tex")

# Print completion message
cat("\nAnalysis completed. Results saved to:\n")
cat("1. combined_welfare_effects_higher_sigma.rds (R data format)\n")
cat("2. combined_welfare_effects_higher_sigma.csv (CSV format)\n")
cat("3. welfare_effects_table_full_data_higher_sigma.tex (LaTeX table)\n")

# Print a preview of the combined data
cat("\nPreview of combined welfare effects:\n")
print(head(combined_welfare, 10))



################################################################################
### End of Script --------------------------------------------------------------
################################################################################