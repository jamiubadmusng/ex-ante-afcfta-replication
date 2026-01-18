################################################################################
# R code for General Equilibrium PPML Analysis of AfCFTA
# This code translates the original Stata GE PPML code by Anderson et al. (2018)
# for the analysis of abolishing borders between AfCFTA countries for Manufacturing
# Trade flows from Structural Gravity Database with sigma = 4, 5, 7, & 10.
# Author: Jamiu Olamilekan Badmus
################################################################################


#################################################################################
# Estimation of GE PPML Model for Structural Gravity Data with Sigma = 4 --------
#################################################################################

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
sigma <- 4

################################################################################
# I. Prepare Data
################################################################################

# Load data
data <- read.dta13("01_input/struc_2016_balanced.dta")
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
saveRDS(dt, "01_input/ge_ppml_data_struc.rds")

################################################################################
# II. GE Analysis in R
################################################################################

# Load prepared data
dt <- readRDS("01_input/ge_ppml_data_struc.rds")

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
saveRDS(dt, "02_output/rds/full_static_all_struc_sigma4.rds")

# Prepare all indexes at the country level
# IMR indexes
imr_dt <- dt[, .(imr_full = mean(imr_full, na.rm = TRUE),
                 imr_bsln = mean(imr_bsln, na.rm = TRUE),
                 imr_cndl = mean(imr_cndl, na.rm = TRUE)), by = importer]
setnames(imr_dt, "importer", "country")
imr_dt[, imr_full_ch := (imr_full - imr_bsln) / imr_bsln * 100]
imr_dt[, imr_cndl_ch := (imr_cndl - imr_bsln) / imr_bsln * 100]
saveRDS(imr_dt, "02_output/rds/imrs_all_struc_sigma4.rds")

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
saveRDS(all_indexes, "02_output/rds/all_indexes_geppml_struc_sigma4.rds")


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
saveRDS(all_indexes_with_row, "02_output/rds/all_indexes_with_row_geppml_struc_sigma4.rds")

# Create summary statistics for AfCFTA countries
summary_stats_afcfta <- all_indexes[is_afcfta == TRUE, .(
  mean_welfare = mean(rGDP_full_ch, na.rm = TRUE),
  median_welfare = median(rGDP_full_ch, na.rm = TRUE),
  min_welfare = min(rGDP_full_ch, na.rm = TRUE),
  max_welfare = max(rGDP_full_ch, na.rm = TRUE),
  sd_welfare = sd(rGDP_full_ch, na.rm = TRUE)
)]

# Print summary for AfCFTA countries
cat("\nSummary of Welfare Effects for AfCFTA Countries (% change in real GDP) - sigma = 4:\n")
print(summary_stats_afcfta)

# Print ROW welfare effect
cat("\nWelfare Effect for Rest of World (ROW) - sigma = 4: ", round(weighted_row_welfare, 4), "%\n")


#################################################################################
# Estimation of GE PPML Model for Structural Gravity Data with Sigma = 5 --------
#################################################################################

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
sigma <- 5

################################################################################
# I. Prepare Data
################################################################################

# Load data
data <- read.dta13("01_input/struc_2016_balanced.dta")
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
saveRDS(dt, "01_input/ge_ppml_data_struc.rds")

################################################################################
# II. GE Analysis in R
################################################################################

# Load prepared data
dt <- readRDS("01_input/ge_ppml_data_struc.rds")

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
saveRDS(dt, "02_output/rds/full_static_all_struc_sigma5.rds")

# Prepare all indexes at the country level
# IMR indexes
imr_dt <- dt[, .(imr_full = mean(imr_full, na.rm = TRUE),
                 imr_bsln = mean(imr_bsln, na.rm = TRUE),
                 imr_cndl = mean(imr_cndl, na.rm = TRUE)), by = importer]
setnames(imr_dt, "importer", "country")
imr_dt[, imr_full_ch := (imr_full - imr_bsln) / imr_bsln * 100]
imr_dt[, imr_cndl_ch := (imr_cndl - imr_bsln) / imr_bsln * 100]
saveRDS(imr_dt, "02_output/rds/imrs_all_struc_sigma5.rds")

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
saveRDS(all_indexes, "02_output/rds/all_indexes_geppml_struc_sigma5.rds")


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
saveRDS(all_indexes_with_row, "02_output/rds/all_indexes_with_row_geppml_struc_sigma5.rds")

# Create summary statistics for AfCFTA countries
summary_stats_afcfta <- all_indexes[is_afcfta == TRUE, .(
  mean_welfare = mean(rGDP_full_ch, na.rm = TRUE),
  median_welfare = median(rGDP_full_ch, na.rm = TRUE),
  min_welfare = min(rGDP_full_ch, na.rm = TRUE),
  max_welfare = max(rGDP_full_ch, na.rm = TRUE),
  sd_welfare = sd(rGDP_full_ch, na.rm = TRUE)
)]

# Print summary for AfCFTA countries
cat("\nSummary of Welfare Effects for AfCFTA Countries (% change in real GDP) - sigma = 5:\n")
print(summary_stats_afcfta)

# Print ROW welfare effect
cat("\nWelfare Effect for Rest of World (ROW) - sigma = 5: ", round(weighted_row_welfare, 4), "%\n")



#################################################################################
# Estimation of GE PPML Model for Structural Gravity Data with Sigma = 7 --------
#################################################################################

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
sigma <- 7

################################################################################
# I. Prepare Data
################################################################################

# Load data
data <- read.dta13("01_input/struc_2016_balanced.dta")
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
saveRDS(dt, "01_input/ge_ppml_data_struc.rds")

################################################################################
# II. GE Analysis in R
################################################################################

# Load prepared data
dt <- readRDS("01_input/ge_ppml_data_struc.rds")

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
saveRDS(dt, "02_output/rds/full_static_all_struc_sigma7.rds")

# Prepare all indexes at the country level
# IMR indexes
imr_dt <- dt[, .(imr_full = mean(imr_full, na.rm = TRUE),
                 imr_bsln = mean(imr_bsln, na.rm = TRUE),
                 imr_cndl = mean(imr_cndl, na.rm = TRUE)), by = importer]
setnames(imr_dt, "importer", "country")
imr_dt[, imr_full_ch := (imr_full - imr_bsln) / imr_bsln * 100]
imr_dt[, imr_cndl_ch := (imr_cndl - imr_bsln) / imr_bsln * 100]
saveRDS(imr_dt, "02_output/rds/imrs_all_struc_sigma7.rds")

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
saveRDS(all_indexes, "02_output/rds/all_indexes_geppml_struc_sigma7.rds")


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
saveRDS(all_indexes_with_row, "02_output/rds/all_indexes_with_row_geppml_struc_sigma7.rds")

# Create summary statistics for AfCFTA countries
summary_stats_afcfta <- all_indexes[is_afcfta == TRUE, .(
  mean_welfare = mean(rGDP_full_ch, na.rm = TRUE),
  median_welfare = median(rGDP_full_ch, na.rm = TRUE),
  min_welfare = min(rGDP_full_ch, na.rm = TRUE),
  max_welfare = max(rGDP_full_ch, na.rm = TRUE),
  sd_welfare = sd(rGDP_full_ch, na.rm = TRUE)
)]

# Print summary for AfCFTA countries
cat("\nSummary of Welfare Effects for AfCFTA Countries (% change in real GDP) - sigma = 7:\n")
print(summary_stats_afcfta)

# Print ROW welfare effect
cat("\nWelfare Effect for Rest of World (ROW) - sigma = 7: ", round(weighted_row_welfare, 4), "%\n")



#################################################################################
# Estimation of GE PPML Model for Structural Gravity Data with Sigma = 10 -------
#################################################################################

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
# I. Prepare Data
################################################################################

# Load data
data <- read.dta13("01_input/struc_2016_balanced.dta")
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
saveRDS(dt, "01_input/ge_ppml_data_struc.rds")

################################################################################
# II. GE Analysis in R
################################################################################

# Load prepared data
dt <- readRDS("01_input/ge_ppml_data_struc.rds")

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
saveRDS(dt, "02_output/rds/full_static_all_struc_sigma10.rds")

# Prepare all indexes at the country level
# IMR indexes
imr_dt <- dt[, .(imr_full = mean(imr_full, na.rm = TRUE),
                 imr_bsln = mean(imr_bsln, na.rm = TRUE),
                 imr_cndl = mean(imr_cndl, na.rm = TRUE)), by = importer]
setnames(imr_dt, "importer", "country")
imr_dt[, imr_full_ch := (imr_full - imr_bsln) / imr_bsln * 100]
imr_dt[, imr_cndl_ch := (imr_cndl - imr_bsln) / imr_bsln * 100]
saveRDS(imr_dt, "02_output/rds/imrs_all_struc_sigma10.rds")

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
saveRDS(all_indexes, "02_output/rds/all_indexes_geppml_struc_sigma10.rds")


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
saveRDS(all_indexes_with_row, "02_output/rds/all_indexes_with_row_geppml_struc_sigma10.rds")

# Create summary statistics for AfCFTA countries
summary_stats_afcfta <- all_indexes[is_afcfta == TRUE, .(
  mean_welfare = mean(rGDP_full_ch, na.rm = TRUE),
  median_welfare = median(rGDP_full_ch, na.rm = TRUE),
  min_welfare = min(rGDP_full_ch, na.rm = TRUE),
  max_welfare = max(rGDP_full_ch, na.rm = TRUE),
  sd_welfare = sd(rGDP_full_ch, na.rm = TRUE)
)]

# Print summary for AfCFTA countries
cat("\nSummary of Welfare Effects for AfCFTA Countries (% change in real GDP) - sigma = 10:\n")
print(summary_stats_afcfta)

# Print ROW welfare effect
cat("\nWelfare Effect for Rest of World (ROW) - sigma = 10: ", round(weighted_row_welfare, 4), "%\n")


######################################################################################
# Combine Welfare Estimates for all the Sectors with sigma = 4, 5, 7, and 10 ---------
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
struc_sigma4 <- readRDS("all_indexes_with_row_geppml_struc_sigma4.rds")
struc_sigma5 <- readRDS("all_indexes_with_row_geppml_struc_sigma5.rds")
struc_sigma7 <- readRDS("all_indexes_with_row_geppml_struc_sigma7.rds")
struc_sigma10 <- readRDS("all_indexes_with_row_geppml_struc_sigma10.rds")

# Step 2: Extract relevant columns from each dataset
struc_sigma04 <- struc_sigma4[, .(country_iso_3 = country, rGDP_full_ch_sigma4 = rGDP_full_ch)]
struc_sigma05 <- struc_sigma5[, .(country_iso_3 = country, rGDP_full_ch_sigma5 = rGDP_full_ch)]
struc_sigma07 <- struc_sigma7[, .(country_iso_3 = country, rGDP_full_ch_sigma7 = rGDP_full_ch)]
struc_sigma010 <- struc_sigma10[, .(country_iso_3 = country, rGDP_full_ch_sigma10 = rGDP_full_ch)]


# Step 3: Merge all datasets
combined_welfare <- struc_sigma04
combined_welfare <- merge(combined_welfare, struc_sigma05, by = "country_iso_3", all = TRUE)
combined_welfare <- merge(combined_welfare, struc_sigma07, by = "country_iso_3", all = TRUE)
combined_welfare <- merge(combined_welfare, struc_sigma010, by = "country_iso_3", all = TRUE)


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
  rGDP_full_ch_sigma4,
  rGDP_full_ch_sigma5,
  rGDP_full_ch_sigma7,
  rGDP_full_ch_sigma10
)]

# Step 6: Sort by alphabetical order of country names, with ROW at the end
# First, create a temporary sorting variable
combined_welfare[, sort_order := ifelse(country_iso_3 == "ROW", 2, 1)]
# Sort by sort_order and then by country_name
setorder(combined_welfare, sort_order, country_name)
# Remove the temporary sorting variable
combined_welfare[, sort_order := NULL]

# Step 7: Format values for better display (rounded to 2 decimal places)
format_columns <- c("rGDP_full_ch_sigma4", "rGDP_full_ch_sigma5", "rGDP_full_ch_sigma7", 
                    "rGDP_full_ch_sigma10")

for (col in format_columns) {
  combined_welfare[, (col) := round(get(col), 2)]
}

# Save the combined dataset
saveRDS(combined_welfare, "combined_welfare_effects_struc_data.rds")
write.csv(combined_welfare, "C:/Users/muham/AfCFTA_UNU_CRIS/02_output/csv/combined_welfare_effects_struc_data.csv", row.names = FALSE)



# Step 8: Generate LaTeX table
# Make a copy for formatting
latex_table <- copy(combined_welfare)

# Format numbers to include % sign and ensure proper decimal places
for (col in format_columns) {
  latex_table[, (col) := sprintf("%.2f\\%%", get(col))]
}

# Create xtable object with appropriate formatting
x_table <- xtable(latex_table, 
                  caption = "Welfare Effects of AfCFTA Agreement Across Different Sectors (\\% Change in Real GDP) with the Structural Gravity Database",
                  label = "table_a_8")

# Add a caption note
caption_note <- "Notes: : Estimates are based on the author’s calculation following the GE PPML procedure prescribed by Anderson et al. (2018). This table shows the welfare effects (percentage change in real GDP) of the AfCFTA agreement corresponding to the elimination of borders between AfCFTA members using the manufacturing trade flows from the Structural Gravity Database. ROW = Rest of World (weighted average for non-AfCFTA countries)."

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
latex_code_with_note <- gsub("rGDP\\_full\\_ch\\_sigma4", "sigma4", latex_code_with_note)
latex_code_with_note <- gsub("rGDP\\_full\\_ch\\_sigma5", "sigma5", latex_code_with_note)
latex_code_with_note <- gsub("rGDP\\_full\\_ch\\_sigma7", "sigma7", latex_code_with_note)
latex_code_with_note <- gsub("rGDP\\_full\\_ch\\_sigma10", "sigma10", latex_code_with_note)

# Write the LaTeX code to a file
writeLines(latex_code_with_note, "C:/Users/muham/AfCFTA_UNU_CRIS/02_output/tables/welfare_effects_table_struc_data.tex")
  
# Print a preview of the combined data
cat("\nPreview of combined welfare effects:\n")
print(head(combined_welfare, 10))






################################################################################
# World Map showing the welfare changes for all countries in all_indexes
################################################################################
# Create world map visualization of welfare effects
################################################################################
cat("\nCreating world map visualization of welfare effects...\n")

# Install required packages if not already installed
if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
  cat("Installing rnaturalearth package...\n")
  install.packages("rnaturalearth")
}
if (!requireNamespace("rnaturalearthdata", quietly = TRUE)) {
  cat("Installing rnaturalearthdata package...\n")
  install.packages("rnaturalearthdata")
}
if (!requireNamespace("sf", quietly = TRUE)) {
  cat("Installing sf package...\n")
  install.packages("sf")
}
if (!requireNamespace("viridis", quietly = TRUE)) {
  cat("Installing viridis package...\n")
  install.packages("viridis")
}

# Load required packages for mapping
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(viridis)  # For alternative color palettes

# Display summary statistics of the welfare effects
cat("\nSummary of all country welfare effects (% change in real GDP):\n")
print(summary(all_indexes$rGDP_full_ch))

# Fix Germany ISO code in all_indexes
cat("Converting Germany's ISO code from ZZZ to DEU...\n")
all_indexes[country == "ZZZ", country := "DEU"]

# Get world map data
cat("Loading world map data...\n")
world <- ne_countries(scale = "medium", returnclass = "sf")

# Prepare the welfare data for mapping
cat("Preparing welfare data for mapping...\n")
welfare_data <- data.table(iso_a3 = all_indexes$country,
                          welfare_change = all_indexes$rGDP_full_ch)

# Print the top countries with the highest welfare gains
cat("\nTop 10 countries with highest welfare gains:\n")
welfare_data[order(-welfare_change)][1:10, .(iso_a3, welfare_change)]

# Print the top countries with the highest welfare losses
cat("\nTop 10 countries with highest welfare losses:\n")
welfare_data[order(welfare_change)][1:10, .(iso_a3, welfare_change)]

# Merge welfare data with map data
cat("Merging welfare data with geographic data...\n")
map_data <- merge(world, welfare_data, by.x = "iso_a3", by.y = "iso_a3", all.x = TRUE)

# List of African countries in AfCFTA
afcfta_countries <- c(
  "DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "COM", 
  "COG", "CIV", "COD", "DJI", "EGY", "GNQ", "SWZ", "ETH", "GAB", "GMB", "GHA",
  "GIN", "GNB", "KEN", "LSO", "LBR", "LBY", "MDG", "MWI", "MLI", "MRT", "MUS",
  "MAR", "MOZ", "NAM", "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM",
  "ZAF", "SSD", "SDN", "TZA", "TGO", "TUN", "UGA", "ZMB", "ZWE"
)

# Calculate summary statistics for AfCFTA countries
cat("\nSummary of welfare effects for AfCFTA member countries (% change in real GDP):\n")
afcfta_welfare <- welfare_data[iso_a3 %in% afcfta_countries]
print(summary(afcfta_welfare$welfare_change))

# Calculate summary statistics for non-AfCFTA countries
cat("\nSummary of welfare effects for non-AfCFTA countries (% change in real GDP):\n")
non_afcfta_welfare <- welfare_data[!iso_a3 %in% afcfta_countries]
print(summary(non_afcfta_welfare$welfare_change))

# Create directory for output if it doesn't exist
if (!dir.exists("output/figures")) {
  dir.create("output/figures", recursive = TRUE)
}

# Add a column to identify AfCFTA countries to map_data for visualization
cat("Adding country labels to maps...\n")
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
saveRDS(all_indexes_with_row, "output/rds/all_indexes_with_row_geppml_struc.rds")

# Create summary statistics for AfCFTA countries
summary_stats_afcfta <- all_indexes[is_afcfta == TRUE, .(
  mean_welfare = mean(rGDP_full_ch, na.rm = TRUE),
  median_welfare = median(rGDP_full_ch, na.rm = TRUE),
  min_welfare = min(rGDP_full_ch, na.rm = TRUE),
  max_welfare = max(rGDP_full_ch, na.rm = TRUE),
  sd_welfare = sd(rGDP_full_ch, na.rm = TRUE)
)]

# Print summary for AfCFTA countries
cat("\nSummary of Welfare Effects for AfCFTA Countries (% change in real GDP):\n")
print(summary_stats_afcfta)

# Print ROW welfare effect
cat("\nWelfare Effect for Rest of World (ROW): ", round(weighted_row_welfare, 4), "%\n")

# Create plots for welfare effects showing AfCFTA countries and ROW
# Color by whether the country is an AfCFTA member
welfare_plot <- ggplot(all_indexes_with_row, 
                      aes(x = reorder(country, rGDP_full_ch), 
                          y = rGDP_full_ch,
                          fill = is_afcfta)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2f%%", rGDP_full_ch)), 
            hjust = ifelse(all_indexes_with_row$rGDP_full_ch > 0, -0.1, 1.1),
            size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("gray60", "steelblue"), 
                    labels = c("Rest of World", "AfCFTA Countries"),
                    name = "") +
  labs(title = "Welfare Effects of AfCFTA Agreement for Manufacturing Trade from SGD",
       subtitle = "% Change in Real GDP",
       x = "Country",
       y = "% Change in Real GDP") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save the plot
ggsave("output/figures/welfare_effects_with_row_geppml_struc.pdf", welfare_plot, width = 10, height = 12, dpi = 300)

# Create another plot showing just AfCFTA countries for detail
welfare_plot_afcfta <- ggplot(all_indexes_with_row[is_afcfta == TRUE], 
                             aes(x = reorder(country, rGDP_full_ch), 
                                 y = rGDP_full_ch)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%.2f%%", rGDP_full_ch)), 
            hjust = ifelse(all_indexes_with_row[is_afcfta == TRUE]$rGDP_full_ch > 0, -0.1, 1.1),
            size = 3.5) +
  coord_flip() +
  labs(title = "Welfare Effects for AfCFTA Member Countries for Manufacturing Trade from SGD",
       subtitle = "% Change in Real GDP",
       x = "Country",
       y = "% Change in Real GDP") +
  theme_minimal()

# Save the AfCFTA-only plot
ggsave("output/figures/welfare_effects_afcfta_only_geppml_struc.pdf", welfare_plot_afcfta, width = 8, height = 12, dpi = 300)

# Add AfCFTA membership flag to the map data
map_data$is_afcfta <- map_data$iso_a3 %in% afcfta_countries

# Create the world map plot
# Define breaks and colors for welfare effects
map_data$welfare_cat <- cut(map_data$welfare_change, 
                           breaks = c(-Inf, -0.5, -0.1, 0, 0.1, 0.5, 1, 2, 10, 20, 50, 100, Inf),
                           labels = c("<-0.5%", "-0.5% to -0.1%", "-0.1% to 0%", 
                                     "0% to 0.1%", "0.1% to 0.5%", "0.5% to 1%", 
                                     "1% to 2%", "2% to 10%", "10% to 20%", "20% to 50%", "50% to 100%", ">100%"))

# Make NA values grey
map_data$welfare_cat[is.na(map_data$welfare_cat)] <- "No data"

# Prepare country label data
# Calculate centroids for country labels
map_data_labels <- map_data %>%
  st_centroid() %>%
  st_coordinates() %>%
  as.data.frame() %>%
  cbind(name = map_data$name, iso_a3 = map_data$iso_a3, 
        is_afcfta = map_data$iso_a3 %in% afcfta_countries,
        welfare_change = map_data$welfare_change)

# Determine which countries to label on the world map
# AfCFTA countries and major economies
major_economies <- c("USA", "CAN", "CHN", "JPN", "DEU", "GBR", "FRA", "IND", "BRA", "RUS", "AUS")
map_data_labels$show_on_world <- map_data_labels$is_afcfta | map_data_labels$iso_a3 %in% major_economies

# Plotting with a custom color scale - Red to Blue colormap
world_map <- ggplot() +
  geom_sf(data = map_data, aes(fill = welfare_cat, color = is_afcfta), 
          size = ifelse(map_data$is_afcfta, 0.5, 0.1)) +
  geom_text(data = map_data_labels[map_data_labels$show_on_world,], 
            aes(X, Y, label = iso_a3, alpha = is_afcfta),
            size = ifelse(map_data_labels$is_afcfta[map_data_labels$show_on_world], 2.5, 2),
            fontface = ifelse(map_data_labels$is_afcfta[map_data_labels$show_on_world], "bold", "plain"),
            check_overlap = TRUE) +
  scale_alpha_manual(values = c("FALSE" = 0.7, "TRUE" = 1), guide = "none") +
  scale_fill_manual(values = c(
    "<-0.5%" = "#d73027",       # Deep red for large negative effects
    "-0.5% to -0.1%" = "#f46d43", # Medium red for moderate negative effects
    "-0.1% to 0%" = "#fdae61",   # Light orange for small negative effects
    "0% to 0.1%" = "#ffffbf",    # Light yellow for minimal effects
    "0.1% to 0.5%" = "#abd9e9",  # Light blue for small positive effects
    "0.5% to 1%" = "#74add1",    # Medium blue for moderate positive effects
    "1% to 2%" = "#4575b4",      # Deep blue for large positive effects
    ">2%" = "#313695",           # Very deep blue for very large positive effects
    "No data" = "grey80"         # Grey for missing data
  ), 
  name = "Welfare change\n(% change in real GDP)") +
  scale_color_manual(values = c("FALSE" = NA, "TRUE" = "black"), 
                    guide = "none") +
  labs(title = "Welfare Effects of AfCFTA Agreement on Global Economies",
       subtitle = "% Change in Real GDP - AfCFTA member countries highlighted with black border",
       caption = "Source: Author's calculations using the Structural Gravity Database") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "right",
    panel.grid = element_blank()
  )

# Save the world map plot
ggsave("output/figures/welfare_effects_world_map_struc.pdf", world_map, width = 12, height = 8, dpi = 300)
ggsave("output/figures/welfare_effects_world_map_struc.png", world_map, width = 12, height = 8, dpi = 300)

# Create an alternative version with viridis color palette (colorblind-friendly)
world_map_viridis <- ggplot() +
  geom_sf(data = map_data, 
          aes(fill = welfare_change, color = is_afcfta), 
          size = ifelse(map_data$is_afcfta, 0.5, 0.1)) +
  geom_text(data = map_data_labels[map_data_labels$show_on_world,], 
            aes(X, Y, label = iso_a3, alpha = is_afcfta),
            size = ifelse(map_data_labels$is_afcfta[map_data_labels$show_on_world], 2.5, 2),
            fontface = ifelse(map_data_labels$is_afcfta[map_data_labels$show_on_world], "bold", "plain"),
            color = "white", # White text for better visibility on viridis background
            check_overlap = TRUE) +
  scale_alpha_manual(values = c("FALSE" = 0.7, "TRUE" = 1), guide = "none") +
  scale_fill_viridis(
    option = "plasma",
    name = "Welfare change\n(% change in real GDP)",
    na.value = "grey80"
  ) +
  scale_color_manual(
    values = c("FALSE" = NA, "TRUE" = "black"), 
    guide = "none"
  ) +
  labs(
    title = "Welfare Effects of AfCFTA Agreement on Global Economies",
    subtitle = "% Change in Real GDP - AfCFTA member countries highlighted with black border",
    caption = "Source: Author's calculations using the Structural Gravity Database"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "right",
    panel.grid = element_blank()
  )

# Save the alternative world map
ggsave("output/figures/welfare_effects_world_map_struc_viridis.pdf", 
       world_map_viridis, width = 12, height = 8, dpi = 300)
ggsave("output/figures/welfare_effects_world_map_struc_viridis.png", 
       world_map_viridis, width = 12, height = 8, dpi = 300)

# Create a map focusing on Africa
# For Africa map, create subset of data for African countries and nearby regions
africa_map_data <- map_data[map_data$continent == "Africa" | 
                           map_data$iso_a3 %in% c("SAU", "IRN", "IRQ", "TUR", "ESP", "ITA", "GRC", "FRA", "PRT"),]

# Prepare label data for Africa map - show all African countries
africa_labels <- map_data_labels[map_data_labels$iso_a3 %in% africa_map_data$iso_a3,]

africa_map <- ggplot() +
  geom_sf(data = africa_map_data, 
          aes(fill = welfare_cat, color = is_afcfta), 
          size = ifelse(africa_map_data$is_afcfta, 0.7, 0.1)) +
  geom_text(data = africa_labels, 
            aes(X, Y, label = iso_a3),
            size = 2.8,
            fontface = ifelse(africa_labels$is_afcfta, "bold", "plain"),
            check_overlap = TRUE) +
  scale_fill_manual(values = c(
    "<-0.5%" = "#d73027",       # Deep red for large negative effects
    "-0.5% to -0.1%" = "#f46d43", # Medium red for moderate negative effects
    "-0.1% to 0%" = "#fdae61",   # Light orange for small negative effects
    "0% to 0.1%" = "#ffffbf",    # Light yellow for minimal effects
    "0.1% to 0.5%" = "#abd9e9",  # Light blue for small positive effects
    "0.5% to 1%" = "#74add1",    # Medium blue for moderate positive effects
    "1% to 2%" = "#4575b4",      # Deep blue for large positive effects
    ">2%" = "#313695",           # Very deep blue for very large positive effects
    "No data" = "grey80"         # Grey for missing data
  ), 
  name = "Welfare change\n(% change in real GDP)") +
  scale_color_manual(values = c("FALSE" = NA, "TRUE" = "black"), 
                    guide = "none") +
  labs(title = "Welfare Effects of AfCFTA Agreement on African Economies",
       subtitle = "% Change in Real GDP - AfCFTA member countries highlighted with black border",
       caption = "Source: Author's calculations using the Structural Gravity Database") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "right",
    panel.grid = element_blank()
  )

# Save the Africa-focused map
ggsave("output/figures/welfare_effects_africa_map_struc.pdf", africa_map, width = 12, height = 8, dpi = 300)
ggsave("output/figures/welfare_effects_africa_map_struc.png", africa_map, width = 12, height = 8, dpi = 300)

# Create a bar chart of the top 20 countries with largest welfare effects (positive and negative)
top10_positive <- welfare_data[order(-welfare_change)][1:10]
top10_negative <- welfare_data[order(welfare_change)][1:10]
top20 <- rbind(top10_positive, top10_negative)
top20$is_afcfta <- top20$iso_a3 %in% afcfta_countries

# Add region information for better visualization
top20$region <- ifelse(top20$iso_a3 %in% afcfta_countries, "AfCFTA Member", "Non-AfCFTA Country")

# Create bar chart
welfare_bar <- ggplot(top20, aes(x = reorder(iso_a3, welfare_change), y = welfare_change, fill = region)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("AfCFTA Member" = "#4575b4", "Non-AfCFTA Country" = "#d73027")) +
  labs(
    title = "Top 20 Countries with Largest Welfare Effects from AfCFTA",
    subtitle = "% Change in Real GDP (10 highest gains and 10 highest losses)",
    x = "Country",
    y = "Welfare change (% change in real GDP)",
    fill = "Region"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom"
  )

# Save the bar chart
ggsave("output/figures/welfare_effects_top20_struc.pdf", welfare_bar, width = 10, height = 8, dpi = 300)
ggsave("output/figures/welfare_effects_top20_struc.png", welfare_bar, width = 10, height = 8, dpi = 300)

cat("\nAnalysis completed. Results saved to output folder.\n")
cat("World map of welfare effects saved as:\n")
cat("  - 'welfare_effects_world_map_struc.pdf/png' (standard color scheme)\n")
cat("  - 'welfare_effects_world_map_struc_viridis.pdf/png' (alternative color scheme)\n")
cat("  - 'welfare_effects_africa_map_struc.pdf/png' (Africa-focused map)\n")
cat("  - 'welfare_effects_top20_struc.pdf/png' (Top 20 countries bar chart)\n")

cat("\nSummary of key findings:\n")
cat("  - AfCFTA member countries see an average welfare change of: ", 
    round(mean(afcfta_welfare$welfare_change, na.rm = TRUE), 4), "%\n")
cat("  - Non-AfCFTA countries see an average welfare change of: ", 
    round(mean(non_afcfta_welfare$welfare_change, na.rm = TRUE), 4), "%\n")
cat("  - Country with highest welfare gain: ", 
    welfare_data[order(-welfare_change)][1, iso_a3], " (", 
    round(welfare_data[order(-welfare_change)][1, welfare_change], 4), "%)\n")
cat("  - Country with highest welfare loss: ", 
    welfare_data[order(welfare_change)][1, iso_a3], " (", 
    round(welfare_data[order(welfare_change)][1, welfare_change], 4), "%)\n")

##################################################################
### End of Script
##################################################################