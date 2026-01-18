################################################################################
######## Data preparation for GE analysis of ex-ante impacts of AfCFTA #########
######################## by Jamiu Olamilekan Badmus ############################
################################################################################

# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Set the working directory ----------------------------------------------------
setwd("C:/Users/muham/AfCFTA_UNU_CRIS")


# Load the necessary libraries -------------------------------------------------
library(tidyverse)
library(readr)
library(haven)

# Load the datasets ------------------------------------------------------------
itpde <- readRDS("01_input/itpder2.rds") # Trade Dataset
dgd <- readRDS("01_input/dgd_2_1.rds") # DGD Dataset

# AfCFTA members excluding Eritrea
afc_countries <- c("DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "COM", 
                   "COG", "CIV", "COD", "DJI", "EGY", "GNQ", "SWZ", "ETH", "GAB", "GMB", "GHA",
                   "GIN", "GNB", "KEN", "LSO", "LBR", "LBY", "MDG", "MWI", "MLI", "MRT", "MUS",
                   "MAR", "MOZ", "NAM", "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM",
                   "ZAF", "SSD", "SDN", "TZA", "TGO", "TUN", "UGA", "ZMB", "ZWE")

# Define wanted_countries to include all AfCFTA countries (excluding Eritrea) 
wanted_countries <- c(
  # Include all AfCFTA countries (already excludes Eritrea)
  afc_countries,
  
  # EU member countries
  "AUT", "BEL", "BGR", "HRV", "CYP", "CZE", "DNK", "EST", "FIN", "FRA", "DEU", "GBR", "GRC", "HUN",
  "IRL", "ITA", "LVA", "LTU", "LUX", "MLT", "NLD", "POL", "PRT", "ROU", "SVK", "SVN", "ESP", "SWE",
  
  # NAFTA countries
  "USA", "CAN", "MEX",
  
  # BRICS (excluding South Africa which is in AfCFTA)
  "BRA", "RUS", "IND", "CHN",
  
  # Other G20 members
  "ARG", "AUS", "IDN", "JPN", "KOR", "SAU", "TUR",
  
  # Major Asian economies
  "HKG", "SGP", "TWN", "THA", "MYS", "VNM", "PHL", "PAK", "BGD", "YEM",

  # Middle East
  "ARE", "BHR", "IRN", "IRQ", "ISR", "JOR", "KWT", "LBN", "OMN", "QAT",
  
  # Europe and others
  "CHE", "NOR", "ISL", "NZL", "CHL", "COL", "PER", "VEN", "ECU",
  "KAZ", "UZB", "UKR", "BLR", "GEO", "ARM", "AZE",
  
  # Other regions
  "JAM", "TTO", "BRB", "BHS", "DOM", "CUB", "GTM", "SLV", "HND", "NIC", "CRI", "PAN",
  "FJI", "PNG", "NPL", "LKA", "MMR", "LAO", "KHM", "BRN", "TJK", "KGZ", "TKM",
  "ALB", "MKD", "MNE", "SRB", "BIH", "MDA"
)

# Remove duplicates (in case any country appears in both lists)
wanted_countries <- unique(wanted_countries)

# Subset ITPDE data for the year 2019 and filter for wanted countries
itpde_2019 <- itpde %>%
  filter(year == 2019,
         exporter %in% wanted_countries,
         importer %in% wanted_countries)

# Function to create balanced squared dataset with proper domestic trade imputation
# This works on ITPDE data BEFORE merging with DGD to avoid missing gravity variables
create_balanced_squared_itpde <- function(df, countries_list, sector_name) {
  cat("Processing", sector_name, "sector...\n")
  
  # Step 1: Filter for the specific sector
  sector_df <- df %>%
    filter(broad_sector == sector_name)
  
  # Step 2: Calculate export-to-internal-sales ratios for this sector
  cat("Computing export-to-internal-sales ratios for", sector_name, "...\n")
  
  # Calculate total exports for each country (excluding domestic trade)
  total_exports <- sector_df %>%
    filter(exporter != importer) %>%
    group_by(exporter) %>%
    summarise(total_exports = sum(trade, na.rm = TRUE), .groups = 'drop')
  
  # Get domestic trade values where available
  domestic_trade <- sector_df %>%
    filter(exporter == importer & trade > 0 & !is.na(trade)) %>%
    select(exporter, domestic_trade = trade)
  
  # Calculate export-to-domestic ratios for countries with both export and domestic data
  export_domestic_ratios <- total_exports %>%
    inner_join(domestic_trade, by = "exporter") %>%
    mutate(export_to_domestic_ratio = total_exports / domestic_trade) %>%
    filter(is.finite(export_to_domestic_ratio) & export_to_domestic_ratio > 0)
  
  # Calculate median ratio for this sector as benchmark
  if (nrow(export_domestic_ratios) > 0) {
    median_ratio <- median(export_domestic_ratios$export_to_domestic_ratio, na.rm = TRUE)
    cat("  Median export-to-domestic ratio:", round(median_ratio, 3), "\n")
  } else {
    # Fallback: use a reasonable default ratio
    median_ratio <- 0.3  # Default assumption: exports = 30% of domestic trade
    cat("No ratio data available, using default ratio:", median_ratio, "\n")
  }
  
  # Step 3: Create all possible combinations of exporter-importer pairs
  all_pairs <- expand.grid(
    exporter = countries_list,
    importer = countries_list,
    stringsAsFactors = FALSE
  )
  
  # Add year 2019 and sector information
  all_pairs$year <- 2019
  all_pairs$broad_sector <- sector_name
  
  # Step 4: Merge with existing sector data
  balanced_df <- all_pairs %>%
    left_join(sector_df, by = c("exporter", "importer", "year", "broad_sector"))
  
  # Step 5: Calculate total exports for all countries in wanted_countries
  all_total_exports <- balanced_df %>%
    filter(exporter != importer & !is.na(trade)) %>%
    group_by(exporter) %>%
    summarise(total_exports = sum(trade, na.rm = TRUE), .groups = 'drop')
  
  # Step 6: Impute missing domestic trade using the ratio approach
  balanced_df <- balanced_df %>%
    left_join(all_total_exports, by = "exporter") %>%
    mutate(
      # For domestic trade (exporter == importer)
      imputed_domestic = ifelse(
        exporter == importer & total_exports > 0,
        total_exports / median_ratio,  # domestic = exports / ratio
        NA
      ),
      
      # Apply imputation logic
      trade = case_when(
        # Keep existing non-missing, non-zero trade values
        !is.na(trade) & trade > 0 ~ trade,
        
        # For domestic trade: use imputed value if available, otherwise set small positive value
        exporter == importer & !is.na(imputed_domestic) & imputed_domestic > 0 ~ imputed_domestic,
        exporter == importer ~ pmax(1, total_exports * 0.1, na.rm = TRUE),  # Fallback: 10% of exports, min 1
        
        # For international trade: set to zero if missing
        TRUE ~ 0
      ),
      
      # Ensure other necessary variables are present
      year = 2019,
      broad_sector = sector_name
    ) %>%
    select(-total_exports, -imputed_domestic) %>%  # Remove temporary variables
    arrange(exporter, importer)
  
  # Step 7: Final check - ensure no domestic trade is zero or missing
  domestic_check <- balanced_df %>%
    filter(exporter == importer) %>%
    summarise(
      min_domestic = min(trade, na.rm = TRUE),
      zero_domestic = sum(trade == 0, na.rm = TRUE),
      na_domestic = sum(is.na(trade))
    )
  
  cat("  Domestic trade check - Min:", round(domestic_check$min_domestic, 2),
      "Zero count:", domestic_check$zero_domestic,
      "NA count:", domestic_check$na_domestic, "\n")
  
  return(balanced_df)
}

# Create balanced squared datasets for each sector from ITPDE data
cat("Creating balanced squared datasets from ITPDE data...\n")

# Create balanced datasets for each sector
afcfta_2019_agri_balanced <- create_balanced_squared_itpde(itpde_2019, wanted_countries, "Agriculture")
afcfta_2019_manu_balanced <- create_balanced_squared_itpde(itpde_2019, wanted_countries, "Manufacturing")
afcfta_2019_mine_balanced <- create_balanced_squared_itpde(itpde_2019, wanted_countries, "Mining and Energy")
afcfta_2019_serv_balanced <- create_balanced_squared_itpde(itpde_2019, wanted_countries, "Services")

# Now merge each balanced dataset with DGD data to get gravity variables
cat("\nMerging with DGD data to add gravity variables...\n")

# Function to merge with DGD and add derived variables
add_gravity_variables <- function(df, sector_name) {
  cat("Adding gravity variables to", sector_name, "...\n")
  
  # Merge with DGD data
  df_with_gravity <- df %>%
    left_join(
      dgd %>% select(year, exporter, importer, cntg, col, lang, dist, rta),
      by = c("year", "exporter", "importer")
    ) %>%
    mutate(
      # Create border variables
      afcfta_brdr = ifelse(exporter %in% afc_countries &
                            importer %in% afc_countries &
                            exporter != importer, 1, 0),
      
      afcfta_row_brdr = ifelse(exporter %in% afc_countries &
                                !(importer %in% afc_countries) &
                                exporter != importer, 1, 0),
      
      row_to_row = ifelse(!(exporter %in% afc_countries) & 
                           !(importer %in% afc_countries) &
                           exporter != importer, 1, 0),
      
      # RTA among AfCFTA countries
      rta_afcfta = ifelse(exporter %in% afc_countries &
                           importer %in% afc_countries &
                           exporter != importer, 1, 0),
      
      # Create ID variables
      exp_year = as.integer(factor(paste(exporter, year, sep = "_"))),
      imp_year = as.integer(factor(paste(importer, year, sep = "_"))),
      pair_id = as.integer(factor(paste(exporter, importer, sep = "_")))
    ) %>%
    arrange(exporter, importer)
  
  return(df_with_gravity)
}

# Add gravity variables to each balanced dataset
afcfta_2019_agri_final <- add_gravity_variables(afcfta_2019_agri_balanced, "Agriculture")
afcfta_2019_manu_final <- add_gravity_variables(afcfta_2019_manu_balanced, "Manufacturing")
afcfta_2019_mine_final <- add_gravity_variables(afcfta_2019_mine_balanced, "Mining and Energy")
afcfta_2019_serv_final <- add_gravity_variables(afcfta_2019_serv_balanced, "Services")
# Save balanced datasets with all variables
cat("\nSaving final balanced datasets...\n")
write_dta(afcfta_2019_agri_final, "01_input/afcfta_2019_agri_balanced.dta")
write_dta(afcfta_2019_manu_final, "01_input/afcfta_2019_manu_balanced.dta")
write_dta(afcfta_2019_mine_final, "01_input/afcfta_2019_mine_balanced.dta")
write_dta(afcfta_2019_serv_final, "01_input/afcfta_2019_serv_balanced.dta")

cat("Balanced datasets saved successfully!\n")

# Verification: Check that each dataset is truly balanced and squared
cat("\nVerification:\n")
datasets <- list(
  "AGRI" = afcfta_2019_agri_final,
  "MANU" = afcfta_2019_manu_final,
  "MINE" = afcfta_2019_mine_final,
  "SERV" = afcfta_2019_serv_final
)

for (sector_name in names(datasets)) {
  df <- datasets[[sector_name]]
  
  n_exporters <- length(unique(df$exporter))
  n_importers <- length(unique(df$importer))
  n_countries <- length(wanted_countries)
  
  # Check domestic trade specifically
  domestic_stats <- df %>%
    filter(exporter == importer) %>%
    summarise(
      n_domestic = n(),
      min_domestic = min(trade, na.rm = TRUE),
      max_domestic = max(trade, na.rm = TRUE),
      mean_domestic = mean(trade, na.rm = TRUE),
      zero_domestic = sum(trade == 0, na.rm = TRUE),
      na_domestic = sum(is.na(trade))
    )
  
  # Check for missing gravity variables
  gravity_missing <- df %>%
    summarise(
      missing_cntg = sum(is.na(cntg)),
      missing_dist = sum(is.na(dist)),
      missing_rta = sum(is.na(rta))
    )
  
  cat(sprintf("%s: %d exporters, %d importers, %d total countries\n", 
              sector_name, n_exporters, n_importers, n_countries))
  cat(sprintf("  Expected rows: %d, Actual rows: %d\n", n_countries^2, nrow(df)))
  cat(sprintf("  Domestic trade: %d rows, Min: %.2f, Max: %.2f, Mean: %.2f\n",
              domestic_stats$n_domestic, domestic_stats$min_domestic, 
              domestic_stats$max_domestic, domestic_stats$mean_domestic))
  cat(sprintf("  Zero/NA domestic: %d zeros, %d NAs\n", 
              domestic_stats$zero_domestic, domestic_stats$na_domestic))
  cat(sprintf("  Missing gravity vars: cntg=%d, dist=%d, rta=%d\n",
              gravity_missing$missing_cntg, gravity_missing$missing_dist, gravity_missing$missing_rta))
  cat(sprintf("  International trade rows: %d\n", sum(df$exporter != df$importer)))
}



##### Robustness check ---------------------------------------------------------

# Load structural gravity database for robustness check
cat("\n=== ROBUSTNESS CHECK: Using Structural Gravity Database ===\n")

# 1. Load the structural gravity data
struc <- readRDS("01_input/struc.rds")
cat("Loaded structural gravity database\n")

# 2. Subset the data to 2016
struc_2016 <- struc %>%
  filter(year == 2016,
         exporter %in% wanted_countries,
         importer %in% wanted_countries)

cat("Filtered structural gravity data for 2016 and wanted countries\n")
cat("Number of observations:", nrow(struc_2016), "\n")

# 3. Create balanced squared dataset with proper domestic trade imputation
# Function for structural gravity data (no broad_sector breakdown)
create_balanced_squared_struc <- function(df, countries_list) {
  cat("Creating balanced squared dataset for structural gravity data...\n")
  
  # Step 1: Calculate export-to-internal-sales ratios
  cat("Computing export-to-internal-sales ratios...\n")
  
  # Calculate total exports for each country (excluding domestic trade)
  total_exports <- df %>%
    filter(exporter != importer) %>%
    group_by(exporter) %>%
    summarise(total_exports = sum(trade, na.rm = TRUE), .groups = 'drop')
  
  # Get domestic trade values where available
  domestic_trade <- df %>%
    filter(exporter == importer & trade > 0 & !is.na(trade)) %>%
    select(exporter, domestic_trade = trade)
  
  # Calculate export-to-domestic ratios for countries with both export and domestic data
  export_domestic_ratios <- total_exports %>%
    inner_join(domestic_trade, by = "exporter") %>%
    mutate(export_to_domestic_ratio = total_exports / domestic_trade) %>%
    filter(is.finite(export_to_domestic_ratio) & export_to_domestic_ratio > 0)
  
  # Calculate median ratio as benchmark
  if (nrow(export_domestic_ratios) > 0) {
    median_ratio <- median(export_domestic_ratios$export_to_domestic_ratio, na.rm = TRUE)
    cat("  Median export-to-domestic ratio:", round(median_ratio, 3), "\n")
  } else {
    # Fallback: use a reasonable default ratio for manufacturing
    median_ratio <- 0.25  # Default assumption: exports = 25% of domestic trade for manufacturing
    cat("  No ratio data available, using default ratio:", median_ratio, "\n")
  }
  
  # Step 2: Create all possible combinations of exporter-importer pairs
  all_pairs <- expand.grid(
    exporter = countries_list,
    importer = countries_list,
    stringsAsFactors = FALSE
  )
  
  # Add year 2016
  all_pairs$year <- 2016
  
  # Step 3: Merge with existing data
  balanced_df <- all_pairs %>%
    left_join(df, by = c("exporter", "importer", "year"))
  
  # Step 4: Calculate total exports for all countries in wanted_countries
  all_total_exports <- balanced_df %>%
    filter(exporter != importer & !is.na(trade)) %>%
    group_by(exporter) %>%
    summarise(total_exports = sum(trade, na.rm = TRUE), .groups = 'drop')
  
  # Step 5: Impute missing domestic trade using the ratio approach
  balanced_df <- balanced_df %>%
    left_join(all_total_exports, by = "exporter") %>%
    mutate(
      # For domestic trade (exporter == importer)
      imputed_domestic = ifelse(
        exporter == importer & total_exports > 0,
        total_exports / median_ratio,  # domestic = exports / ratio
        NA
      ),
      
      # Apply imputation logic
      trade = case_when(
        # Keep existing non-missing, non-zero trade values
        !is.na(trade) & trade > 0 ~ trade,
        
        # For domestic trade: use imputed value if available, otherwise set small positive value
        exporter == importer & !is.na(imputed_domestic) & imputed_domestic > 0 ~ imputed_domestic,
        exporter == importer ~ pmax(1, total_exports * 0.1, na.rm = TRUE),  # Fallback: 10% of exports, min 1
        
        # For international trade: set to zero if missing
        TRUE ~ 0
      ),
      
      # Ensure year is set
      year = 2016
    ) %>%
    select(-total_exports, -imputed_domestic) %>%  # Remove temporary variables
    arrange(exporter, importer)
  
  # Step 6: Final check - ensure no domestic trade is zero or missing
  domestic_check <- balanced_df %>%
    filter(exporter == importer) %>%
    summarise(
      min_domestic = min(trade, na.rm = TRUE),
      zero_domestic = sum(trade == 0, na.rm = TRUE),
      na_domestic = sum(is.na(trade))
    )
  
  cat("  Domestic trade check - Min:", round(domestic_check$min_domestic, 2),
      "Zero count:", domestic_check$zero_domestic,
      "NA count:", domestic_check$na_domestic, "\n")
  
  return(balanced_df)
}

# Create balanced squared dataset
struc_2016_balanced <- create_balanced_squared_struc(struc_2016, wanted_countries)

# 4. Merge with standard gravity covariates dataset and add derived variables
cat("\nMerging with robustness gravity dataset...\n")

# Load robustness gravity dataset
dgd_rob <- readRDS("01_input/dgd_2_1_rob.rds")

# Merge with gravity variables and add derived variables
struc_2016_final <- struc_2016_balanced %>%
  left_join(
    dgd_rob %>% select(year, exporter, importer, cntg, col, lang, dist, rta),
    by = c("year", "exporter", "importer")
  ) %>%
  mutate(
    # Create border variables
    afcfta_brdr = ifelse(exporter %in% afc_countries &
                          importer %in% afc_countries &
                          exporter != importer, 1, 0),
    
    afcfta_row_brdr = ifelse(exporter %in% afc_countries &
                              !(importer %in% afc_countries) &
                              exporter != importer, 1, 0),
    
    row_to_row = ifelse(!(exporter %in% afc_countries) & 
                         !(importer %in% afc_countries) &
                         exporter != importer, 1, 0),
    
    # RTA among AfCFTA countries
    rta_afcfta = ifelse(exporter %in% afc_countries &
                         importer %in% afc_countries &
                         exporter != importer, 1, 0),
    
    # Create ID variables
    exp_year = as.integer(factor(paste(exporter, year, sep = "_"))),
    imp_year = as.integer(factor(paste(importer, year, sep = "_"))),
    pair_id = as.integer(factor(paste(exporter, importer, sep = "_")))
  ) %>%
  arrange(exporter, importer)

# 5. Save as "struc_2016_balanced.dta"
cat("Saving robustness dataset...\n")
write_dta(struc_2016_final, "01_input/struc_2016_balanced.dta")

# Verification for robustness dataset
cat("\nRobustness Dataset Verification:\n")

n_exporters <- length(unique(struc_2016_final$exporter))
n_importers <- length(unique(struc_2016_final$importer))
n_countries <- length(wanted_countries)

# Check domestic trade specifically
domestic_stats <- struc_2016_final %>%
  filter(exporter == importer) %>%
  summarise(
    n_domestic = n(),
    min_domestic = min(trade, na.rm = TRUE),
    max_domestic = max(trade, na.rm = TRUE),
    mean_domestic = mean(trade, na.rm = TRUE),
    zero_domestic = sum(trade == 0, na.rm = TRUE),
    na_domestic = sum(is.na(trade))
  )

# Check for missing gravity variables
gravity_missing <- struc_2016_final %>%
  summarise(
    missing_cntg = sum(is.na(cntg)),
    missing_dist = sum(is.na(dist)),
    missing_rta = sum(is.na(rta))
  )

cat(sprintf("STRUC 2016: %d exporters, %d importers, %d total countries\n", 
            n_exporters, n_importers, n_countries))
cat(sprintf("  Expected rows: %d, Actual rows: %d\n", n_countries^2, nrow(struc_2016_final)))
cat(sprintf("  Domestic trade: %d rows, Min: %.2f, Max: %.2f, Mean: %.2f\n",
            domestic_stats$n_domestic, domestic_stats$min_domestic, 
            domestic_stats$max_domestic, domestic_stats$mean_domestic))
cat(sprintf("  Zero/NA domestic: %d zeros, %d NAs\n", 
            domestic_stats$zero_domestic, domestic_stats$na_domestic))
cat(sprintf("  Missing gravity vars: cntg=%d, dist=%d, rta=%d\n",
            gravity_missing$missing_cntg, gravity_missing$missing_dist, gravity_missing$missing_rta))
cat(sprintf("  International trade rows: %d\n", sum(struc_2016_final$exporter != struc_2016_final$importer)))

cat("Robustness dataset preparation completed successfully!\n")

################################################################################
# End of the code
################################################################################