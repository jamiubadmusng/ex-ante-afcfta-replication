################################################################################
####### R Code for robustness checks of the ex-ante analysis of AfCFTA #########
###################### by Jamiu Olamilekan Badmus ##############################
################################################################################

# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Set the working directory ----------------------------------------------------
setwd("C:/Users/muham/AfCFTA_UNU_CRIS")


# Load the necessary libraries -------------------------------------------------
library(tidyverse)
library(here)
library(fixest)
library(flextable)
library(huxtable)
library(readr)
library(msm)
library(car)
library(haven)


# Load the datasets ------------------------------------------------------------
struc <- readRDS("01_input/struc.rds") # Trade Dataset
dgd_2_1_rob <- readRDS("01_input/dgd_2_1_rob.rds") # DGD Dataset

#limit struc data to 2000-2016
struc <- struc %>%
  filter(year >= 2000 & year <= 2016)

# Merge datasets
afcfta_data <- left_join(
  struc,
  dgd_2_1_rob %>%
    select(year, exporter, importer, cntg, col, lang, dist, rta, wto, pta, eia),
  by = c("year", "exporter", "importer")
)


# construct exporter-time, importer-time, and asymmetric pair fixed effects
afcfta_data <- afcfta_data %>%
  mutate(
    exp_year = as.integer(factor(paste(exporter, year, sep = "_"))),
    imp_year = as.integer(factor(paste(importer, year, sep = "_"))),
    pair_id = as.integer(factor(paste(exporter, importer, sep = "_")))
  )


################################################################################
### Robustness check models ----------------------------------------------------
################################################################################

# AfCFTA members excluding Eritrea
afc_countries <- c("DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "COM", 
                   "COG", "CIV", "COD", "DJI", "EGY", "GNQ", "SWZ", "ETH", "GAB", "GMB", "GHA",
                   "GIN", "GNB", "KEN", "LSO", "LBR", "LBY", "MDG", "MWI", "MLI", "MRT", "MUS",
                   "MAR", "MOZ", "NAM", "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM",
                   "ZAF", "SSD", "SDN", "TZA", "TGO", "TUN", "UGA", "ZMB", "ZWE")


# create time-varying border variables between AfCFTA countries and ROW
afcfta_data <- afcfta_data %>%
  mutate(
    # indicator for intra-AfCFTA international trade
    afcfta_brdr = ifelse(exporter %in% afc_countries &
                           importer %in% afc_countries &
                           exporter != importer, 1, 0),
    
    # indicator for AfCFTA exports to the ROW
    afcfta_row_brdr = ifelse(exporter %in% afc_countries &
                               !(importer %in% afc_countries) &
                               exporter != importer, 1, 0),
    
    # indicator for ROW exports to ROW
    row_to_row = ifelse(!(exporter %in% afc_countries) &
                          !(importer %in% afc_countries) &
                          exporter != importer, 1, 0)
  ) %>%
  mutate(
    afcfta_brdr = as.integer(factor(paste(afcfta_brdr, year, sep = "_"))),
    afcfta_row_brdr = as.integer(factor(paste(afcfta_row_brdr, year, sep = "_"))),
    row_to_row = as.integer(factor(paste(row_to_row, year, sep = "_")))
  )


## With only rta
rta = fepois(trade ~ cntg + col + dist + lang + rta + row_to_row
             + afcfta_brdr + afcfta_row_brdr
             | exp_year + imp_year,
             data = afcfta_data,
             vcov = cluster ~ pair_id)
summary(rta)

## With rta + wto
wto = fepois(trade ~ cntg + col + dist + lang + rta + wto + row_to_row
             + afcfta_brdr + afcfta_row_brdr
             | exp_year + imp_year,
             data = afcfta_data,
             vcov = cluster ~ pair_id)
summary(wto)

## With rta + wto + pta
pta = fepois(trade ~ cntg + col + dist + lang + rta + wto + pta + row_to_row
             + afcfta_brdr + afcfta_row_brdr
             | exp_year + imp_year,
             data = afcfta_data,
             vcov = cluster ~ pair_id)
summary(pta)


## With rta + wto + pta + eia
eia = fepois(trade ~ cntg + col + dist + lang + rta + wto + pta + eia + row_to_row
             + afcfta_brdr + afcfta_row_brdr
             | exp_year + imp_year,
             data = afcfta_data,
             vcov = cluster ~ pair_id)
summary(eia)



# AfCFTA robustness results  ---------------------------------------------------
table_a_7 <- huxreg("RTA" = rta,
                         "+ WTO" = wto,
                         "+ PTA" = pta,
                         "+ EIA" = eia,
                         coefs = c("CNTG" = "cntg",
                                   "COLN" = "col",
                                   "DIST" = "dist",
                                   "LANG" = "lang",
                                   "RTA" = "rta",
                                   "WTO" = "wto",
                                   "PTA" = "pta",
                                   "EIA" = "eia",
                                   "INTL_BRDR" = "row_to_row",
                                   "BRDR_AFC" = "afcfta_brdr",
                                   "BRDR_AFC_ROW" = "afcfta_row_brdr"),
                         stars = c("***"=0.001, "**"=0.01, "*"=0.05),
                         note = "Notes: Statistics based on the author’s calculations. The table presents the robustness checks of the gravity estimates of the AfCFTA border on intra-African trade under different specifications. The estimates are obtained using the PPML estimator with nominal bilateral trade flows, exporter-time and importer-time fixed effects, and clustered standard errors (reported in parentheses) by country pair. *** p < 0.01, ** p < 0.05, * p < 0.10.") %>%
  insert_row("", "(1)", "(2)", "(3)", "(4)", after = 0) %>%
  set_top_border(1, everywhere, 1) %>%
  set_align(1, everywhere, "center") %>%
  set_col_width(c(0.32, 0.17, 0.17, 0.17, 0.17)) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Robustness Checks for Gravity Estimates of the AfCFTA Border Effects on Sectoral Intra-African Trade") %>%
  set_label("table_a_7")
width(table_a_7) <- 1.1
cat(to_latex(table_a_7), file = "02_output/tables/table_a_7.tex")

# Save workspace --------------------------------------------------------------
save.image("afcfta_robustness_checks.RData")

################################################################################
#### End of script #############################################################
################################################################################
