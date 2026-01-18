################################################################################
### R Code for the sectoral-level analysis of the ex-ante analysis of AfCFTA ###
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


# Load the datasets ------------------------------------------------------------
itpder2 <- readRDS("01_input/itpder2.rds") # Trade Dataset
dgd_2_1 <- readRDS("01_input/dgd_2_1.rds") # DGD Dataset

# Merge datasets
afcfta_data <- left_join(
  itpder2,
  dgd_2_1 %>%
    select(year, exporter, importer, cntg, col, lang, dist, rta),
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
### Table 2 analysis with full trade data using PPML ---------------------------
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


## Total Trade Estimation for AfCFTA --------------------------------------------
tot_obj = fepois(trade ~ cntg + col + dist + lang + rta + row_to_row
                  + afcfta_brdr + afcfta_row_brdr
                  | exp_year + imp_year,
                  data = afcfta_data,
                  vcov = cluster ~ pair_id)
summary(tot_obj)

# Delta method for Total Trade
expr_afcfta <- paste0("(exp(afcfta_brdr) - 1) * 100")
delta_afcfta_tot <- deltaMethod(object = coef(tot_obj), vcov. = vcov(tot_obj), g = expr_afcfta, parameterNames = names(coef(tot_obj)))
print(delta_afcfta_tot)

## Agriculture trade
agri_obj = fepois(trade ~ cntg + col + dist + lang + rta + row_to_row
                  + afcfta_brdr + afcfta_row_brdr
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Agriculture"),
                  vcov = cluster ~ pair_id)
summary(agri_obj)

# Delta method for Agriculture
delta_afcfta_agri <- deltaMethod(object = coef(agri_obj), vcov. = vcov(agri_obj), g = expr_afcfta, parameterNames = names(coef(agri_obj)))
print(delta_afcfta_agri)


## Manufacturing trade
manu_obj = fepois(trade ~ cntg + col + dist + lang + rta + row_to_row
                  + afcfta_brdr + afcfta_row_brdr
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Manufacturing"),
                  vcov = cluster ~ pair_id)
summary(manu_obj)

# Delta method for Manufacturing
delta_afcfta_manu <- deltaMethod(object = coef(manu_obj), vcov. = vcov(manu_obj), g = expr_afcfta, parameterNames = names(coef(manu_obj)))
print(delta_afcfta_manu)

## Mining and Energy Trade
mine_obj = fepois(trade ~ cntg + col + dist + lang + rta + row_to_row
                  + afcfta_brdr + afcfta_row_brdr
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Mining and Energy"),
                  vcov = cluster ~ pair_id)
summary(mine_obj)

# Delta method for Mining and Energy
delta_afcfta_mine <- deltaMethod(object = coef(mine_obj), vcov. = vcov(mine_obj), g = expr_afcfta, parameterNames = names(coef(mine_obj)))
print(delta_afcfta_mine)

# Service trade
serv_obj = fepois(trade ~ cntg + col + dist + lang + rta + row_to_row
                  + afcfta_brdr + afcfta_row_brdr
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(broad_sector == "Services"),
                  vcov = cluster ~ pair_id)
summary(serv_obj)

# Delta method for Services
delta_afcfta_serv <- deltaMethod(object = coef(serv_obj), vcov. = vcov(serv_obj), g = expr_afcfta, parameterNames = names(coef(serv_obj)))
print(delta_afcfta_serv)

# AfCFTA Latex results ---------------------------------------------------------
table_2 <- huxreg("Tot" = tot_obj,
                     "Agri" = agri_obj,
                     "Manu" = manu_obj,
                     "Mine" = mine_obj,
                     "Serv" = serv_obj,
                     coefs = c("CNTG" = "cntg",
                               "COLN" = "col",
                               "DIST" = "dist",
                               "LANG" = "lang",
                               "RTA" = "rta",
                               "INTL_BRDR" = "row_to_row",
                               "BRDR_AFC" = "afcfta_brdr",
                               "BRDR_AFC_ROW" = "afcfta_row_brdr"),
                     stars = c("***"=0.001, "**"=0.01, "*"=0.05),
                     note = "Notes: Statistics based on the author's calculations. The table presents the gravity estimates of the AfCFTA Border Effects on intra-African Total (Tot), Agriculture (Agri), Manufacturing (Manu), Mining & Energy (Mine), and Service (Serv) trade flows. The estimates are obtained using the PPML estimator with nominal bilateral trade flows (including zero trade flows), two-way (exporter-time and importer-time) fixed effects, and clustered standard errors (reported in parentheses) by country pair. *** p < 0.01, ** p < 0.05, * p < 0.10.") %>%
  insert_row("", "(1)", "(2)", "(3)", "(4)", "(5)", after = 0) %>%
  set_top_border(1, everywhere, 1) %>%
  set_align(1, everywhere, "center") %>%
  set_col_width(c(0.25, 0.15, 0.15, 0.15, 0.15, 0.15)) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Gravity Estimates of the AfCFTA Border Effects on Sectoral Intra-African Trade") %>%
  set_label("table_2")
width(table_2) <- 1
cat(to_latex(table_2), file = "02_output/tables/table_2.tex")



################################################################################
### Table 4 analysis with only trade > 0 using PPML ----------------------------
################################################################################

## Total Trade Estimation for AfCFTA --------------------------------------------
tot_obj = fepois(trade ~ cntg + col + dist + lang + rta + row_to_row
                  + afcfta_brdr + afcfta_row_brdr
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(trade > 0),
                  vcov = cluster ~ pair_id)
summary(tot_obj)

# Delta method for Total Trade
expr_afcfta <- paste0("(exp(afcfta_brdr) - 1) * 100")
delta_afcfta_tot <- deltaMethod(object = coef(tot_obj), vcov. = vcov(tot_obj), g = expr_afcfta, parameterNames = names(coef(tot_obj)))
print(delta_afcfta_tot)

## Agriculture trade
agri_obj = fepois(trade ~ cntg + col + dist + lang + rta + row_to_row
                  + afcfta_brdr + afcfta_row_brdr
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(trade > 0 & broad_sector == "Agriculture"),
                  vcov = cluster ~ pair_id)
summary(agri_obj)

# Delta method for Agriculture
delta_afcfta_agri <- deltaMethod(object = coef(agri_obj), vcov. = vcov(agri_obj), g = expr_afcfta, parameterNames = names(coef(agri_obj)))
print(delta_afcfta_agri)


## Manufacturing trade
manu_obj = fepois(trade ~ cntg + col + dist + lang + rta + row_to_row
                  + afcfta_brdr + afcfta_row_brdr
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(trade > 0 & broad_sector == "Manufacturing"),
                  vcov = cluster ~ pair_id)
summary(manu_obj)

# Delta method for Manufacturing
delta_afcfta_manu <- deltaMethod(object = coef(manu_obj), vcov. = vcov(manu_obj), g = expr_afcfta, parameterNames = names(coef(manu_obj)))
print(delta_afcfta_manu)

## Mining and Energy Trade
mine_obj = fepois(trade ~ cntg + col + dist + lang + rta + row_to_row
                  + afcfta_brdr + afcfta_row_brdr
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(trade > 0 & broad_sector == "Mining and Energy"),
                  vcov = cluster ~ pair_id)
summary(mine_obj)

# Delta method for Mining and Energy
delta_afcfta_mine <- deltaMethod(object = coef(mine_obj), vcov. = vcov(mine_obj), g = expr_afcfta, parameterNames = names(coef(mine_obj)))
print(delta_afcfta_mine)

# Service trade
serv_obj = fepois(trade ~ cntg + col + dist + lang + rta + row_to_row
                  + afcfta_brdr + afcfta_row_brdr
                  | exp_year + imp_year,
                  data = afcfta_data %>% filter(trade > 0 & broad_sector == "Services"),
                  vcov = cluster ~ pair_id)
summary(serv_obj)

# Delta method for Services
delta_afcfta_serv <- deltaMethod(object = coef(serv_obj), vcov. = vcov(serv_obj), g = expr_afcfta, parameterNames = names(coef(serv_obj)))
print(delta_afcfta_serv)

# AfCFTA Latex results ---------------------------------------------------------
table_a_4 <- huxreg("Tot" = tot_obj,
                     "Agri" = agri_obj,
                     "Manu" = manu_obj,
                     "Mine" = mine_obj,
                     "Serv" = serv_obj,
                     coefs = c("CNTG" = "cntg",
                               "COLN" = "col",
                               "DIST" = "dist",
                               "LANG" = "lang",
                               "RTA" = "rta",
                               "INTL_BRDR" = "row_to_row",
                               "BRDR_AFC" = "afcfta_brdr",
                               "BRDR_AFC_ROW" = "afcfta_row_brdr"),
                     stars = c("***"=0.001, "**"=0.01, "*"=0.05),
                     note = "Notes: Statistics based on the author's calculations. The table presents the gravity estimates of the AfCFTA Border Effects on intra-African Total (Tot), Agriculture (Agri), Manufacturing (Manu), Mining & Energy (Mine), and Service (Serv) trade flows. The estimates are obtained using the PPML estimator with nominal bilateral trade flows (excluding zero trade flows), two-way (exporter-time and importer-time) fixed effects, and clustered standard errors (reported in parentheses) by country pair. *** p < 0.01, ** p < 0.05, * p < 0.10.") %>%
  insert_row("", "(1)", "(2)", "(3)", "(4)", "(5)", after = 0) %>%
  set_top_border(1, everywhere, 1) %>%
  set_align(1, everywhere, "center") %>%
  set_col_width(c(0.25, 0.15, 0.15, 0.15, 0.15, 0.15)) %>%
  set_align(everywhere,-1,"center") %>%
  set_caption("Gravity Estimates of the AfCFTA Border Effects on Sectoral Intra-African Trade (Excluding Zero Trade Flows)") %>%
  set_label("table_a_4")
width(table_a_4) <- 1
cat(to_latex(table_a_4), file = "02_output/tables/table_a_4.tex")


################################################################################
### End of script --------------------------------------------------------------
################################################################################