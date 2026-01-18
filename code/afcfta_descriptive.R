################################################################################
### R Code for the descriptive analysis the ex-ante quantification of AfCFTA ###
###################### by Jamiu Olamilekan Badmus ##############################
################################################################################

# clear the workspace ----------------------------------------------------------
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
library(dplyr)
library(kableExtra)
library(dplyr)
library(ggplot2)
library(scales) 


# Load the datasets ------------------------------------------------------------
itpder2 <- readRDS("01_input/itpder2.rds") # Trade Dataset

dgd1 <- read_csv("01_input/release_2.1_2000_2019/release_2.1_2000_2004.csv") # DGD Dataset 2000-2004
dgd2 <- read_csv("01_input/release_2.1_2000_2019/release_2.1_2005_2009.csv") # DGD Dataset 2005-2009
dgd3 <- read_csv("01_input/release_2.1_2000_2019/release_2.1_2010_2014.csv") # DGD Dataset 2010-2014
dgd4 <- read_csv("01_input/release_2.1_2000_2019/release_2.1_2015_2019.csv") # DGD Dataset 2015-2019

dgd <- rbind(dgd1, dgd2, dgd3, dgd4) # Combine DGD datasets

# subset relevant columns from DGD dataset
dgd_sub <- dgd %>%
  select(year, iso3_o, iso3_d, colony_ever, contiguity, distance, common_language,agree_fta, gdp_pwt_cur_o, gdp_pwt_cur_d)

# rename iso3_o to exporter and iso3_d to importer to match trade dataset
dgd_sub <- dgd_sub %>%
  rename(exporter = iso3_o,
         importer = iso3_d,
         cntg = contiguity,
         col = colony_ever,
         dist = distance,
         lang = common_language,
         rta = agree_fta,
         gdp_exp = gdp_pwt_cur_o,
         gdp_imp = gdp_pwt_cur_d)


# Merge datasets
afcfta_data <- left_join(
  itpder2,
  dgd_sub %>%
    select(year, exporter, importer, cntg, col, lang, dist, rta, gdp_exp, gdp_imp),
  by = c("year", "exporter", "importer")
)

# List of African countries
afc_countries <- c(
  "DZA", "AGO", "BEN", "BWA", "BFA", "BDI", "CPV", "CMR", "CAF", "TCD", "COM", 
  "COG", "CIV", "COD", "DJI", "EGY", "GNQ", "SWZ", "ETH", "GAB", "GMB", "GHA",
  "GIN", "GNB", "KEN", "LSO", "LBR", "LBY", "MDG", "MWI", "MLI", "MRT", "MUS",
  "MAR", "MOZ", "NAM", "NER", "NGA", "RWA", "STP", "SEN", "SYC", "SLE", "SOM",
  "ZAF", "SSD", "SDN", "TZA", "TGO", "TUN", "UGA", "ZMB", "ZWE"
)

# Create intra-Africa trade dataset
intra_africa <- afcfta_data %>%
  filter(exporter %in% afc_countries, importer %in% afc_countries, exporter != importer)



# 1. Trade summary statistics by sector
trade_by_sector <- intra_africa %>%
  group_by(broad_sector) %>%
  summarise(
    Observations = n(),
    Pairs = n_distinct(paste(exporter, importer, sep = "-")),
    Mean = mean(trade, na.rm = TRUE),
    Median = median(trade, na.rm = TRUE),
    SD = sd(trade, na.rm = TRUE),
    Min = min(trade, na.rm = TRUE),
    Max = max(trade, na.rm = TRUE)
  ) %>%
  mutate(across(where(is.numeric), round, 2))

# 2. Overall summary statistics for all other variables
overall_summary <- intra_africa %>%
  summarise(
    # Distance
    Dist_N = sum(!is.na(dist)),
    Dist_Mean = mean(dist, na.rm = TRUE),
    Dist_Median = median(dist, na.rm = TRUE),
    Dist_SD = sd(dist, na.rm = TRUE),
    Dist_Min = min(dist, na.rm = TRUE),
    Dist_Max = max(dist, na.rm = TRUE),
    
    # Contiguity
    Cntg_N = sum(!is.na(cntg)),
    Cntg_Mean = mean(cntg, na.rm = TRUE),
    Cntg_Median = median(cntg, na.rm = TRUE),
    Cntg_SD = sd(cntg, na.rm = TRUE),
    Cntg_Min = min(cntg, na.rm = TRUE),
    Cntg_Max = max(cntg, na.rm = TRUE),
    
    # Colony
    Col_N = sum(!is.na(col)),
    Col_Mean = mean(col, na.rm = TRUE),
    Col_Median = median(col, na.rm = TRUE),
    Col_SD = sd(col, na.rm = TRUE),
    Col_Min = min(col, na.rm = TRUE),
    Col_Max = max(col, na.rm = TRUE),
    
    # Language
    Lang_N = sum(!is.na(lang)),
    Lang_Mean = mean(lang, na.rm = TRUE),
    Lang_Median = median(lang, na.rm = TRUE),
    Lang_SD = sd(lang, na.rm = TRUE),
    Lang_Min = min(lang, na.rm = TRUE),
    Lang_Max = max(lang, na.rm = TRUE),
    
    # RTA
    RTA_N = sum(!is.na(rta)),
    RTA_Mean = mean(rta, na.rm = TRUE),
    RTA_Median = median(rta, na.rm = TRUE),
    RTA_SD = sd(rta, na.rm = TRUE),
    RTA_Min = min(rta, na.rm = TRUE),
    RTA_Max = max(rta, na.rm = TRUE),
    
    # GDP Exporter
    GDP_Exp_N = sum(!is.na(gdp_exp)),
    GDP_Exp_Mean = mean(gdp_exp, na.rm = TRUE),
    GDP_Exp_Median = median(gdp_exp, na.rm = TRUE),
    GDP_Exp_SD = sd(gdp_exp, na.rm = TRUE),
    GDP_Exp_Min = min(gdp_exp, na.rm = TRUE),
    GDP_Exp_Max = max(gdp_exp, na.rm = TRUE),
    
    # GDP Importer
    GDP_Imp_N = sum(!is.na(gdp_imp)),
    GDP_Imp_Mean = mean(gdp_imp, na.rm = TRUE),
    GDP_Imp_Median = median(gdp_imp, na.rm = TRUE),
    GDP_Imp_SD = sd(gdp_imp, na.rm = TRUE),
    GDP_Imp_Min = min(gdp_imp, na.rm = TRUE),
    GDP_Imp_Max = max(gdp_imp, na.rm = TRUE)
  )

# 3. Format overall summary into a table
overall_summary_table <- data.frame(
  Variable = c("Distance", "Contiguity", "Colony", "Common Language", 
               "RTA", "GDP (Exporter)", "GDP (Importer)"),
  Observations = c(overall_summary$Dist_N,
                   overall_summary$Cntg_N,
                   overall_summary$Col_N,
                   overall_summary$Lang_N,
                   overall_summary$RTA_N,
                   overall_summary$GDP_Exp_N,
                   overall_summary$GDP_Imp_N),
  Mean = c(overall_summary$Dist_Mean,
           overall_summary$Cntg_Mean,
           overall_summary$Col_Mean,
           overall_summary$Lang_Mean,
           overall_summary$RTA_Mean,
           overall_summary$GDP_Exp_Mean,
           overall_summary$GDP_Imp_Mean),
  Median = c(overall_summary$Dist_Median,
             overall_summary$Cntg_Median,
             overall_summary$Col_Median,
             overall_summary$Lang_Median,
             overall_summary$RTA_Median,
             overall_summary$GDP_Exp_Median,
             overall_summary$GDP_Imp_Median),
  SD = c(overall_summary$Dist_SD,
         overall_summary$Cntg_SD,
         overall_summary$Col_SD,
         overall_summary$Lang_SD,
         overall_summary$RTA_SD,
         overall_summary$GDP_Exp_SD,
         overall_summary$GDP_Imp_SD),
  Min = c(overall_summary$Dist_Min,
          overall_summary$Cntg_Min,
          overall_summary$Col_Min,
          overall_summary$Lang_Min,
          overall_summary$RTA_Min,
          overall_summary$GDP_Exp_Min,
          overall_summary$GDP_Imp_Min),
  Max = c(overall_summary$Dist_Max,
          overall_summary$Cntg_Max,
          overall_summary$Col_Max,
          overall_summary$Lang_Max,
          overall_summary$RTA_Max,
          overall_summary$GDP_Exp_Max,
          overall_summary$GDP_Imp_Max)
) %>%
  mutate(across(where(is.numeric), round, 2))

# 4. Create LaTeX table for trade by sector
trade_sector_table <- kable(trade_by_sector, 
                            format = "latex", 
                            booktabs = TRUE,
                            col.names = c("Broad Sector", "Observations", "Pairs", 
                                          "Mean", "Median", "Std. Dev.", "Min", "Max"),
                            caption = "Summary Statistics of Trade by Broad Sector",
                            label = "tab:trade_sector") %>%
  kable_styling(latex_options = c("hold_position", "striped"))

# 5. Create LaTeX table for overall variables
overall_var_table <- kable(overall_summary_table, 
                           format = "latex", 
                           booktabs = TRUE,
                           col.names = c("Variable", "Observations", "Mean", 
                                         "Median", "Std. Dev.", "Min", "Max"),
                           caption = "Summary Statistics of Gravity Covariates (Overall Dataset)",
                           label = "tab:overall_vars") %>%
  kable_styling(latex_options = c("hold_position", "striped"))

# 6. Save tables
# Create output directory if it doesn't exist
if (!dir.exists("02_output/tables")) {
  dir.create("02_output/tables", recursive = TRUE)
}

save_kable(trade_sector_table, file = "02_output/tables/table_1a.tex")
save_kable(overall_var_table, file = "02_output/tables/table_1b.tex")


###### Trade and Tariff

## Compute trade/GDP
intra_african <- intra_africa %>%
  mutate(trade = trade / gdp_exp)


# Plot intra-African trade flows and tariff over time by broad sector ----------
# Summarize trade by year and broad sector
intra_africa_year <- intra_african %>%
  group_by(year, broad_sector) %>%
  summarise(total_trade = sum(trade, na.rm = TRUE), .groups = "drop")

# Load intra-Africa_tariff data
intra_africa_tariff <- read_csv("01_input/intra_africa_tariff.csv")

# Filter AHS tariff data
intra_africa_tariff <- intra_africa_tariff %>%
  filter(DutyType == "AHS") %>%
  select('Reporter Name', Product, 'Product Name', 'Partner Name', 'Tariff Year', DutyType, 'Weighted Average')

# Subset the year to 2000 to 2019
intra_africa_tariff <- intra_africa_tariff %>%
  filter(`Tariff Year` >= 2000 & `Tariff Year` <= 2019)

# Prepare the ahs_data for sectoral mapping
intra_africa_tariff <- intra_africa_tariff %>%
  mutate(hs6 = substr(Product, 1, 6),
         hs2 = substr(hs6, 1, 2),
         sector = case_when(
           hs2 %in% sprintf("%02d", 1:24) ~ "Agriculture",
           hs2 %in% c("25", "26", "27") ~ "Mining and Energy",
           hs2 %in% sprintf("%02d", 28:99) ~ "Manufacturing",
           TRUE ~ NA_character_
         ))

# Compute weighted sectoral averages
intra_africa_tariff <- intra_africa_tariff %>%
  filter(sector %in% c("Agriculture", "Manufacturing", "Mining and Energy")) %>%
  group_by(`Tariff Year`, sector) %>%
  summarise(
    avg_ahs_tariff = mean(`Weighted Average`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(year = `Tariff Year`, broad_sector = sector)

# Merge trade and tariff data for sectors with tariff data
intra_africa_trade_tariff <- intra_africa_year %>%
  filter(broad_sector %in% c("Agriculture", "Manufacturing", "Mining and Energy")) %>%
  left_join(intra_africa_tariff, by = c("year", "broad_sector"))


# Function to create individual sector plots with dual y-axes
create_sector_plot <- function(data, sector_name, subplot_label) {
  # Filter data for specific sector
  sector_data <- data %>% filter(broad_sector == sector_name)
  
  # Calculate scaling factor for secondary axis
  trade_max <- max(sector_data$total_trade, na.rm = TRUE)
  tariff_max <- max(sector_data$avg_ahs_tariff, na.rm = TRUE)
  scale_factor <- trade_max / tariff_max
  
  # Create the plot
  ggplot(sector_data, aes(x = year)) +
    # Trade line - primary y-axis (left)
    geom_line(aes(y = total_trade, color = "Export (% of GDP)"), 
              size = 1.2, alpha = 0.8) +
    geom_point(aes(y = total_trade, color = "Export (% of GDP)"), 
               size = 2, alpha = 0.8) +
    
    # Tariff line - secondary y-axis (right)
    geom_line(aes(y = avg_ahs_tariff * scale_factor, color = "Average Tariff"), 
              size = 1.2, linetype = "dashed", alpha = 0.8) +
    geom_point(aes(y = avg_ahs_tariff * scale_factor, color = "Average Tariff"), 
               size = 2, alpha = 0.8) +
    
    # Primary y-axis for trade
    scale_y_continuous(
      name = "Export (% of GDP)",
      labels = comma_format(),
      # Secondary y-axis for tariffs
      sec.axis = sec_axis(~ . / scale_factor, 
                          name = "AHS Weighted Average Tariff (%)",
                          labels = comma_format())
    ) +
    
    # Customize colors
    scale_color_manual(values = c("Export (% of GDP)" = "steelblue", 
                                  "Average Tariff" = "darkred")) +
    
    # Labels
    labs(
      title = paste0("(", subplot_label, ") ", sector_name),
      x = "Year",
      color = NULL
    ) +
    
    # Theme customization
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      axis.title.y.left = element_text(color = "steelblue", size = 10),
      axis.title.y.right = element_text(color = "darkred", size = 10),
      axis.text.y.left = element_text(color = "steelblue"),
      axis.text.y.right = element_text(color = "darkred")
    ) +
    
    scale_x_continuous(breaks = seq(2000, 2019, 4))
}

# Function to create Services plot (trade only)
create_services_plot <- function(data, sector_name, subplot_label) {
  # Filter data for Services sector
  sector_data <- data %>% filter(broad_sector == sector_name)
  
  # Create the plot
  ggplot(sector_data, aes(x = year)) +
    # Trade line only
    geom_line(aes(y = total_trade, color = "Export (% of GDP)"), 
              size = 1.2, alpha = 0.8) +
    geom_point(aes(y = total_trade, color = "Export (% of GDP)"), 
               size = 2, alpha = 0.8) +
    
    # Y-axis for trade only
    scale_y_continuous(
      name = "Export (% of GDP)",
      labels = comma_format()
    ) +
    
    # Customize colors
    scale_color_manual(values = c("Export (% of GDP)" = "steelblue")) +
    
    # Labels
    labs(
      title = paste0("(", subplot_label, ") ", sector_name),
      x = "Year",
      color = NULL
    ) +
    
    # Theme customization
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(color = "steelblue", size = 10),
      axis.text.y = element_text(color = "steelblue")
    ) +
    
    scale_x_continuous(breaks = seq(2000, 2019, 4))
}


# Create individual plots
library(gridExtra)  # Load gridExtra package for grid.arrange

# Plot (a): Agriculture
plot_agri <- create_sector_plot(intra_africa_trade_tariff, "Agriculture", "a")

# Plot (b): Manufacturing  
plot_manu <- create_sector_plot(intra_africa_trade_tariff, "Manufacturing", "b")

# Plot (c): Mining and Energy
plot_mine <- create_sector_plot(intra_africa_trade_tariff, "Mining and Energy", "c")

# Plot (d): Services (trade only)
services_data <- intra_africa_year %>% filter(broad_sector == "Services")
plot_serv <- create_services_plot(intra_africa_year, "Services", "d")

# Save individual plots as PDF
ggsave("02_output/figures/figure_2a.pdf", 
       plot = plot_agri, 
       width = 8, height = 6, dpi = 300)

ggsave("02_output/figures/figure_2b.pdf", 
       plot = plot_manu, 
       width = 8, height = 6, dpi = 300)

ggsave("02_output/figures/figure_2c.pdf", 
       plot = plot_mine, 
       width = 8, height = 6, dpi = 300)

ggsave("02_output/figures/figure_2d.pdf", 
       plot = plot_serv, 
       width = 8, height = 6, dpi = 300)

# Create combined plot with all four subplots
combined_plot <- grid.arrange(plot_agri, plot_manu, plot_mine, plot_serv, 
                              nrow = 2, ncol = 2)

# Save combined plot
ggsave("02_output/figures/figure_2.pdf", 
       plot = combined_plot, 
       width = 16, height = 12, dpi = 300)


###################################################################################
# End of the script ---------------------------------------------------------------
###################################################################################