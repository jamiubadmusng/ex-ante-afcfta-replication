################################################################################
### R Code for the industry-level analysis of the ex-ante analysis of AfCFTA ###
###################### by Jamiu Olamilekan Badmus ##############################
################################################################################

################################################################################
### Industry-level analysis with full trade data using PPML --------------------
################################################################################

# ---- clear environment -------
rm(list=ls())

# --- Packages ---
library(dplyr)     # data wrangling
library(fixest)    # fepois (PPML with FEs)
library(purrr)     # map / loops
library(broom)     # tidy(), conf.int
library(ggplot2)   # plotting
library(readr)     # write_csv
library(tidyr)     # pivot_longer


# --- Paths ---
base_dir  <- "C:/Users/muham/AfCFTA_UNU_CRIS/01_input"
grav_file <- file.path(base_dir, "dgd_2_1.rds")   # <-- .rds gravity file

# --- AfCFTA members (excl. Eritrea) ---
afc_countries <- c(
  "DZA","AGO","BEN","BWA","BFA","BDI","CPV","CMR","CAF","TCD","COM",
  "COG","CIV","COD","DJI","EGY","GNQ","SWZ","ETH","GAB","GMB","GHA",
  "GIN","GNB","KEN","LSO","LBR","LBY","MDG","MWI","MLI","MRT","MUS",
  "MAR","MOZ","NAM","NER","NGA","RWA","STP","SEN","SYC","SLE","SOM",
  "ZAF","SSD","SDN","TZA","TGO","TUN","UGA","ZMB","ZWE"
)

# --- Data prep function (merge + time-varying border vars) ---
prepare_data <- function(itpde_file) {
  
  # read ITPDE industry data (.rds), then rename to match gravity file
  itpde <- readRDS(itpde_file) %>%
    rename(exporter = exporter_iso3,
           importer = importer_iso3)
  
  # read gravity covariates (.rds) with exporter/importer columns
  grav <- readRDS(grav_file)
  
  # merge
  df <- itpde %>%
    left_join(grav, by = c("exporter", "importer", "year"))
  
  # construct time-varying border variables (as provided)
  df <- df %>%
    mutate(
      # indicator dummies
      afcfta_brdr = ifelse(exporter %in% afc_countries &
                             importer %in% afc_countries &
                             exporter != importer, 1, 0),
      
      afcfta_row_brdr = ifelse(exporter %in% afc_countries &
                                 !(importer %in% afc_countries) &
                                 exporter != importer, 1, 0),
      
      row_to_row = ifelse(!(exporter %in% afc_countries) &
                            !(importer %in% afc_countries) &
                            exporter != importer, 1, 0)
    ) %>%
    mutate(
      # make them time-varying via year-indexed factor codes (your original approach)
      afcfta_brdr     = as.integer(factor(paste(afcfta_brdr, year, sep = "_"))),
      afcfta_row_brdr = as.integer(factor(paste(afcfta_row_brdr, year, sep = "_"))),
      row_to_row      = as.integer(factor(paste(row_to_row, year, sep = "_"))),
      
      # ids for FE and clustering
      pair_id  = paste(exporter, importer, sep = "_"),
      exp_year = paste(exporter, year, sep = "_"),
      imp_year = paste(importer, year, sep = "_")
    )
  
  return(df)
}

# --- Estimation function ---
estimate_model <- function(df, industry_id, industry_descr) {
  m <- fepois(
    trade ~ cntg + col + dist + lang + rta + row_to_row +
      afcfta_brdr + afcfta_row_brdr | exp_year + imp_year,
    data = df,
    vcov = ~ pair_id
  )
  
  # grab only afcfta_brdr
  tidy(m, conf.int = TRUE) %>%
    filter(term == "afcfta_brdr") %>%
    mutate(industry_id = industry_id,
           industry_descr = industry_descr) %>%
    select(industry_id, industry_descr, term, estimate, std.error, conf.low, conf.high, p.value)
}

# --- Loop over 170 industries ---
itpde_files <- file.path(base_dir, sprintf("itpder2_%d.rds", 1:170))

results <- map_dfr(seq_along(itpde_files), function(i) {
  df <- prepare_data(itpde_files[i])
  # Get industry description from the first row (assuming it's the same for all rows in each file)
  industry_descr <- df$industry_descr[1]
  estimate_model(df, i, industry_descr)
})


# --- Classify industries by sector ---
results_with_sector <- results %>%
  filter(!industry_id %in% c(161, 167)) %>%  # Remove outlier industries
  mutate(sector = case_when(
    industry_id >= 1 & industry_id <= 28 ~ "Agriculture",
    industry_id >= 29 & industry_id <= 35 ~ "Mining and Energy",
    industry_id >= 36 & industry_id <= 156 ~ "Manufacturing",
    industry_id >= 157 & industry_id <= 170 ~ "Services",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(sector))  # Remove any industries not in defined ranges

# --- Save results (afcfta_brdr only) ---
out_csv <- file.path("C:/Users/muham/AfCFTA_UNU_CRIS/02_output/csv/table_a_3.csv")
write_csv(results_with_sector, out_csv)

# --- Function to create rank plot by sector ---
create_sector_plot <- function(data, sector_name) {
  # Filter for specific sector
  sector_data <- data %>%
    filter(sector == sector_name) %>%
    mutate(across(c(estimate, conf.low, conf.high), as.numeric)) %>%
    filter(!is.na(estimate), !is.na(conf.low), !is.na(conf.high)) %>%
    arrange(estimate) %>%
    mutate(rank = row_number())
  
  # Convert to long format for plotting
  plot_long <- sector_data %>%
    select(rank, estimate, conf.low, conf.high) %>%
    pivot_longer(cols = c(estimate, conf.low, conf.high),
                 names_to = "series", values_to = "value") %>%
    mutate(series = factor(series,
                           levels = c("estimate", "conf.low", "conf.high"),
                           labels = c("Estimate", "Lower 95 CL", "Upper 95 CL")))
  
  # Create plot
  p <- ggplot(plot_long, aes(x = rank, y = value, color = series, linetype = series)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("Estimate" = "blue",
                                  "Lower 95 CL" = "red",
                                  "Upper 95 CL" = "green")) +
    scale_linetype_manual(values = c("Estimate" = "solid",
                                     "Lower 95 CL" = "dashed",
                                     "Upper 95 CL" = "dashed")) +
    labs(x = "Rank of Estimate in Terms of Size",
         y = "AfCFTA Border Effect (PPML Estimate)",
         title = paste("AfCFTA Border Effects:", sector_name),
         color = NULL, linetype = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "horizontal",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

# --- Create and save plots for each sector ---
sectors <- c("Agriculture", "Mining and Energy", "Manufacturing", "Services")

for (sector in sectors) {
  p <- create_sector_plot(results_with_sector, sector)
  
  # Create filename (replace spaces with underscores)
  filename <- paste0("figure_3_", tolower(gsub(" ", "_", sector)), ".pdf")
  out_pdf <- file.path("C:/Users/muham/AfCFTA_UNU_CRIS/02_output/figures", filename)
  
  ggsave(out_pdf, plot = p, width = 8, height = 5)
  cat("Saved plot for", sector, "to", filename, "\n")
}

# --- OPTIONAL: Create combined plot with facets ---
# Prepare data for all sectors
all_sectors_plot <- results_with_sector %>%
  mutate(across(c(estimate, conf.low, conf.high), as.numeric)) %>%
  filter(!is.na(estimate), !is.na(conf.low), !is.na(conf.high)) %>%
  group_by(sector) %>%
  arrange(estimate, .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Convert to long format
all_sectors_long <- all_sectors_plot %>%
  select(sector, rank, estimate, conf.low, conf.high) %>%
  pivot_longer(cols = c(estimate, conf.low, conf.high),
               names_to = "series", values_to = "value") %>%
  mutate(series = factor(series,
                         levels = c("estimate", "conf.low", "conf.high"),
                         labels = c("Estimate", "Lower 95 CL", "Upper 95 CL")))

# Create faceted plot
p_combined <- ggplot(all_sectors_long, aes(x = rank, y = value, color = series, linetype = series)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ sector, scales = "free", ncol = 2) +
  scale_color_manual(values = c("Estimate" = "blue",
                                "Lower 95 CL" = "red",
                                "Upper 95 CL" = "green")) +
  scale_linetype_manual(values = c("Estimate" = "solid",
                                   "Lower 95 CL" = "dashed",
                                   "Upper 95 CL" = "dashed")) +
  labs(x = "Rank of Estimate in Terms of Size",
       y = "AfCFTA Border Effect (PPML Estimate)",
       color = NULL, linetype = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        strip.text = element_text(face = "bold"))

# Save combined plot
out_pdf_combined <- file.path("C:/Users/muham/AfCFTA_UNU_CRIS/02_output/figures/figure_3.pdf")
ggsave(out_pdf_combined, plot = p_combined, width = 10, height = 8)
cat("Saved combined faceted plot to figure_3.pdf\n")


################################################################################
### Industry-level analysis with only trade > 0 using PPML ---------------------
################################################################################

# ---- clear environment -------
rm(list=ls())

# --- Packages ---
library(dplyr)     # data wrangling
library(fixest)    # fepois (PPML with FEs)
library(purrr)     # map / loops
library(broom)     # tidy(), conf.int
library(ggplot2)   # plotting
library(readr)     # write_csv
library(tidyr)     # pivot_longer


# --- Paths ---
base_dir  <- "C:/Users/muham/AfCFTA_UNU_CRIS/01_input"
grav_file <- file.path(base_dir, "dgd_2_1.rds")   # <-- .rds gravity file

# --- AfCFTA members (excl. Eritrea) ---
afc_countries <- c(
  "DZA","AGO","BEN","BWA","BFA","BDI","CPV","CMR","CAF","TCD","COM",
  "COG","CIV","COD","DJI","EGY","GNQ","SWZ","ETH","GAB","GMB","GHA",
  "GIN","GNB","KEN","LSO","LBR","LBY","MDG","MWI","MLI","MRT","MUS",
  "MAR","MOZ","NAM","NER","NGA","RWA","STP","SEN","SYC","SLE","SOM",
  "ZAF","SSD","SDN","TZA","TGO","TUN","UGA","ZMB","ZWE"
)

# --- Data prep function (merge + time-varying border vars) ---
prepare_data <- function(itpde_file) {
  
  # read ITPDE industry data (.rds), then rename to match gravity file
  itpde <- readRDS(itpde_file) %>%
    rename(exporter = exporter_iso3,
           importer = importer_iso3)
  
  # read gravity covariates (.rds) with exporter/importer columns
  grav <- readRDS(grav_file)
  
  # merge
  df <- itpde %>%
    left_join(grav, by = c("exporter", "importer", "year"))
  
  # construct time-varying border variables (as provided)
  df <- df %>%
    mutate(
      # indicator dummies
      afcfta_brdr = ifelse(exporter %in% afc_countries &
                             importer %in% afc_countries &
                             exporter != importer, 1, 0),
      
      afcfta_row_brdr = ifelse(exporter %in% afc_countries &
                                 !(importer %in% afc_countries) &
                                 exporter != importer, 1, 0),
      
      row_to_row = ifelse(!(exporter %in% afc_countries) &
                            !(importer %in% afc_countries) &
                            exporter != importer, 1, 0)
    ) %>%
    mutate(
      # make them time-varying via year-indexed factor codes (your original approach)
      afcfta_brdr     = as.integer(factor(paste(afcfta_brdr, year, sep = "_"))),
      afcfta_row_brdr = as.integer(factor(paste(afcfta_row_brdr, year, sep = "_"))),
      row_to_row      = as.integer(factor(paste(row_to_row, year, sep = "_"))),
      
      # ids for FE and clustering
      pair_id  = paste(exporter, importer, sep = "_"),
      exp_year = paste(exporter, year, sep = "_"),
      imp_year = paste(importer, year, sep = "_")
    )
  
  return(df)
}

# --- Estimation function ---
estimate_model <- function(df, industry_id, industry_descr) {
  m <- fepois(
    trade ~ cntg + col + dist + lang + rta + row_to_row +
      afcfta_brdr + afcfta_row_brdr | exp_year + imp_year,
    data = df %>% filter(trade > 0),
    vcov = ~ pair_id
  )
  
  # grab only afcfta_brdr
  tidy(m, conf.int = TRUE) %>%
    filter(term == "afcfta_brdr") %>%
    mutate(industry_id = industry_id,
           industry_descr = industry_descr) %>%
    select(industry_id, industry_descr, term, estimate, std.error, conf.low, conf.high, p.value)
}

# --- Loop over 170 industries ---
itpde_files <- file.path(base_dir, sprintf("itpder2_%d.rds", 1:170))

results <- map_dfr(seq_along(itpde_files), function(i) {
  df <- prepare_data(itpde_files[i])
  # Get industry description from the first row (assuming it's the same for all rows in each file)
  industry_descr <- df$industry_descr[1]
  estimate_model(df, i, industry_descr)
})


# --- Classify industries by sector ---
results_with_sector <- results %>%
  filter(!industry_id %in% c(17, 18, 30, 35)) %>%  # Remove outlier industries
  mutate(sector = case_when(
    industry_id >= 1 & industry_id <= 28 ~ "Agriculture",
    industry_id >= 29 & industry_id <= 35 ~ "Mining and Energy",
    industry_id >= 36 & industry_id <= 156 ~ "Manufacturing",
    industry_id >= 157 & industry_id <= 170 ~ "Services",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(sector))  # Remove any industries not in defined ranges


# --- Save results (afcfta_brdr only) ---
out_csv <- file.path("C:/Users/muham/AfCFTA_UNU_CRIS/02_output/csv/table_a_5.csv")
write_csv(results_with_sector, out_csv)

# --- Function to create rank plot by sector ---
create_sector_plot <- function(data, sector_name) {
  # Filter for specific sector
  sector_data <- data %>%
    filter(sector == sector_name) %>%
    mutate(across(c(estimate, conf.low, conf.high), as.numeric)) %>%
    filter(!is.na(estimate), !is.na(conf.low), !is.na(conf.high)) %>%
    arrange(estimate) %>%
    mutate(rank = row_number())
  
  # Convert to long format for plotting
  plot_long <- sector_data %>%
    select(rank, estimate, conf.low, conf.high) %>%
    pivot_longer(cols = c(estimate, conf.low, conf.high),
                 names_to = "series", values_to = "value") %>%
    mutate(series = factor(series,
                           levels = c("estimate", "conf.low", "conf.high"),
                           labels = c("Estimate", "Lower 95 CL", "Upper 95 CL")))
  
  # Create plot
  p <- ggplot(plot_long, aes(x = rank, y = value, color = series, linetype = series)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("Estimate" = "blue",
                                  "Lower 95 CL" = "red",
                                  "Upper 95 CL" = "green")) +
    scale_linetype_manual(values = c("Estimate" = "solid",
                                     "Lower 95 CL" = "dashed",
                                     "Upper 95 CL" = "dashed")) +
    labs(x = "Rank of Estimate in Terms of Size",
         y = "AfCFTA Border Effect (PPML Estimate)",
         title = paste("AfCFTA Border Effects:", sector_name),
         color = NULL, linetype = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.box = "horizontal",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

# --- Create and save plots for each sector ---
sectors <- c("Agriculture", "Mining and Energy", "Manufacturing", "Services")

for (sector in sectors) {
  p <- create_sector_plot(results_with_sector, sector)
  
  # Create filename (replace spaces with underscores)
  filename <- paste0("figure_a1_", tolower(gsub(" ", "_", sector)), ".pdf")
  out_pdf <- file.path("C:/Users/muham/AfCFTA_UNU_CRIS/02_output/figures", filename)
  
  ggsave(out_pdf, plot = p, width = 8, height = 5)
  cat("Saved plot for", sector, "to", filename, "\n")
}

# --- OPTIONAL: Create combined plot with facets ---
# Prepare data for all sectors
all_sectors_plot <- results_with_sector %>%
  mutate(across(c(estimate, conf.low, conf.high), as.numeric)) %>%
  filter(!is.na(estimate), !is.na(conf.low), !is.na(conf.high)) %>%
  group_by(sector) %>%
  arrange(estimate, .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

# Convert to long format
all_sectors_long <- all_sectors_plot %>%
  select(sector, rank, estimate, conf.low, conf.high) %>%
  pivot_longer(cols = c(estimate, conf.low, conf.high),
               names_to = "series", values_to = "value") %>%
  mutate(series = factor(series,
                         levels = c("estimate", "conf.low", "conf.high"),
                         labels = c("Estimate", "Lower 95 CL", "Upper 95 CL")))

# Create faceted plot
p_combined <- ggplot(all_sectors_long, aes(x = rank, y = value, color = series, linetype = series)) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ sector, scales = "free", ncol = 2) +
  scale_color_manual(values = c("Estimate" = "blue",
                                "Lower 95 CL" = "red",
                                "Upper 95 CL" = "green")) +
  scale_linetype_manual(values = c("Estimate" = "solid",
                                   "Lower 95 CL" = "dashed",
                                   "Upper 95 CL" = "dashed")) +
  labs(x = "Rank of Estimate in Terms of Size",
       y = "AfCFTA Border Effect (PPML Estimate)",
       color = NULL, linetype = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        strip.text = element_text(face = "bold"))

# Save combined plot
out_pdf_combined <- file.path("C:/Users/muham/AfCFTA_UNU_CRIS/02_output/figures/figure_a1.pdf")
ggsave(out_pdf_combined, plot = p_combined, width = 10, height = 8)
cat("Saved combined faceted plot to figure_a1.pdf\n")


################################################################################
### End of script --------------------------------------------------------------
################################################################################