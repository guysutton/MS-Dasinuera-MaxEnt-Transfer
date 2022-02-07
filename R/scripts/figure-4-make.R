# -----------------------------------------------------------------------------
# Script: - Reproduce Figure 4
#
# AUTHOR: Guy F. Sutton
# AFFILIATION: Centre for Biological Control, Rhodes University, South Africa
# DATE MODIFIED: 03/02/2021
# CONTACT: g.sutton@ru.ac.za
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Session setup
# -----------------------------------------------------------------------------

# Load required packages
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  tidyverse,
  spocc,
  ggmap,
  maptools,
  maps,
  ggplot2,
  scrubr,
  mapr,
  tidyr,
  stringr,
  rnaturalearth,
  rnaturalearthdata,
  rlang,
  sf,
  ggspatial,
  raster,
  here,
  tidytext,
  ggtext,
  tidymodels,
  themis,
  ranger,
  yardstick,
  vip,
  DHARMa,
  rms
)

# ----------------------------------------------------------------------------
# Evaluate performance of default settings model: 
# ----------------------------------------------------------------------------

# Import Maxent prediction raster 
clim_map <- raster::raster("./models/rasters/raster_default.asc")

# Import invaded-range GPS points 
inv_data <- readr::read_csv("./data/data_clean/dasi_rubi_invaded_range_gps_clean.csv")
head(inv_data)

# Keep just the columns of interest
inv_clean <- inv_data %>%
  dplyr::select(lon = lon_converted,
                lat = lat_converted,
                established) %>%
  mutate(lat = lat * -1)
inv_clean

# Keep only the sites with establishment data
inv_clean <- inv_clean %>%
  tidyr::drop_na(established)

# Extract the maxent scores for each site
clim_data <- raster::extract(clim_map,
                             inv_clean[, 1:2]) %>%
  as.data.frame()

# Bind the maxent scores to the site data
clim_data <- dplyr::bind_cols(inv_clean, clim_data) %>%
  dplyr::rename(clim_match = 4) %>%
  tidyr::drop_na(clim_match) %>%
  dplyr::mutate(establsihed = as.numeric(established))
head(clim_data)

# Run logistic regression
# - Regress insect establishment (1/0) as a function of 
#   the MaxEnt climate matching score 

# Have to store our data as an 'rms' object for the model to extract predictions later. 
dd <- datadist(clim_data)
options(datadist = "dd")

# Fit model
mod_lrm_default <- lrm(established ~ clim_match,
                       data = clim_data,
                       x = FALSE,
                       y = FALSE)

# Calculate uncertainty
set.seed(48572)
boot_vals_m1 <- boot_val(data = clim_data,  
                         formula = as.formula("established ~ clim_match"))
boot_vals_m1 <- as.data.frame(boot_vals_m1, row.names = NULL)
head(boot_vals_m1)

# Extract AUC data
AUC_data_default <- boot_vals_m1 %>%
  dplyr::slice(2) %>%
  tidyr::pivot_longer(cols = 1:500,
                      names_to = "iteration",
                      values_to = "AUC_scores") %>%
  dplyr::mutate(iteration = 1:500)
head(AUC_data_default)

# ----------------------------------------------------------------------------
# Evaluate performance of AUCtest model: 
# ----------------------------------------------------------------------------

# Import Maxent prediction raster 
clim_map <- raster::raster("./models/rasters/raster_AUCtest.asc")

# Import invaded-range GPS points 
inv_data <- readr::read_csv("./data/data_clean/dasi_rubi_invaded_range_gps_clean.csv")
head(inv_data)

# Keep just the columns of interest
inv_clean <- inv_data %>%
  dplyr::select(lon = lon_converted,
                lat = lat_converted,
                established) %>%
  mutate(lat = lat * -1)
inv_clean

# Keep only the sites with establishment data
inv_clean <- inv_clean %>%
  tidyr::drop_na(established)

# Extract the maxent scores for each site
clim_data <- raster::extract(clim_map,
                             inv_clean[, 1:2]) %>%
  as.data.frame()

# Bind the maxent scores to the site data
clim_data <- dplyr::bind_cols(inv_clean, clim_data) %>%
  dplyr::rename(clim_match = 4) %>%
  tidyr::drop_na(clim_match) %>%
  dplyr::mutate(establsihed = as.numeric(established))
head(clim_data)

# Run logistic regression
# - Regress insect establishment (1/0) as a function of 
#   the MaxEnt climate matching score 

# Have to store our data as an 'rms' object for the model to extract predictions later. 
dd <- datadist(clim_data)
options(datadist = "dd")

# Fit model
mod_lrm_AUCtest <- lrm(established ~ clim_match,
                       data = clim_data,
                       x = FALSE,
                       y = FALSE)

# Calculate uncertainty
set.seed(48572)
boot_vals_m1 <- boot_val(data = clim_data,  
                         formula = as.formula("established ~ clim_match"))
boot_vals_m1 <- as.data.frame(boot_vals_m1, row.names = NULL)
head(boot_vals_m1)

# Extract AUC data
AUC_data_AUCtest <- boot_vals_m1 %>%
  dplyr::slice(2) %>%
  tidyr::pivot_longer(cols = 1:500,
                      names_to = "iteration",
                      values_to = "AUC_scores") %>%
  dplyr::mutate(iteration = 1:500)
head(AUC_data_AUCtest)

# ----------------------------------------------------------------------------
# Evaluate performance of AUCdiff model:
# ----------------------------------------------------------------------------

# Import Maxent prediction raster
clim_map <- raster::raster("./models/rasters/raster_AUCdiff.asc")

# Import invaded-range GPS points
inv_data <-
  readr::read_csv("./data/data_clean/dasi_rubi_invaded_range_gps_clean.csv")
head(inv_data)

# Keep just the columns of interest
inv_clean <- inv_data %>%
  dplyr::select(lon = lon_converted,
                lat = lat_converted,
                established) %>%
  mutate(lat = lat * -1)
inv_clean

# Keep only the sites with establishment data
inv_clean <- inv_clean %>%
  tidyr::drop_na(established)

# Extract the maxent scores for each site
clim_data <- raster::extract(clim_map,
                             inv_clean[, 1:2]) %>%
  as.data.frame()

# Bind the maxent scores to the site data
clim_data <- dplyr::bind_cols(inv_clean, clim_data) %>%
  dplyr::rename(clim_match = 4) %>%
  tidyr::drop_na(clim_match) %>%
  dplyr::mutate(establsihed = as.numeric(established))
head(clim_data)

# Run logistic regression
# - Regress insect establishment (1/0) as a function of
#   the MaxEnt climate matching score

# Have to store our data as an 'rms' object for the model to extract predictions later.
dd <- datadist(clim_data)
options(datadist = "dd")

# Fit model
mod_lrm_AUCdiff <- lrm(established ~ clim_match,
                       data = clim_data,
                       x = FALSE,
                       y = FALSE)

# Calculate uncertainty
set.seed(48572)
boot_vals_m1 <- boot_val(data = clim_data,
                         formula = as.formula("established ~ clim_match"))
boot_vals_m1 <- as.data.frame(boot_vals_m1, row.names = NULL)
head(boot_vals_m1)

# Extract AUC data
AUC_data_AUCdiff <- boot_vals_m1 %>%
  dplyr::slice(2) %>%
  tidyr::pivot_longer(cols = 1:500,
                      names_to = "iteration",
                      values_to = "AUC_scores") %>%
  dplyr::mutate(iteration = 1:500)
head(AUC_data_AUCdiff)

# ----------------------------------------------------------------------------
# Evaluate performance of OR10 model:
# ----------------------------------------------------------------------------

# Import Maxent prediction raster
clim_map <- raster::raster("./models/rasters/raster_OR10.asc")

# Import invaded-range GPS points
inv_data <-
  readr::read_csv("./data/data_clean/dasi_rubi_invaded_range_gps_clean.csv")
head(inv_data)

# Keep just the columns of interest
inv_clean <- inv_data %>%
  dplyr::select(lon = lon_converted,
                lat = lat_converted,
                established) %>%
  mutate(lat = lat * -1)
inv_clean

# Keep only the sites with establishment data
inv_clean <- inv_clean %>%
  tidyr::drop_na(established)

# Extract the maxent scores for each site
clim_data <- raster::extract(clim_map,
                             inv_clean[, 1:2]) %>%
  as.data.frame()

# Bind the maxent scores to the site data
clim_data <- dplyr::bind_cols(inv_clean, clim_data) %>%
  dplyr::rename(clim_match = 4) %>%
  tidyr::drop_na(clim_match) %>%
  dplyr::mutate(establsihed = as.numeric(established))
head(clim_data)

# Run logistic regression
# - Regress insect establishment (1/0) as a function of
#   the MaxEnt climate matching score

# Have to store our data as an 'rms' object for the model to extract predictions later.
dd <- datadist(clim_data)
options(datadist = "dd")

# Fit model
mod_lrm_OR10 <- lrm(established ~ clim_match,
                    data = clim_data,
                    x = FALSE,
                    y = FALSE)

# Calculate uncertainty
set.seed(48572)
boot_vals_m1 <- boot_val(data = clim_data,
                         formula = as.formula("established ~ clim_match"))
boot_vals_m1 <- as.data.frame(boot_vals_m1, row.names = NULL)
head(boot_vals_m1)

# Extract AUC data
AUC_data_OR10 <- boot_vals_m1 %>%
  dplyr::slice(2) %>%
  tidyr::pivot_longer(cols = 1:500,
                      names_to = "iteration",
                      values_to = "AUC_scores") %>%
  dplyr::mutate(iteration = 1:500)
head(AUC_data_OR10)

# ----------------------------------------------------------------------------
# Evaluate performance of AICc model:
# ----------------------------------------------------------------------------

# Import Maxent prediction raster
clim_map <- raster::raster("./models/rasters/raster_AICc.asc")

# Import invaded-range GPS points
inv_data <-
  readr::read_csv("./data/data_clean/dasi_rubi_invaded_range_gps_clean.csv")
head(inv_data)

# Keep just the columns of interest
inv_clean <- inv_data %>%
  dplyr::select(lon = lon_converted,
                lat = lat_converted,
                established) %>%
  mutate(lat = lat * -1)
inv_clean

# Keep only the sites with establishment data
inv_clean <- inv_clean %>%
  tidyr::drop_na(established)

# Extract the maxent scores for each site
clim_data <- raster::extract(clim_map,
                             inv_clean[, 1:2]) %>%
  as.data.frame()

# Bind the maxent scores to the site data
clim_data <- dplyr::bind_cols(inv_clean, clim_data) %>%
  dplyr::rename(clim_match = 4) %>%
  tidyr::drop_na(clim_match) %>%
  dplyr::mutate(establsihed = as.numeric(established))
head(clim_data)

# Run logistic regression
# - Regress insect establishment (1/0) as a function of
#   the MaxEnt climate matching score

# Have to store our data as an 'rms' object for the model to extract predictions later.
dd <- datadist(clim_data)
options(datadist = "dd")

# Fit model
mod_lrm_AICc <- lrm(established ~ clim_match,
                    data = clim_data,
                    x = FALSE,
                    y = FALSE)

# Calculate uncertainty
set.seed(48572)
boot_vals_m1 <- boot_val(data = clim_data,
                         formula = as.formula("established ~ clim_match"))
boot_vals_m1 <- as.data.frame(boot_vals_m1, row.names = NULL)
head(boot_vals_m1)

# Extract AUC data
AUC_data_AICc <- boot_vals_m1 %>%
  dplyr::slice(2) %>%
  tidyr::pivot_longer(cols = 1:500,
                      names_to = "iteration",
                      values_to = "AUC_scores") %>%
  dplyr::mutate(iteration = 1:500)
head(AUC_data_AICc)


# ----------------------------------------------------------------------------
# Make Figure 4:
# ----------------------------------------------------------------------------

# Combine all AUC datasets 
auc_all <- dplyr::bind_rows(AUC_data_default,
                            AUC_data_AUCtest,
                            AUC_data_AUCdiff,
                            AUC_data_OR10,
                            AUC_data_AICc,
                            .id = "id") %>%
  dplyr::mutate(id = dplyr::case_when(
    id == 1 ~ "Default",
    id == 2 ~ "AUCtest",
    id == 3 ~ "AUCdiff",
    id == 4 ~ "OR10",
    id == 5 ~ "AICc")) %>%
  dplyr::mutate(id = as.factor(id))
head(auc_all)

# Plot
auc_all %>%
  dplyr::mutate(id = dplyr::case_when(
    id == "Default" ~ "Default MaxEnt settings",
    id == "AUCtest" ~ "AUC<sub>*test*</sub>",
    id == "AUCdiff" ~ "AUC<sub>*diff*</sub>",
    id == "OR10"    ~ "OR<sub>10</sub>",
    id == "AICc"    ~ "*AIC<sub>c</sub>*"
  )) %>%
  dplyr::mutate(id = forcats::fct_relevel(id, 
                                          "Default MaxEnt settings", 
                                          "AUC<sub>*test*</sub>",
                                          "AUC<sub>*diff*</sub>",
                                          "OR<sub>10</sub>", 
                                          "*AIC<sub>c</sub>*")) %>%
  ggplot(data = ., aes(x = id,
                       y = AUC_scores,
                       group = id)) +
  geom_boxplot(aes(fill = id)) +
  geom_hline(yintercept = 0.75, linetype = "dashed") +
  scale_fill_manual(values = c("grey50",
                               "grey60",
                               "grey70",
                               "grey80",
                               "grey90")) +
  labs(x = "Model settings",
       y = "AUC (out-of-sample)") +
  scale_y_continuous(limits = c(0.4, 1.0),
                     breaks = seq(0.4, 1.0, 0.1)) +
  # Use ggtext to specify that the x-axis text and legend text should be interpreted as markdown text
  theme(
    axis.text.x = element_markdown()
  )

# Save
ggsave("./ms_body/ms_figs/fig_4_auc.png",
       dpi = 600,
       height = 4, 
       width = 6)
