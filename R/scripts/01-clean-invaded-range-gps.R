###########################################################################
# Script: 01 - Clean invaded range records  -------------------------------
#
# AUTHOR: Guy F. Sutton
# AFFILIATION: Centre for Biological Control, Rhodes University, South Africa
# DATE MODIFIED: 03/02/2021
# CONTACT: g.sutton@ru.ac.za
###########################################################################

###########################################################################
# Setup -------------------------------------------------------------------
###########################################################################

# Load required packages
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(tidyverse,
               tidyr,
               DHARMa,
               here,
               readxl,
               glmmTMB)

# Change ggplot theme
theme_set(
  theme_classic() +
    theme(
      panel.border = element_rect(colour = "black",
                                  fill = NA),
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(margin = unit(c(2, 0, 0, 0),
                                                "mm")),
      axis.title.y = element_text(margin = unit(c(0, 4, 0, 0),
                                                "mm")),
      legend.position = "none"
    )
)

# Load source functions
# - Load in function to process GPS records
source(here::here("./R/functions/function_convert_GPS.R"))

###########################################################################
# Load data ---------------------------------------------------------------
###########################################################################

# Load invaded-range GPS points
raw_data <-
  readxl::read_xlsx(
    here::here(
      "./data/data_raw/dasineura_rubiformis_invaded_range_database_cleaned.xlsx"
    )
  )
head(raw_data)

###########################################################################
# Process data ------------------------------------------------------------
###########################################################################

# Remove NA values in latitude column
raw_data <- raw_data %>%
  tidyr::drop_na(lat2)
raw_data

# Need to convert the decimal minutes into decimal degrees
data_proc <- raw_data %>%
  dplyr::mutate(
    lat_converted = dg2dec(
      varb = raw_data$lat2,
      Dg = "°",
      Min = "'"
    ),
    lon_converted = dg2dec(
      varb = raw_data$lon2,
      Dg = "°",
      Min = "'"
    )
  )
head(data_proc)

# How many invaded-range points do we have?
nrow(data_proc)

# Round nicely (4 decimal points)
data_proc <- data_proc %>%
  dplyr::mutate(
    lat_converted = base::round(lat_converted, digits = 4),
    lon_converted = base::round(lon_converted, digits = 4)
  )

# Write the processed data to a .csv file to use later
readr::write_csv(
  data_proc,
  here::here("./data/data_clean/dasi_rubi_invaded_range_gps_clean.csv")
)
