# -----------------------------------------------------------------------------
# Script: 04 - MESS analysis  
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
  dismo,
  raster,
  here,
  corrplot,
  Hmisc,
  patchwork,
  ecospat,
  gridSVG,
  gridExtra,
  grid,
  ENMeval,
  spThin,
  viridis,
  viridisLite,
  mapdata,
  maptools,
  scales,
  geosphere,
  rgdal,
  ggtext,
  rJava,
  rgeos,
  sp,
  sf,
  ggspatial,
  ecospat,
  rnaturalearth,
  rnaturalearthdata,
  rms
)

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

# Set the theme for the maps
theme_opts <- list(
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_rect(fill = 'white', colour = NA),
    plot.background = element_rect(),
    axis.line = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    plot.title = element_text(colour = "black"),
    panel.border = element_rect(fill = NA),
    legend.key = element_blank()
  )
)


# -----------------------------------------------------------------------------
# Load WORLDCLIM layers:
#
# The standard 19 WORLDCLIM layers are imported here to be used as covariates.
# -----------------------------------------------------------------------------

# Load environmental layers
files <-
  list.files(
    path = here::here("./data/environmental_layers/"),
    pattern = 'bil',
    full.names = TRUE
  )

# Check raster layers imported correctly by plotting first layer
mat <-
  raster::raster(here::here("./data/environmental_layers/bio1.bil"))
plot(mat)

# Create a raster stack of all environmental variables
predictors <-
  raster::stack(list.files(
    here::here("./data/environmental_layers/"),
    full.names = TRUE,
    pattern = '.bil'
  ))

# These are the column names of env predictors which were not colinear
climPred <- c("bio1", "bio2", "bio3", "bio6", "bio8",
              "bio12", "bio14", "bio19")

# -----------------------------------------------------------------------------
# Calculate MESS map:
# -----------------------------------------------------------------------------

# Get map of South Africa to project our model over
world <- rnaturalearth::ne_countries(scale = "medium",
                                     returnclass = "sf") %>%
  dplyr::filter(name == "South Africa") %>%
  st_as_sf() %>%
  st_transform(., 4326)

# Subset environmental vars to just layers used to make models
mess_pred <- raster::subset(predictors, climPred)
plot(mess_pred)

# Create a new raster stack for climate in Africa
mess_pred <- raster::crop(mess_pred, world)

# Import the thinned GPS records
species <-
  readr::read_csv(here::here("./data/data_clean/dasi_rubi_native_thinned.csv"))

# Combine thinned records with environmental data
speciesEnv <- base::data.frame(raster::extract(predictors,
                                               cbind(species$lon,
                                                     species$lat)))

# Combine thinned records with environmental data
speciesSwd <- cbind(species, speciesEnv)

# Create points for MESS
reference_points <- speciesSwd[, climPred]

# Select (in order, the vars used to build the model)
reference_points <- reference_points %>%
  base::as.data.frame() %>%
  tidyr::drop_na() %>%
  dplyr::select(climPred)

# Run MESS analysis
mss <- dismo::mess(x = mess_pred,
                   v = reference_points,
                   full = FALSE)

# Change the MESS(-) regions into '0' and MESS(+) into '1'
r.mess.mask <- mss>0
plot(r.mess.mask)

# MEss map
clim_map <- r.mess.mask

# Load invaded-range GPS points 
inv_data <- readr::read_csv("./data/data_clean/dasi_rubi_invaded_range_gps_clean.csv")

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
nrow(inv_clean)

# Extract the MESS 0/1 scores for each site 
clim_data <- raster::extract(clim_map,
                             inv_clean[,1:2])
clim_data 
clim_data <- as.data.frame(clim_data)

# Bind the MESS 0/1 scores to the site data 
clim_data <- cbind(inv_clean, clim_data)
clim_data

# Drop the NA values 
clim_data <- clim_data %>%
  dplyr::select(lon,
                lat,
                established,
                mess_match = clim_data) %>%
  tidyr::drop_na(mess_match)
clim_data
clim_mess <- clim_data

# --------------------------------------------------------------------
# Run GLM:
# --------------------------------------------------------------------

# Have to store our data as an 'rms' object for the model to extract predictions later. 
dd <- datadist(clim_data)
options(datadist = "dd")
head(dd)

# Refit model
mod_lrm_mess <- lrm(established ~ mess_match,
                    data = clim_data,
                    x = TRUE,
                    y = TRUE)
print(mod_lrm_mess)


