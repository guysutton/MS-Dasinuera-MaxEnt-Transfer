# -----------------------------------------------------------------------------
# Script: 02 - MaxEnt models - native-range (Australia)  
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
  kuenm,
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
  rnaturalearthdata
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

# Load source functions
# - Load in function to process GPS records
source(here::here("./R/functions/function_convert_GPS.R"))
source(here::here("./R/functions/correlation_wheel.R"))

# Start the MaxEnt java application (after pushing system memory up)
memory.limit(memory.limit()*2^30)
dismo::maxent()

# -----------------------------------------------------------------------------
# Load species occurence records:
# -----------------------------------------------------------------------------

# Import cleaned native-range data
species <-
  readr::read_csv(here::here("./data/data_clean/dasi_rubi_native_range_gps_clean.csv"))

# Keep only columns of interest 
species <- species %>%
  dplyr::select(lat = lat_converted,
                lon = lon_converted,
                state,
                plant_host = acacia_sampled)
species

# Need to make the latitude values into S, not N
species <- species %>%
  dplyr::mutate(lat = lat * -1)

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

# -----------------------------------------------------------------------------
# Remove duplicate records:
#
# In this section, we remove duplicate records and thin occurrences to one gps
# point per raster cell.
# -----------------------------------------------------------------------------

# Remove duplicated data based on latitude and longitude
dups <- base::duplicated(species[c("lat", "lon")])
occ_unique <- species[!dups,]
cat(nrow(species) - nrow(occ_unique), "records are removed")
species <- occ_unique


# -----------------------------------------------------------------------------
# Select background points:
#
# In this section, we will select background points (pseudo-absences) from a cropped 
# area surrounding the gps points.
# Define background area by only Koppen-Geiger zones with at least one GPS record. 
# -----------------------------------------------------------------------------

# Load Koppen-Geiger layer 
eco <- rgdal::readOGR(here::here("./data/shapefiles"), 
                      "WC05_1975H_Koppen", 
                      verbose = FALSE)

# Plot to make sure KG layer imported correctly 
sp::plot(eco)

# Coerce GPS records into SPDF
recordsSpatial <- sp::SpatialPointsDataFrame(
  coords = cbind(species$lon, species$lat),
  data = species,
  proj4string = CRS(
    '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  )
)

# Reproject KG-layer
geo_proj = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
eco = spTransform(eco, geo_proj)

# Plot to check we can overlay points on KG map
sp::plot(eco)
points(species$lon,
       species$lat,
       pch = 21,
       bg = 'mediumseagreen')

# Select KG ecoregions in which there is at least one GPS record
ecoContain <- eco[recordsSpatial, ]

# Plot regions containing at least one record
sp::plot(ecoContain)
sp::plot(eco, 
         add = TRUE, 
         col = 'gray70')
# Fill the ecoregions with a khaki colour if they contain at least 1 GPS point
sp::plot(ecoContain, 
         add = TRUE, 
         col = 'khaki')
# Overlay GPS points 
points(species$lon, 
       species$lat, 
       pch = 21, 
       bg = 'mediumseagreen', 
       cex = 1)

# Define background area by KG regions 
# - Will throw 'initialize' errors when code is run multiple
#   times in the same R session.
# - Don't be alarmed. 
# - Just plot the studyArea raster after the errors to make sure it worked.
# - Plot should be the WORLDclim layers masked to keep only the 
#   khaki zones from the Koppen-Geiger map.
# - This code will crop the WORLDCLIM layers to only the KG ecoregions 
#   containing at least 1 GPS records (i.e. khaki in maps above )
studyArea <- raster::mask(predictors, ecoContain)  
plot(studyArea[[1]]) # See, it does work despite the warning messages 

# Select background points from the masked area
# - We return the background points in spatial points format using 'sp = TRUE` 
bgPoints <- raster::sampleRandom(
  x = studyArea,
  size = 10000,
  na.rm = TRUE,
  sp = TRUE
) 

# Plot background extent with background points and species occurrences
plot(studyArea[[1]])

# Add the background points to the plotted raster
plot(bgPoints,
     add = T,
     col = "blue")

# Convert background points into a data frame
bgPoints <- raster::as.data.frame(bgPoints)
head(bgPoints)

# --------------------------------------------------------------------------
# Reduce multicollinearity between environmental predictors:
#
# In this section, we make sure that we do not produce overfit models 
# due to colinearity between environmental layers. 
# We will use correlation co-efficient (r |>0.70|), and remove any predictors above. 
# --------------------------------------------------------------------------

# Extract WORLDCLIM values at each species GPS point used to make models
speciesEnv <- base::data.frame(raster::extract(predictors,
                                               cbind(species$lon,
                                                     species$lat)))

# Plot correlation wheel for our potential predictor variables 
corrWheel( speciesEnv[ , c("bio1", "bio2", "bio3", "bio4", "bio5", 
                            "bio6", "bio7", "bio8", "bio9", "bio10", "bio11",
                            "bio12", "bio13", "bio14", "bio15", "bio16", 
                            "bio17", "bio18", "bio19")], 
           # Specify the R2 value we want as our lowest limit here
           0.70)

# Subset just our potential predictor variables where R2 < 0.70
redPred <- c("bio1", "bio2", "bio3", "bio6", "bio8", 
              "bio12", "bio14", "bio19")

# Must remove NA values first
corData <- speciesEnv[ , redPred] %>%
  tidyr::drop_na()

# Make correlation matrix
corData <- stats::cor(corData)
head(corData)

# Plot the correlation wheel nicely to visualise colinearity
# - Check that no values are greater than 0.7 (or whatever R2 you choose above)
corrplot::corrplot(
  corData,
  type = "upper",
  method = "number",
  order = "original",
  tl.col = "black",
  tl.srt = 45
)

# -----------------------------------------------------------------------------
# Spatial thinning and spatial filtering:
#
# In this section, we perform spatial filtering to remove records 
# that demonstrate spatial autocorrelation. 
# -----------------------------------------------------------------------------

# Must remove NA values first
corData <- speciesEnv[ , redPred] %>%
  tidyr::drop_na()

# Combine background points with species GPS data 
spatialData <- dplyr::bind_cols(corData,
                                species) 
head(spatialData)

# At what distance does spatial autocorrelation occur?
#spatialCorAnalysis <- ecospat::ecospat.mantel.correlogram(
#  # Data frame with environmental variables
#  dfvar = spatialData,
#  # Columns containing lat/long
#  colxy = 9:10,
#  # Number of random occurrences
#  n = 100,
#  # Climate variables
#  colvar = 1:8,
#  max = 1000,
#  nclass = 100,
#  nperm = 100
#)

# Add a species column to the GPS data
#species <- species %>%
#  dplyr::mutate(species = "D. rubiformis")

# Thin by spatial autocorrelation value
#speciesThinned <- spThin::thin(
#  loc.data = species,
#  lat.col = "lat",
#  long.col = "lon",
#  spec.col = "species",
#  # Km unit of correlogram
#  thin.par = 0.2,
#  reps = 100,
#  max.files = 1,
#  out.dir = here::here("./data/data_clean/")
#)

# Import the thinned GPS records
speciesThinned <-
  readr::read_csv(here::here("./data/data_clean/dasi_rubi_native_thinned.csv"))
head(speciesThinned)

# How many records were removed?
nrow(species)
nrow(speciesThinned)

# Just for ease, make species_thinned = species
species <- speciesThinned

# --------------------------------------------------------------------
# Set up data to run MaxEnt models:
#
# Tell R where to store the MaxEnt outputs, and do final pre-processing
# of the GPS data and raster layers before running models 
# --------------------------------------------------------------------

# Create a new directory to store model output
# Use getwd() to find working directory
dir.create(here::here("./models/"), 
           recursive = TRUE) 

# Combine thinned records with environmental data
speciesEnv <- base::data.frame(raster::extract(predictors,
                                               cbind(species$lon,
                                                     species$lat)))

# Combine thinned records with environmental data 
speciesSwd <- cbind(species, speciesEnv)

# Visualise distribution of predictor variables at GPS points
# - Check that there are no outliers in the covariates 
hist(speciesSwd$bio1)
hist(speciesSwd$bio2)
hist(speciesSwd$bio3)
hist(speciesSwd$bio6)
hist(speciesSwd$bio8)
hist(speciesSwd$bio12)
hist(speciesSwd$bio14)
hist(speciesSwd$bio19)

# These are the column names of env predictors which were not colinear
climPred <- c("bio1", "bio2", "bio3", "bio6", "bio8",
              "bio12", "bio14", "bio19")

# Create dummy training data (only WORLDCLIM extracted values) for presence/absence 
dummyTraindata <- rbind(speciesSwd[, climPred], 
                         bgPoints[, climPred]) 

# Add column to identify presence/absence
(n_data <- nrow(speciesSwd))
(n_bg <- nrow(bgPoints))
dummyTraindata$presence <- c(replicate(n_data, "present"), 
                              replicate(n_bg, "absent"))
head(dummyTraindata)

# Make the vector of 1s and 0s, same length as there are rows in training data
presentBg <- ifelse(dummyTraindata$presence=='present', 1, 0)

# Create actual training data for MaxEnt model
trainData <- rbind(speciesSwd[ , climPred], 
                   bgPoints[ , climPred])

# -----------------------------------------------------------------------------
# Tune models using the 'ENMeval' package
# -----------------------------------------------------------------------------

# Tune MaxEnt models by RM and feature classes
# occ - lon, lat (two columns)
# env - raster of background study extent
# bg.coords - lon, lat (two columns of bg points)

# Subset environmental vars to just layers used to make models
reducedPred <- raster::subset(predictors, climPred)

# Run evaluation
# - This can take 30 minutes on my terrible PC
set.seed(2012)
modelTuning <- ENMevaluate(
  occ = species[, 2:3],
  env = reducedPred,
  bg.coords = bgPoints[, 20:21],
  RMvalues = c(1.0, 1.5, 2.0, 3.0,
               4.0, 5.0, 6.0, 8.0,
               10.0),
  clamp = TRUE,
  fc = c("HLQPT", "HLQ", "H"),
  method = "randomkfold",
  kfolds = c(5),
  rasterPreds = TRUE,
  algorithm = "maxent.jar"
)

# Extract the model tuning results 
# - Pretty messy code but meh, it works. 
settings <- modelTuning@results$settings
feat.class <- modelTuning@results$features
features <- modelTuning@results$features
rm <- modelTuning@results$rm
AUC.calibrate <- modelTuning@results$train.AUC
AUC.test <- modelTuning@results$avg.test.AUC
AUC.diff <- modelTuning@results$avg.diff.AUC
OR.10 <- modelTuning@results$avg.test.or10pct
OR.min <- modelTuning@results$avg.test.orMTP
AICc <- modelTuning@results$AICc
delta.AIC <- modelTuning@results$delta.AICc
parameters <-  modelTuning@results$parameters
modelEvaluation <- as.data.frame(cbind(settings, 
                                        rm, 
                                        AUC.calibrate, 
                                        AUC.test, 
                                        AUC.diff, 
                                        OR.10, 
                                        OR.min, 
                                        AICc, 
                                        delta.AIC,
                                        parameters))
modelEvaluation$feat.class <- as.factor(feat.class)
modelEvaluation$rm <- as.numeric(rm)
modelEvaluation$AUC.calibrate <- as.numeric(AUC.calibrate)
modelEvaluation$AUC.test <- as.numeric(AUC.test)
modelEvaluation$AUC.diff <- as.numeric(AUC.diff)
modelEvaluation$OR.10 <- as.numeric(OR.10)
modelEvaluation$delta.AIC <- as.numeric(delta.AIC)

# Save model tuning results to file for later use
readr::write_csv(
  x = modelEvaluation,
  file = here::here("./models/model_tuning/model_tuning_results.csv")
)

# Plot AUCtest
plot_aucTest <- ggplot(data = modelEvaluation, aes(x = rm,
                                                   y = AUC.test,
                                                   group = feat.class)) +
  geom_point(aes(colour = feat.class)) +
  geom_line(aes(colour = feat.class)) +
  scale_colour_grey(start = 0.4, end = 0.8) +
  #scale_y_continuous(breaks = seq(0.80, 1.00, 0.04),
  #                   limits = c(0.80, 1.00)) +
  labs(x = "Regularisation multiplier",
       y = expression(paste(AUC[test])),
       subtitle = "(a)")
plot_aucTest

# Plot AUCdiff
plot_aucDiff <- ggplot(data = modelEvaluation, aes(x = rm,
                                                   y = AUC.diff,
                                                   group = feat.class)) +
  geom_point(aes(colour = feat.class)) +
  geom_line(aes(colour = feat.class)) +
  scale_colour_grey(start = 0.4, end = 0.8) +
  #scale_y_continuous(breaks = seq(0.80, 1.00, 0.04),
  #                   limits = c(0.80, 1.00)) +
  labs(
    x = "Regularisation multiplier",
    y = expression(paste(AUC[diff])),
    subtitle = "(b)",
    colour = "Feature class"
  ) +
  theme(legend.position = c(0.8, 0.7))
plot_aucDiff

# Plot 10th percentile omission rates (OR10)
plot_or10 <- ggplot(data = modelEvaluation, aes(x = rm,
                                                y = OR.10,
                                                group = feat.class)) +
  geom_point(aes(colour = feat.class)) +
  geom_line(aes(colour = feat.class)) +
  scale_colour_grey(start = 0.4, end = 0.8) +
  geom_hline(yintercept = 0.10, linetype = "dashed") +
  labs(x = "Regularisation multiplier",
       y = expression(paste(OR[10])),
       subtitle = "(c)")
plot_or10

# Plot DeltaAICc
plot_aic <- ggplot(data = modelEvaluation, aes(x = rm,
                                               y = delta.AIC,
                                               group = feat.class)) +
  geom_point(aes(colour = feat.class)) +
  geom_line(aes(colour = feat.class)) +
  scale_colour_grey(start = 0.4, end = 0.8) +
  labs(x = "Regularisation multiplier",
       y = "AICc",
       subtitle = "(d)") +
  scale_y_continuous(breaks = seq(0, 800, 200),
                     limits = c(0, 800))
plot_aic

# Put the model tuning plots together
plots_modelTuning <- cowplot::plot_grid(plot_aucTest,
                                        plot_aucDiff,
                                        plot_or10,
                                        plot_aic,
                                        nrow = 2)
plots_modelTuning

# Save optimal model settings to file 
# Extract settings that maximised AUCtest
eval_tab_AUCtest <- modelEvaluation %>%
  dplyr::arrange(desc(AUC.test)) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(model = "AUCtest") %>%
  dplyr::select(model, rm, feat.class, AUC.test, AUC.diff, OR.10, delta.AIC) 
eval_tab_AUCtest

# Extract settings that minimised AUCdiff
eval_tab_AUCdiff <- modelEvaluation %>%
  dplyr::arrange(AUC.diff) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(model = "AUCdiff") %>%
  dplyr::select(model, rm, feat.class, AUC.test, AUC.diff, OR.10, delta.AIC) 
eval_tab_AUCdiff

# Extract settings that best approximated OR10
eval_tab_OR10 <- modelEvaluation %>%
  dplyr::arrange(OR.10) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(model = "OR10") %>%
  dplyr::select(model, rm, feat.class, AUC.test, AUC.diff, OR.10, delta.AIC) 
eval_tab_OR10

# Extract settings that minimised AICc
eval_tab_AICc <- modelEvaluation %>%
  dplyr::arrange(delta.AIC) %>%
  dplyr::slice(1) %>%
  dplyr::mutate(model = "AICc") %>%
  dplyr::select(model, rm, feat.class, AUC.test, AUC.diff, OR.10, delta.AIC)
eval_tab_AICc

# Default settings (default model for n > 80 GPS points has 
# rm = 1 and HLQPT features)
eval_tab_default <-modelEvaluation %>%
  dplyr::filter(rm == 1 & feat.class == "HLQPT") %>%
  dplyr::mutate(model = "Default settings") %>%
  dplyr::select(model, rm, feat.class, AUC.test, AUC.diff, OR.10, delta.AIC)
eval_tab_default

# Put the settings configurations together in a table
settings_table <- dplyr::bind_rows(
  eval_tab_default,
  eval_tab_AUCtest,
  eval_tab_AUCdiff,
  eval_tab_OR10,
  eval_tab_AICc
)
settings_table

# -----------------------------------------------------------------------------
# Run full MaxEnt models:
# -----------------------------------------------------------------------------

# Model 1: Default MaxEnt settings
# - The default model for n > 80 has rm = 1 and HLQPT features.
model_Default <- dismo::maxent(
  x = trainData,
  p = presentBg,
  path = paste(getwd(),
               './models/models/optimal_settings_default',
               sep = ''),
  args = c(
    'betamultiplier=1.0',
    'linear=true',
    'quadratic=true',
    'product=true',
    'threshold=true',
    'hinge=true',
    'threads=2',
    'doclamp=true',
    #'fadebyclamp=true',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=false',
    'responsecurves=true',
    'writemess=true',
    'writeplotdata=true',
    'writebackgroundpredictions=true'
  )
)

# Model 2: Optimise AUC.test
# - The top performing model by AUC.test had rm = 6 and hinge features only.
# Fit model with optimally tuned settings based on AICc
model_AUCtest <- dismo::maxent(
  x = trainData,
  p = presentBg,
  path = paste(getwd(),
               './models/models/optimal_settings_AUCtest',
               sep = ''),
  args = c(
    'betamultiplier=6.0',
    'linear=false',
    'quadratic=false',
    'product=false',
    'threshold=false',
    'hinge=true',
    'threads=2',
    'doclamp=true',
    #'fadebyclamp=true',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=false',
    'responsecurves=true',
    'writemess=true',
    'writeplotdata=true',
    'writebackgroundpredictions=true'
  )
)

# Model 3: Optimise AUC.diff
# - The top performing model by AUC.diff had rm = 4 and HLQPT features.
model_AUCdiff <- dismo::maxent(
  x = trainData,
  p = presentBg,
  path = paste(getwd(),
               './models/models/optimal_settings_AUCdiff',
               sep = ''),
  args = c(
    'betamultiplier=4',
    'linear=true',
    'quadratic=true',
    'product=true',
    'threshold=true',
    'hinge=true',
    'threads=2',
    'doclamp=true',
    #'fadebyclamp=true',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=false',
    'responsecurves=true',
    'writemess=true',
    'writeplotdata=true',
    'writebackgroundpredictions=true'
  )
)

# Model 4: Optimise OR10
# - The top performing model by OR10 had rm = 6 and only HLQPT features.
model_OR10 <- dismo::maxent(
  x = trainData,
  p = presentBg,
  path = paste(getwd(),
               './models/models/optimal_settings_OR10',
               sep = ''),
  args = c(
    'betamultiplier=6.0',
    'linear=true',
    'quadratic=true',
    'product=true',
    'threshold=true',
    'hinge=true',
    'threads=2',
    'doclamp=true',
    #'fadebyclamp=true',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=false',
    'responsecurves=true',
    'writemess=true',
    'writeplotdata=true',
    'writebackgroundpredictions=true'
  )
)

# Model 5: Optimise AICc
# - The top performing model by AICc had rm = 5 and only H features.
# Fit model with optimally tuned settings based on AICc
model_AICc <- dismo::maxent(
  x = trainData,
  p = presentBg,
  path = paste(getwd(),
               './models/models/optimal_settings_AICc',
               sep = ''),
  args = c(
    'betamultiplier=5.0',
    'linear=false',
    'quadratic=false',
    'product=false',
    'threshold=false',
    'hinge=true',
    'threads=2',
    'doclamp=true',
    #'fadebyclamp=true',
    'responsecurves=true',
    'jackknife=true',
    'askoverwrite=false',
    'responsecurves=true',
    'writemess=true',
    'writeplotdata=true',
    'writebackgroundpredictions=true'
  )
)

# -----------------------------------------------------------------------------
# Plot map - Default settings
# -----------------------------------------------------------------------------

# Get map of South Africa to project our model over
world <- rnaturalearth::ne_countries(scale = "medium",
                                     returnclass = "sf") %>%
  dplyr::filter(name == "South Africa")

# Extract MaxEnt raster output
# - Project MaxEnt scores over South Africa
afrExtent <- raster::mask(reducedPred, world)
predictMaxent <- raster::predict(model_Default, afrExtent)

# MaxEnt scores are in a raster layer above, but we need the MaxEnt scores
# in a dataframe. 
# - Below, create dataframe of MaxEnt model projection/scores
rcp <- raster::rasterToPoints(predictMaxent)
rcpdf <- base::data.frame(rcp)
colnames(rcpdf) <- c("Longitude", "Latitude", "Species")

# Plot map
final_map_default <- ggplot(data=rcpdf) +
  # Plot MaxEnt scores
  geom_tile(data = rcpdf, aes(x = Longitude,
                              y = Latitude, 
                              fill = Species)) +
  # Colour each grid cell according the MaxEnt scores
  # Scores close to 1 = High climatic suitability (good)
  # Scores close to 0 = Low climatic suitability (bad)
  scale_fill_gradientn(colours=c("white", "blue", "lightgreen",
                                 "yellow", "orange", "red"),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                       limits = c(0, 1)) +
  # Plot the shapefile of Africa
  geom_sf(data = world, fill = NA) +
  #scale_x_continuous(limits = c(10, 40)) + 
  #scale_y_continuous(limits = c(-35, -28)) +
  # Crops Africa to just the geographic extent of South Africa
  coord_sf(xlim = c(15, 33.5), 
           ylim = c(-35, -21.5), 
           crs = 4326, 
           expand = FALSE) +
  # Create title for the legend
  labs(fill = "Climatic similarity") + 
  # Add scale bar to bottom-right of map
  annotation_scale(location = "br", # 'br' = bottom right
                   style = "ticks", 
                   width_hint = 0.2) +
  # Add north arrow
  annotation_north_arrow(location = "br", 
                         which_north = "true", 
                         pad_x = unit(0.135, "in"), 
                         pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  # Apply the theme for the map we defined above. 
  theme_opts +
  theme(legend.position = "right") +
  # Change appearance of the legend
  guides(fill = guide_colorbar(ticks = FALSE),
         colour = guide_legend(order = 2))
final_map_default

# Save the plot to your PC 
ggsave("./models/figures/map_maxent_default.png",
       width = 8,
       height = 6,
       dpi = 600)

# Save raster
raster::writeRaster(
  predictMaxent,
  filename = "./models/rasters/raster_default",
  options = "interleave=band",
  overwrite = T
)
raster::writeRaster(
  predictMaxent,
  filename = "./models/rasters/raster_default",
  format = "GTiff",
  overwrite = T
)
raster::writeRaster(
  predictMaxent,
  filename = "./models/rasters/raster_default",
  format = "ascii",
  overwrite = T
)

# -----------------------------------------------------------------------------
# Plot map - AUCtest
# -----------------------------------------------------------------------------

# Extract MaxEnt raster output
# - Project MaxEnt scores over South Africa
predictMaxent <- raster::predict(model_AUCtest, afrExtent)

# MaxEnt scores are in a raster layer above, but we need the MaxEnt scores
# in a dataframe. 
# - Below, create dataframe of MaxEnt model projection/scores
rcp <- raster::rasterToPoints(predictMaxent)
rcpdf <- base::data.frame(rcp)
colnames(rcpdf) <- c("Longitude", "Latitude", "Species")

# Plot map
final_map_AUCtest <- ggplot(data=rcpdf) +
  # Plot MaxEnt scores
  geom_tile(data = rcpdf, aes(x = Longitude,
                              y = Latitude, 
                              fill = Species)) +
  # Colour each grid cell according the MaxEnt scores
  # Scores close to 1 = High climatic suitability (good)
  # Scores close to 0 = Low climatic suitability (bad)
  scale_fill_gradientn(colours=c("white", "blue", "lightgreen",
                                 "yellow", "orange", "red"),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                       limits = c(0, 1)) +
  # Plot the shapefile of Africa
  geom_sf(data = world, fill = NA) +
  #scale_x_continuous(limits = c(10, 40)) + 
  #scale_y_continuous(limits = c(-35, -28)) +
  # Crops Africa to just the geographic extent of South Africa
  coord_sf(xlim = c(15, 33.5), 
           ylim = c(-35, -21.5), 
           crs = 4326, 
           expand = FALSE) +
  # Create title for the legend
  labs(fill = "Climatic similarity") + 
  # Add scale bar to bottom-right of map
  annotation_scale(location = "br", # 'br' = bottom right
                   style = "ticks", 
                   width_hint = 0.2) +
  # Add north arrow
  annotation_north_arrow(location = "br", 
                         which_north = "true", 
                         pad_x = unit(0.135, "in"), 
                         pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  # Apply the theme for the map we defined above. 
  theme_opts +
  theme(legend.position = "right") +
  # Change appearance of the legend
  guides(fill = guide_colorbar(ticks = FALSE),
         colour = guide_legend(order = 2))
final_map_AUCtest

# Save the plot to your PC 
ggsave("./models/figures/map_maxent_AUCtest.png",
       width = 8,
       height = 6,
       dpi = 600)

# Save raster
raster::writeRaster(
  predictMaxent,
  filename = "./models/rasters/raster_AUCtest",
  options = "interleave=band",
  overwrite = T
)
raster::writeRaster(
  predictMaxent,
  filename = "./models/rasters/raster_AUCtest",
  format = "GTiff",
  overwrite = T
)
raster::writeRaster(
  predictMaxent,
  filename = "./models/rasters/raster_AUCtest",
  format = "ascii",
  overwrite = T
)

# -----------------------------------------------------------------------------
# Plot map - AUCdiff
# -----------------------------------------------------------------------------

# Extract MaxEnt raster output
# - Project MaxEnt scores over South Africa
predictMaxent <- raster::predict(model_AUCdiff, afrExtent)

# MaxEnt scores are in a raster layer above, but we need the MaxEnt scores
# in a dataframe. 
# - Below, create dataframe of MaxEnt model projection/scores
rcp <- raster::rasterToPoints(predictMaxent)
rcpdf <- base::data.frame(rcp)
colnames(rcpdf) <- c("Longitude", "Latitude", "Species")

# Plot map
final_map_AUCdiff <- ggplot(data=rcpdf) +
  # Plot MaxEnt scores
  geom_tile(data = rcpdf, aes(x = Longitude,
                              y = Latitude, 
                              fill = Species)) +
  # Colour each grid cell according the MaxEnt scores
  # Scores close to 1 = High climatic suitability (good)
  # Scores close to 0 = Low climatic suitability (bad)
  scale_fill_gradientn(colours=c("white", "blue", "lightgreen",
                                 "yellow", "orange", "red"),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                       limits = c(0, 1)) +
  # Plot the shapefile of Africa
  geom_sf(data = world, fill = NA) +
  #scale_x_continuous(limits = c(10, 40)) + 
  #scale_y_continuous(limits = c(-35, -28)) +
  # Crops Africa to just the geographic extent of South Africa
  coord_sf(xlim = c(15, 33.5), 
           ylim = c(-35, -21.5), 
           crs = 4326, 
           expand = FALSE) +
  # Create title for the legend
  labs(fill = "Climatic similarity") + 
  # Add scale bar to bottom-right of map
  annotation_scale(location = "br", # 'br' = bottom right
                   style = "ticks", 
                   width_hint = 0.2) +
  # Add north arrow
  annotation_north_arrow(location = "br", 
                         which_north = "true", 
                         pad_x = unit(0.135, "in"), 
                         pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  # Apply the theme for the map we defined above. 
  theme_opts +
  theme(legend.position = "right") +
  # Change appearance of the legend
  guides(fill = guide_colorbar(ticks = FALSE),
         colour = guide_legend(order = 2))
final_map_AUCdiff

# Save the plot to your PC 
ggsave("./models/figures/map_maxent_AUCdiff.png",
       width = 8,
       height = 6,
       dpi = 600)

# Save raster
raster::writeRaster(
  predictMaxent,
  filename = "./models/rasters/raster_AUCdiff",
  options = "interleave=band",
  overwrite = T
)
raster::writeRaster(
  predictMaxent,
  filename = "./models/rasters/raster_AUCdiff",
  format = "GTiff",
  overwrite = T
)
raster::writeRaster(
  predictMaxent,
  filename = "./models/rasters/raster_AUCdiff",
  format = "ascii",
  overwrite = T
)

# -----------------------------------------------------------------------------
# Plot map - OR10
# -----------------------------------------------------------------------------

# Extract MaxEnt raster output
# - Project MaxEnt scores over South Africa
predictMaxent <- raster::predict(model_OR10, afrExtent)

# MaxEnt scores are in a raster layer above, but we need the MaxEnt scores
# in a dataframe. 
# - Below, create dataframe of MaxEnt model projection/scores
rcp <- raster::rasterToPoints(predictMaxent)
rcpdf <- base::data.frame(rcp)
colnames(rcpdf) <- c("Longitude", "Latitude", "Species")

# Plot map
final_map_OR10 <- ggplot(data=rcpdf) +
  # Plot MaxEnt scores
  geom_tile(data = rcpdf, aes(x = Longitude,
                              y = Latitude, 
                              fill = Species)) +
  # Colour each grid cell according the MaxEnt scores
  # Scores close to 1 = High climatic suitability (good)
  # Scores close to 0 = Low climatic suitability (bad)
  scale_fill_gradientn(colours=c("white", "blue", "lightgreen",
                                 "yellow", "orange", "red"),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                       limits = c(0, 1)) +
  # Plot the shapefile of Africa
  geom_sf(data = world, fill = NA) +
  #scale_x_continuous(limits = c(10, 40)) + 
  #scale_y_continuous(limits = c(-35, -28)) +
  # Crops Africa to just the geographic extent of South Africa
  coord_sf(xlim = c(15, 33.5), 
           ylim = c(-35, -21.5), 
           crs = 4326, 
           expand = FALSE) +
  # Create title for the legend
  labs(fill = "Climatic similarity") + 
  # Add scale bar to bottom-right of map
  annotation_scale(location = "br", # 'br' = bottom right
                   style = "ticks", 
                   width_hint = 0.2) +
  # Add north arrow
  annotation_north_arrow(location = "br", 
                         which_north = "true", 
                         pad_x = unit(0.135, "in"), 
                         pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  # Apply the theme for the map we defined above. 
  theme_opts +
  theme(legend.position = "right") +
  # Change appearance of the legend
  guides(fill = guide_colorbar(ticks = FALSE),
         colour = guide_legend(order = 2))
final_map_OR10

# Save the plot to your PC 
ggsave("./models/figures/map_maxent_OR10.png",
       width = 8,
       height = 6,
       dpi = 600)

# Save raster
raster::writeRaster(
  predictMaxent,
  filename = "./models/rasters/raster_OR10",
  options = "interleave=band",
  overwrite = T
)
raster::writeRaster(
  predictMaxent,
  filename = "./models/rasters/raster_OR10",
  format = "GTiff",
  overwrite = T
)
raster::writeRaster(
  predictMaxent,
  filename = "./models/rasters/raster_OR10",
  format = "ascii",
  overwrite = T
)

# -----------------------------------------------------------------------------
# Plot map - AICc
# -----------------------------------------------------------------------------

# Extract MaxEnt raster output
# - Project MaxEnt scores over South Africa
predictMaxent <- raster::predict(model_AICc, afrExtent)

# MaxEnt scores are in a raster layer above, but we need the MaxEnt scores
# in a dataframe. 
# - Below, create dataframe of MaxEnt model projection/scores
rcp <- raster::rasterToPoints(predictMaxent)
rcpdf <- base::data.frame(rcp)
colnames(rcpdf) <- c("Longitude", "Latitude", "Species")

# Plot map
final_map_AICc <- ggplot(data=rcpdf) +
  # Plot MaxEnt scores
  geom_tile(data = rcpdf, aes(x = Longitude,
                              y = Latitude, 
                              fill = Species)) +
  # Colour each grid cell according the MaxEnt scores
  # Scores close to 1 = High climatic suitability (good)
  # Scores close to 0 = Low climatic suitability (bad)
  scale_fill_gradientn(colours=c("white", "blue", "lightgreen",
                                 "yellow", "orange", "red"),
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                       limits = c(0, 1)) +
  # Plot the shapefile of Africa
  geom_sf(data = world, fill = NA) +
  #scale_x_continuous(limits = c(10, 40)) + 
  #scale_y_continuous(limits = c(-35, -28)) +
  # Crops Africa to just the geographic extent of South Africa
  coord_sf(xlim = c(15, 33.5), 
           ylim = c(-35, -21.5), 
           crs = 4326, 
           expand = FALSE) +
  # Create title for the legend
  labs(fill = "Climatic similarity") + 
  # Add scale bar to bottom-right of map
  annotation_scale(location = "br", # 'br' = bottom right
                   style = "ticks", 
                   width_hint = 0.2) +
  # Add north arrow
  annotation_north_arrow(location = "br", 
                         which_north = "true", 
                         pad_x = unit(0.135, "in"), 
                         pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  # Apply the theme for the map we defined above. 
  theme_opts +
  theme(legend.position = "right") +
  # Change appearance of the legend
  guides(fill = guide_colorbar(ticks = FALSE),
         colour = guide_legend(order = 2))
final_map_AICc

# Save the plot to your PC 
ggsave("./models/figures/map_maxent_AICc.png",
       width = 8,
       height = 6,
       dpi = 600)

# Save raster
raster::writeRaster(
  predictMaxent,
  filename = "./models/rasters/raster_AICc",
  options = "interleave=band",
  overwrite = T
)
raster::writeRaster(
  predictMaxent,
  filename = "./models/rasters/raster_AICc",
  format = "GTiff",
  overwrite = T
)
raster::writeRaster(
  predictMaxent,
  filename = "./models/rasters/raster_AICc",
  format = "ascii",
  overwrite = T
)