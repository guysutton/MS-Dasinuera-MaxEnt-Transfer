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
spatialCorAnalysis <- ecospat::ecospat.mantel.correlogram(
  # Data frame with environmental variables
  dfvar = spatialData,
  # Columns containing lat/long
  colxy = 9:10,
  # Number of random occurrences
  n = 100,
  # Climate variables
  colvar = 1:8,
  max = 1000,
  nclass = 100,
  nperm = 100
)

# Add a species column to the GPS data
species <- species %>%
  dplyr::mutate(species = "D. rubiformis")

# Thin by spatial autocorrelation value
speciesThinned <- spThin::thin(
  loc.data = species,
  lat.col = "lat",
  long.col = "lon",
  spec.col = "species",
  # Km unit of correlogram
  thin.par = 0.2,
  reps = 100,
  max.files = 1,
  out.dir = here::here("./data/data_clean/")
)

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
# - This can take > 1 hour on my terrible PC
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
  labs(x = "Regularisation multiplier",
       y = expression(paste(AUC[diff])),
       subtitle = "(b)",
       colour = "Feature class") +
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
plots_modelTuning <- cowplot::plot_grid(
  plot_aucTest,
  plot_aucDiff,
  plot_or10,
  plot_aic,
  nrow = 2
)
plots_modelTuning


