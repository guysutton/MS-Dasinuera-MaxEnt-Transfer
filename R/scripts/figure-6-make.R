# -----------------------------------------------------------------------------
# Script: - Reproduce Figure 6
#
# AUTHOR: Guy F. Sutton
# AFFILIATION: Centre for Biological Control, Rhodes University, South Africa
# DATE MODIFIED: 03/02/2021
# CONTACT: g.sutton@ru.ac.za
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Session setup
# -----------------------------------------------------------------------------

# Load packages
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
  patchwork,
  rnaturalearth,
  rnaturalearthdata,
  rlang,
  sf,
  ggspatial,
  raster,
  here,
  tidytext,
  ggtext,
  dismo,
  tidymodels,
  themis,
  ranger,
  viridis,
  yardstick,
  vip,
  DHARMa,
  rms
)

# Set the theme for the figures
theme_set(
  theme_classic() +
    theme(
      panel.border = element_rect(colour = "black", fill = NA),
      axis.text = element_text(colour = "black"),
      axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
      axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
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
(mess <- plot(mss))

# -----------------------------------------------------------------------------
# Plot MESS map:
# -----------------------------------------------------------------------------

# Now plot MESS map nicely
# Raster should be the MESS map raster
test <- mss
test <- raster::mask(test, world)
rcp <- raster::rasterToPoints(test)
rcpdf <- base::data.frame(rcp)
colnames(rcpdf) <- c("Longitude", "Latitude", "Species")

mess_map <- ggplot(data=rcpdf) +
  geom_tile(data = rcpdf, aes(x = Longitude,
                              y = Latitude, 
                              fill=Species)) +
  geom_sf(data = world, fill = NA) +
  # Colour each grid cell according the MaxEnt scores
  scale_fill_viridis(limits = c(0, 1),
                     breaks = seq(0, 1, 0.25)) +
  #scale_x_continuous(limits = c(10, 40)) + 
  #scale_y_continuous(limits = c(-35, -28)) +
  coord_sf(xlim = c(-20, 52), 
           ylim = c(-37, 40), 
           crs = 4326, 
           expand = FALSE) +
  #geom_point(data = occ_gps, 
  #           aes(x = lon, y = lat),
  #           size = 1.5) + 
  #geom_text_repel(data = occ_gps, 
  #                aes(x = lon, 
  #                    y = lat, 
  #                    label = gen), 
  #                force = 1.0, 
  #                nudge_y = 1.0,
  #                nudge_x = 1.0) +
  labs(fill = "MESS") + 
  annotation_scale(location = "bl", # 'br' = bottom right
                   style = "ticks", 
                   width_hint = 0.2) +
  annotation_north_arrow(location = "bl", 
                         which_north = "true", 
                         pad_x = unit(0.135, "in"), 
                         pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_opts +
  guides(fill = guide_colorbar(ticks = FALSE),
         colour = guide_legend(order = 2)) +
  theme(legend.text = element_text(size = 15))
mess_map

# Change the MESS(-) regions into '0' and MESS(+) into '1'
r.mess.mask <- mss>0
plot(r.mess.mask)

# Raster should be the MESS map raster
test <- r.mess.mask
test <- raster::crop(test, world)
test <- raster::mask(test, world)

# Convert from raster to data frame
rcp <- rasterToPoints(test)
rcpdf <- data.frame(rcp)
colnames(rcpdf) <- c("Longitude", "Latitude", "Species")
rcpdf$Species <- as.factor(rcpdf$Species)

# Plot binary map
fig6a <- ggplot(data=rcpdf) +
  geom_tile(data = rcpdf, aes(x = Longitude,
                              y = Latitude, 
                              fill=Species)) +
  scale_fill_manual(values=c("gray90", "gray30"),
                    labels = c("Extrapolate", "Interpolate")) +
  geom_sf(data = world, fill = NA) +
  coord_sf(
    xlim = c(16, 33.5), 
    ylim = c(-21.5, -35), 
    crs = 4326, 
    expand = FALSE) +
  labs(fill = "MESS",
       subtitle = "(a)") + 
  annotation_scale(location = "br", # 'br' = bottom right
                   style = "ticks", 
                   width_hint = 0.2) +
  annotation_north_arrow(location = "br", 
                         which_north = "true", 
                         pad_x = unit(0.3, "in"), 
                         pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_opts +
  theme(legend.text = element_text(size = 10),
        legend.position = c(0.15, 0.9))
fig6a

# -----------------------------------------------------------------------------
# Make results plot:
# -----------------------------------------------------------------------------

# Mess map
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

# Run basic model 
mod1 <- glm(established ~ mess_match, 
            data = clim_data,
            family = binomial)
summary(mod1)

# Extract the maxent scores for each site 
clim_data <- raster::extract(r.mess.mask,
                             inv_clean[,1:2])
head(clim_data)
clim_data <- as.data.frame(clim_data)
head(clim_data)

# Bind the maxent scores to the site data 
clim_data <- cbind(inv_clean, clim_data)
clim_data

# Drop the NA values 
clim_data <- clim_data %>%
  dplyr::select(lon,
                lat,
                established,
                clim_match = clim_data) %>%
  tidyr::drop_na(clim_match)
head(clim_data)

# Calculate cross-tabulation of establishment or not versus MESS score
mess_tab <- clim_data %>%
  dplyr::group_by(clim_match, established) %>%
  dplyr::count()
mess_tab

# Add n and calculate proportion
total <- c(54, 54, 24, 24)
mess_tab$total <- total
mess_tab <- mess_tab %>%
  ungroup() %>%
  mutate(prop = n/total * 100,
         established = as.factor(established),
         mess_match = as.factor(clim_match),
         established = recode(established, 
                              "1" = "Yes",
                              "0" = "No"),
         mess_match = recode(mess_match, 
                             "1" = "Interpolation",
                             "0" = "Extrapolation"))
mess_tab 


# Plot
fig6b <- ggplot(data = mess_tab, aes(x = mess_match,
                                     y = prop,
                                     group = established)) +
  geom_bar(aes(fill = established),
           position = "dodge",
           stat = "identity") +
  scale_fill_manual(values = c("grey80", "grey60")) +
  #geom_text(aes(label = n),
  #          position = position_dodge(0.9),
  #          vjust = -0.5) +
  labs(x = "MESS scores \n(climate analogous to native range)",
       y = expression(paste(
         "Proportion of sites with ",
         italic("D. rubiformis "),
         "establishment"
       )),
       fill = "Established",
       subtitle = "(b)") +
  theme(legend.position = "right")
fig6b

# ----------------------------------------------------------------------------
# Make Figure 6
# ----------------------------------------------------------------------------

# Put all the maps together in one figure
all_maps <- fig6a + fig6b + 
  plot_layout(nrow = 1, byrow = TRUE)
all_maps

# Save the plot to your PC
ggsave(
  "./ms_body/ms_figs/fig_6_mess.png",
  width = 14,
  height = 16,
  dpi = 600
)
