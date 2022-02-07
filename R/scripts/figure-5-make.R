# -----------------------------------------------------------------------------
# Script: - Reproduce Figure 5
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

# ----------------------------------------------------------------------------
# Plot default map: 
# ----------------------------------------------------------------------------

# Import Maxent prediction raster
clim_map <- raster::raster("./models/rasters/raster_default.asc")

# Get map of South Africa to project our model over
world <- rnaturalearth::ne_countries(scale = "medium",
                                     returnclass = "sf") %>%
  dplyr::filter(name == "South Africa")

# Mask MaxEnt projection to South Africa
predictMaxent <- raster::mask(clim_map, world)

# MaxEnt scores are in a raster layer above, but we need the MaxEnt scores
# in a dataframe.
rcp <- raster::rasterToPoints(predictMaxent)
rcpdf <- base::data.frame(rcp)
colnames(rcpdf) <- c("Longitude", "Latitude", "Species")

# Plot map
final_map_default <- ggplot(data = rcpdf) +
  # Plot MaxEnt scores
  geom_tile(data = rcpdf, aes(x = Longitude,
                              y = Latitude,
                              fill = Species)) +
  # Colour each grid cell according the MaxEnt scores
  scale_fill_viridis(limits = c(0, 1),
                     breaks = seq(0, 1, 0.25)) +
  # Plot the shapefile of Africa
  geom_sf(data = world, fill = NA) +
  # Crops Africa to just the geographic extent of South Africa
  coord_sf(
    xlim = c(15, 33.5),
    ylim = c(-35,-21.5),
    crs = 4326,
    expand = FALSE
  ) +
  # Create title for the legend
  labs(fill = "Climatic suitability",
       subtitle = "(a) Default MaxEnt settings") +
  # Add scale bar to bottom-right of map
  annotation_scale(location = "br",
                   # 'br' = bottom right
                   style = "ticks",
                   width_hint = 0.2) +
  # Add north arrow
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(0.1, "in"),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  # Apply the theme for the map we defined above.
  theme_opts +
  theme(legend.position = "right") +
  # Change appearance of the legend
  guides(fill = guide_colorbar(ticks = FALSE),
         colour = guide_legend(order = 2))
final_map_default

# ----------------------------------------------------------------------------
# Plot AUCtest map: 
# ----------------------------------------------------------------------------

# Import Maxent prediction raster
clim_map <- raster::raster("./models/rasters/raster_AUCtest.asc")

# Mask MaxEnt projection to South Africa
predictMaxent <- raster::mask(clim_map, world)

# MaxEnt scores are in a raster layer above, but we need the MaxEnt scores
# in a dataframe.
rcp <- raster::rasterToPoints(predictMaxent)
rcpdf <- base::data.frame(rcp)
colnames(rcpdf) <- c("Longitude", "Latitude", "Species")

# Plot map
final_map_AUCtest <- ggplot(data = rcpdf) +
  # Plot MaxEnt scores
  geom_tile(data = rcpdf, aes(x = Longitude,
                              y = Latitude,
                              fill = Species)) +
  # Colour each grid cell according the MaxEnt scores
  scale_fill_viridis(limits = c(0, 1),
                     breaks = seq(0, 1, 0.25)) +
  # Plot the shapefile of Africa
  geom_sf(data = world, fill = NA) +
  # Crops Africa to just the geographic extent of South Africa
  coord_sf(
    xlim = c(15, 33.5),
    ylim = c(-35,-21.5),
    crs = 4326,
    expand = FALSE
  ) +
  # Create title for the legend
  labs(fill = "Climatic suitability",
       subtitle = "(b) AUC<sub>*test*</sub>") +
  # Add scale bar to bottom-right of map
  annotation_scale(location = "br",
                   # 'br' = bottom right
                   style = "ticks",
                   width_hint = 0.2) +
  # Add north arrow
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(0.1, "in"),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  # Apply the theme for the map we defined above.
  theme_opts +
  theme(legend.position = "right",
        plot.subtitle = element_markdown()) +
  # Change appearance of the legend
  guides(fill = guide_colorbar(ticks = FALSE),
         colour = guide_legend(order = 2))
final_map_AUCtest

# ----------------------------------------------------------------------------
# Plot AUCdiff map: 
# ----------------------------------------------------------------------------

# Import Maxent prediction raster
clim_map <- raster::raster("./models/rasters/raster_AUCdiff.asc")

# Mask MaxEnt projection to South Africa
predictMaxent <- raster::mask(clim_map, world)

# MaxEnt scores are in a raster layer above, but we need the MaxEnt scores
# in a dataframe.
rcp <- raster::rasterToPoints(predictMaxent)
rcpdf <- base::data.frame(rcp)
colnames(rcpdf) <- c("Longitude", "Latitude", "Species")

# Plot map
final_map_AUCdiff <- ggplot(data = rcpdf) +
  # Plot MaxEnt scores
  geom_tile(data = rcpdf, aes(x = Longitude,
                              y = Latitude,
                              fill = Species)) +
  # Colour each grid cell according the MaxEnt scores
  scale_fill_viridis(limits = c(0, 1),
                     breaks = seq(0, 1, 0.25)) +
  # Plot the shapefile of Africa
  geom_sf(data = world, fill = NA) +
  # Crops Africa to just the geographic extent of South Africa
  coord_sf(
    xlim = c(15, 33.5),
    ylim = c(-35,-21.5),
    crs = 4326,
    expand = FALSE
  ) +
  # Create title for the legend
  labs(fill = "Climatic suitability",
       subtitle = "(c) AUC<sub>*diff*</sub>") +
  # Add scale bar to bottom-right of map
  annotation_scale(location = "br",
                   # 'br' = bottom right
                   style = "ticks",
                   width_hint = 0.2) +
  # Add north arrow
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(0.1, "in"),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  # Apply the theme for the map we defined above.
  theme_opts +
  theme(legend.position = "right",
        plot.subtitle = element_markdown()) +
  # Change appearance of the legend
  guides(fill = guide_colorbar(ticks = FALSE),
         colour = guide_legend(order = 2))
final_map_AUCdiff

# ----------------------------------------------------------------------------
# Plot OR10 map: 
# ----------------------------------------------------------------------------

# Import Maxent prediction raster
clim_map <- raster::raster("./models/rasters/raster_OR10.asc")

# Mask MaxEnt projection to South Africa
predictMaxent <- raster::mask(clim_map, world)

# MaxEnt scores are in a raster layer above, but we need the MaxEnt scores
# in a dataframe.
rcp <- raster::rasterToPoints(predictMaxent)
rcpdf <- base::data.frame(rcp)
colnames(rcpdf) <- c("Longitude", "Latitude", "Species")

# Plot map
final_map_OR10 <- ggplot(data = rcpdf) +
  # Plot MaxEnt scores
  geom_tile(data = rcpdf, aes(x = Longitude,
                              y = Latitude,
                              fill = Species)) +
  # Colour each grid cell according the MaxEnt scores
  scale_fill_viridis(limits = c(0, 1),
                     breaks = seq(0, 1, 0.25)) +
  # Plot the shapefile of Africa
  geom_sf(data = world, fill = NA) +
  # Crops Africa to just the geographic extent of South Africa
  coord_sf(
    xlim = c(15, 33.5),
    ylim = c(-35,-21.5),
    crs = 4326,
    expand = FALSE
  ) +
  # Create title for the legend
  labs(fill = "Climatic suitability",
       subtitle = "(d) OR<sub>10</sub>") +
  # Add scale bar to bottom-right of map
  annotation_scale(location = "br",
                   # 'br' = bottom right
                   style = "ticks",
                   width_hint = 0.2) +
  # Add north arrow
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(0.1, "in"),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  # Apply the theme for the map we defined above.
  theme_opts +
  theme(legend.position = "right",
        plot.subtitle = element_markdown()) +
  # Change appearance of the legend
  guides(fill = guide_colorbar(ticks = FALSE),
         colour = guide_legend(order = 2))
final_map_OR10

# ----------------------------------------------------------------------------
# Plot AICc map: 
# ----------------------------------------------------------------------------

# Import Maxent prediction raster
clim_map <- raster::raster("./models/rasters/raster_AICc.asc")

# Mask MaxEnt projection to South Africa
predictMaxent <- raster::mask(clim_map, world)

# MaxEnt scores are in a raster layer above, but we need the MaxEnt scores
# in a dataframe.
rcp <- raster::rasterToPoints(predictMaxent)
rcpdf <- base::data.frame(rcp)
colnames(rcpdf) <- c("Longitude", "Latitude", "Species")

# Plot map
final_map_AICc <- ggplot(data = rcpdf) +
  # Plot MaxEnt scores
  geom_tile(data = rcpdf, aes(x = Longitude,
                              y = Latitude,
                              fill = Species)) +
  # Colour each grid cell according the MaxEnt scores
  scale_fill_viridis(limits = c(0, 1),
                     breaks = seq(0, 1, 0.25)) +
  # Plot the shapefile of Africa
  geom_sf(data = world, fill = NA) +
  # Crops Africa to just the geographic extent of South Africa
  coord_sf(
    xlim = c(15, 33.5),
    ylim = c(-35,-21.5),
    crs = 4326,
    expand = FALSE
  ) +
  # Create title for the legend
  labs(fill = "Climatic suitability",
       subtitle = "(e) AIC<sub>c</sub>") +
  # Add scale bar to bottom-right of map
  annotation_scale(location = "br",
                   # 'br' = bottom right
                   style = "ticks",
                   width_hint = 0.2) +
  # Add north arrow
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(0.1, "in"),
    pad_y = unit(0.3, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  # Apply the theme for the map we defined above.
  theme_opts +
  theme(legend.position = "right",
        plot.subtitle = element_markdown()) +
  # Change appearance of the legend
  guides(fill = guide_colorbar(ticks = FALSE),
         colour = guide_legend(order = 2))
final_map_AICc

# ----------------------------------------------------------------------------
# Save all the maps as one figure
# ----------------------------------------------------------------------------

# Put all the maps together in one figure
all_maps <- final_map_default + 
  final_map_AUCtest + 
  final_map_AUCdiff +
  final_map_OR10 +
  final_map_AICc +
  plot_layout(nrow = 2, byrow = TRUE)
all_maps

# Save the plot to your PC
ggsave(
  "./ms_body/ms_figs/fig_5_maxent_maps.png",
  width = 18,
  height = 12,
  dpi = 600
)
