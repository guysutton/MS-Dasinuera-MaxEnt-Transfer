# -----------------------------------------------------------------------------
# Script: - Reproduce Figure 7
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
# Plot map: 
# ----------------------------------------------------------------------------

# Import Maxent prediction raster
clim_map <- raster::raster("./models/rasters/raster_AICc.asc")

# Threshold map to where MaxEnt scores were 0.30 (>75% probability
# of establishment, see Fig. 3)
clim_map <- clim_map > 0.30

# Mask MaxEnt projection to South Africa
predictMaxent <- raster::mask(clim_map, world)
plot(predictMaxent)

# MaxEnt scores are in a raster layer above, but we need the MaxEnt scores
# in a dataframe.
rcp <- raster::rasterToPoints(predictMaxent)
rcpdf <- base::data.frame(rcp)
colnames(rcpdf) <- c("Longitude", "Latitude", "Species")
rcpdf$Species <- as.factor(rcpdf$Species)

# Plot map
ggplot(data = rcpdf) +
  # Plot MaxEnt scores
  geom_tile(data = rcpdf, aes(x = Longitude,
                              y = Latitude,
                              fill = Species)) +
  # Colour each grid cell according the MaxEnt scores
  scale_fill_manual(
    values = c("gray90", "gray30"),
    labels = c("Unsuitable", "Suitable")
  ) +
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
  labs(fill = "Climatic suitability") +
  annotation_scale(location = "br", # 'br' = bottom right
                   style = "ticks", 
                   width_hint = 0.2) +
  annotation_north_arrow(location = "br", 
                         which_north = "true", 
                         pad_x = unit(0.55, "in"), 
                         pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  # Apply the theme for the map we defined above.
  theme_opts +
  theme(legend.position = c(0.15, 0.9),
        plot.subtitle = element_markdown())

# Save the plot to your PC
ggsave(
  "./ms_body/ms_figs/fig_7_biocontrol.png",
  width = 10,
  height = 12,
  dpi = 600
)