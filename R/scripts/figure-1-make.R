# -----------------------------------------------------------------------------
# Script: - Reproduce Figure 1
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
  cowplot,
  tidymodels,
  themis,
  ranger,
  yardstick,
  vip,
  DHARMa
)

# Set theme for ggplot
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

# ----------------------------------------------------------------------------
# Make invaded-range map:
# ----------------------------------------------------------------------------

# Load invaded-range GPS points
data <-
  readr::read_csv("./data/data_clean/dasi_rubi_invaded_range_gps_clean.csv")

# Clean records
data <- data %>%
  janitor::clean_names() %>%
  dplyr::select(lat_converted, lon_converted, established) %>%
  # Make lat negative to be in southern hemisphere
  dplyr::mutate(lat_converted = lat_converted * -1) %>%
  dplyr::filter(lat_converted < -23)
head(data)


# Get world map
world <- rnaturalearth::ne_countries(scale = "medium",
                                     returnclass = "sf")

# Keep only South Africa
world <- world %>%
  dplyr::filter(name %in% c("South Africa"))

# Classify sites as established (1) or absent (0)
data <- data %>%
  dplyr::mutate(established = dplyr::case_when(established == 1 ~ "Established",
                                               TRUE ~ "Absent")) %>%
  dplyr::mutate(established = forcats::fct_relevel(established, "Absent", "Established"))

# Reorder Y variable to make established points plot over absent points 
data <- data[order(as.numeric(factor(data$established))),]

# Make map
set.seed(2021)
options(repr.plot.width = 1, repr.plot.height = 0.75)
p <- ggplot() +
  # Provincial boundaries
  geom_sf(data = world,
          aes(fill = NA),
          # Colour of the outline
          colour = "black",
          # Width of the province border lines
          size = 0.35) +
  scale_fill_brewer(palette = "YlOrBr") +
  geom_point(data = data,
             aes(x = lon_converted,
                 y = lat_converted,
                 colour = established)) +
  scale_colour_grey(start = 0.7, end = 0.3) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "(b) South Africa (introduced range)",
    colour = NULL
  ) +
  coord_sf(
    xlim = c(15.5, 33.5),
    ylim = c(-35,-21.75),
    expand = FALSE
  ) +
  annotation_scale(location = "br",
                   #
                   style = "ticks",
                   width_hint = 0.150) +
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(0.11, "in"),
    pad_y = unit(0.3, "in"),
    style = ggspatial::north_arrow_fancy_orienteering
  ) +
  theme(legend.position = c(0.15, 0.9))
p

# Save the plot to your PC
ggsave(
  "./models/figures/map_invaded_distribution.png",
  width = 10,
  height = 8,
  dpi = 600
)

# ----------------------------------------------------------------------------
# Make native-range map:
# ----------------------------------------------------------------------------

# Load native-range GPS points 
data <- readr::read_csv("./data/data_clean/dasi_rubi_native_range_gps_clean.csv")

# Clean records
data <- data %>%
  janitor::clean_names() %>%
  dplyr::select(lat_converted, lon_converted) %>%
  # Make lat negative to be in southern hemisphere 
  dplyr::mutate(lat_converted = lat_converted * -1) 
head(data)

# Get world map
world <- rnaturalearth::ne_countries(scale = "medium",
                                     returnclass = "sf")

# Keep only Australia
world <- world %>%
  dplyr::filter(name %in% c("Australia"))

# Make map
set.seed(2021)
options(repr.plot.width = 1, repr.plot.height = 0.75)
p1 <- ggplot() +
  # Provincial boundaries
  geom_sf(data = world, aes(fill = NA),
          # Colour of the outline
          colour = "black",
          # Width of the province border lines
          size = 0.35) +
  scale_fill_brewer(palette = "YlOrBr") +
  geom_point(data = data, aes(x = lon_converted,
                              y = lat_converted)) +
  #scale_colour_grey(start = 0.7, end = 0.3)+ 
  labs(x = "Longitude",
       y = "Latitude",
       title = "(a) Australia (native range)",
       colour = NULL) +
  #ylim(-45, -10) +
  #xlim(110, 155) +
  coord_sf(xlim = c(110, 155),
           ylim = c(-45, -8),
           expand = FALSE) +
  annotation_scale(location = "bl", #
                   style = "ticks",
                   width_hint = 0.150) +
  annotation_north_arrow(location = "bl",
                         which_north = "true",
                         pad_x = unit(0.165, "in"),
                         pad_y = unit(0.3, "in"),
                         style = ggspatial::north_arrow_fancy_orienteering) 
p1


# Save the plot to your PC 
ggsave("./models/figures/map_native_distribution.png",
       width = 10,
       height = 8,
       dpi = 600)

# ----------------------------------------------------------------------------
# Make Figure 1:
# ----------------------------------------------------------------------------

# Combine plots
plots <- plot_grid(p1, p, nrow = 1, align = "h")
plots

# Save the plot to your PC 
ggsave("./ms_body/ms_figs/fig_1_distributions.png",
       width = 12,
       height = 8,
       dpi = 600)


