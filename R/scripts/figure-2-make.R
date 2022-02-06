# -----------------------------------------------------------------------------
# Script: - Reproduce Figure 2
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
# Make Figure 2:

# The raw data required to reproduce Figure 2 is produced in script 
# '02-maxent-native-range.R'. Below, we import this raw data and then make 
# the final figure.
# ----------------------------------------------------------------------------

# Load raw data 
modelEvaluation <- readr::read_csv("./models/model_tuning/model_tuning_results.csv")
head(modelEvaluation)

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
       subtitle = "(d)") 
plot_aic

# Put the model tuning plots together
plots_modelTuning <- cowplot::plot_grid(plot_aucTest,
                                        plot_aucDiff,
                                        plot_or10,
                                        plot_aic,
                                        nrow = 2)
plots_modelTuning


# Save the plot to your PC 
ggsave("./ms_body/ms_figs/fig_2_model_tuning.png",
       width = 8,
       height = 6,
       dpi = 600)









