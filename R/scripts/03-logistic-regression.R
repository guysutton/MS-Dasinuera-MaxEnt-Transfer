# -----------------------------------------------------------------------------
# Script: 03 - Logistic regression
#            - Do MaxEnt scores correlate with higher insect establishment?
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
  rms,
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

# Define function to calculate uncertainty in out-of-sample metrics
boot_val <- function(data, formula, ...) {
  out <- list()
  for (i in 1:500) {
    df <- sample_n(data, nrow(data), replace = TRUE)
    md <- glm(formula, data = df, family = binomial)
    out[[i]] <- val.prob(predict(md, type = "response"),
                         as.numeric(df$established) - 1,
                         pl = FALSE) %>% round(3)
  }
  return(out)
}

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
print(mod_lrm_default)
summary(mod_lrm_default)

# Calculate beta for 0.01 unit increase in clim_match
beta = 7.53 # Input the effect from clim_match row here 
exp(beta * 0.01)

# Validate model on out-of-sample AUC
d <- clim_data %>%
  dplyr::select(clim_match, established)
y <- clim_data %>%
  dplyr::select(established)
pred.logit <- predict(mod_lrm_default, d)
phat <- 1/(1+exp(-pred.logit))
v <- val.prob(phat, 
              y[1:77, ])
print(v)

# Calculate uncertainty
set.seed(48572)
boot_vals_m1 <- boot_val(data = clim_data,  
                         formula = as.formula("established ~ clim_match"))
boot_vals_m1 <- as.data.frame(boot_vals_m1, row.names = NULL)
head(boot_vals_m1)

# Extract AUC data
AUC_data <- boot_vals_m1 %>%
  dplyr::slice(2) %>%
  tidyr::pivot_longer(cols = 1:500,
                      names_to = "iteration",
                      values_to = "AUC_scores") %>%
  dplyr::mutate(iteration = 1:500)
head(AUC_data)

# Plot the predictions from the model 
mod <- glm(established ~ clim_match, 
                   data = clim_data,
                   family = binomial(link = "logit"))

# Get the link function for the binomial model
# - We need this to back-transform the predictions 
ilink <- family(mod)$linkinv

# Create some data to predict at: 100 values over the range of clim_match
ndata <- expand.grid(# Range of female mass to predict over 
  clim_match = seq(0, 1, 
                   # How many points?
                   l = 100))

# Add the fitted values by predicting from the model for the new data
ndata <- add_column(ndata, 
                    fit = predict(mod, 
                                  newdata = ndata,
                                  type = 'response'))

# Add predictions from the model to a new df
ndata <- bind_cols(ndata, setNames(as_tibble(predict(mod, 
                                                     ndata,
                                                     se.fit = TRUE)[1:2]),
                                   c("fit_link", "se_link")))
head(ndata)


# Create the confidence interval and back-transform
ndata <- ndata %>%
  dplyr::mutate(
    fit_resp = ilink(fit_link),
    right_upr = ilink(fit_link + (2 * se_link)),
    right_lwr = ilink(fit_link - (2 * se_link))
  )
head(ndata)

# Make the plot
logit_default <- ggplot(data = ndata) +
  geom_ribbon(aes(x = clim_match,
                  ymin = right_lwr,
                  ymax = right_upr),
              fill = "gray40",
              alpha = 0.2) +
  geom_line(aes(x = clim_match, 
                y = fit_resp)) +
  geom_rug(data = subset(clim_data, established == 1),
           aes(x = clim_match),
           sides = "t") +
  geom_rug(data = subset(clim_data, established == 0),
           aes(x = clim_match),
           sides = "b") +
  coord_cartesian(clip = "off") + 
  labs(x = "Maxent suitability score \n(climatic similarity to native range)",
       y = expression(paste("Probability of ",italic("D. rubiformis "), establishment)),
       subtitle = "(a) Default MaxEnt settings")
logit_default

# Save
ggsave("./models/figures/probability_of_establishment_default.png",
       dpi = 600,
       height = 6, 
       width = 7)

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
print(mod_lrm_AUCtest)
summary(mod_lrm_AUCtest)

# Calculate beta for 0.01 unit increase in clim_match
beta = 9.02 # Input the effect from clim_match row here 
exp(beta * 0.01)

# Validate model on out-of-sample AUC
d <- clim_data %>%
  dplyr::select(clim_match, established)
y <- clim_data %>%
  dplyr::select(established)
pred.logit <- predict(mod_lrm_AUCtest, d)
phat <- 1/(1+exp(-pred.logit))
v <- val.prob(phat, 
              y[1:77, ])
print(v)

# Calculate uncertainty
set.seed(48572)
boot_vals_m1 <- boot_val(data = clim_data,  
                         formula = as.formula("established ~ clim_match"))
boot_vals_m1 <- as.data.frame(boot_vals_m1, row.names = NULL)
head(boot_vals_m1)

# Extract AUC data
AUC_data <- boot_vals_m1 %>%
  dplyr::slice(2) %>%
  tidyr::pivot_longer(cols = 1:500,
                      names_to = "iteration",
                      values_to = "AUC_scores") %>%
  dplyr::mutate(iteration = 1:500)
head(AUC_data)

# Plot the predictions from the model 
mod <- glm(established ~ clim_match, 
           data = clim_data,
           family = binomial(link = "logit"))

# Get the link function for the binomial model
# - We need this to back-transform the predictions 
ilink <- family(mod)$linkinv

# Create some data to predict at: 100 values over the range of clim_match
ndata <- expand.grid(# Range of female mass to predict over 
  clim_match = seq(0, 1, 
                   # How many points?
                   l = 100))

# Add the fitted values by predicting from the model for the new data
ndata <- add_column(ndata, 
                    fit = predict(mod, 
                                  newdata = ndata,
                                  type = 'response'))

# Add predictions from the model to a new df
ndata <- bind_cols(ndata, setNames(as_tibble(predict(mod, 
                                                     ndata,
                                                     se.fit = TRUE)[1:2]),
                                   c("fit_link", "se_link")))
head(ndata)


# Create the confidence interval and back-transform
ndata <- ndata %>%
  dplyr::mutate(
    fit_resp = ilink(fit_link),
    right_upr = ilink(fit_link + (2 * se_link)),
    right_lwr = ilink(fit_link - (2 * se_link))
  )
head(ndata)

# Make the plot
logit_AUCtest <- ggplot(data = ndata) +
  geom_ribbon(aes(x = clim_match,
                  ymin = right_lwr,
                  ymax = right_upr),
              fill = "gray40",
              alpha = 0.2) +
  geom_line(aes(x = clim_match, 
                y = fit_resp)) +
  geom_rug(data = subset(clim_data, established == 1),
           aes(x = clim_match),
           sides = "t") +
  geom_rug(data = subset(clim_data, established == 0),
           aes(x = clim_match),
           sides = "b") +
  coord_cartesian(clip = "off") + 
  labs(x = "Maxent suitability score \n(climatic similarity to native range)",
       y = expression(paste("Probability of ",italic("D. rubiformis "), establishment)),
       subtitle = "(b) AUC<sub>*test*</sub>") +
  theme(plot.subtitle = element_markdown())
logit_AUCtest

# Save
ggsave("./models/figures/probability_of_establishment_AUCtest.png",
       dpi = 600,
       height = 6, 
       width = 7)

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
print(mod_lrm_AUCdiff)
summary(mod_lrm_AUCdiff)

# Calculate beta for 0.01 unit increase in clim_match
beta = 2.09 # Input the effect from clim_match row here
exp(beta * 0.01)

# Validate model on out-of-sample AUC
d <- clim_data %>%
  dplyr::select(clim_match, established)
y <- clim_data %>%
  dplyr::select(established)
pred.logit <- predict(mod_lrm_AUCdiff, d)
phat <- 1 / (1 + exp(-pred.logit))
v <- val.prob(phat,
              y[1:77,])
print(v)

# Calculate uncertainty
set.seed(48572)
boot_vals_m1 <- boot_val(data = clim_data,
                         formula = as.formula("established ~ clim_match"))
boot_vals_m1 <- as.data.frame(boot_vals_m1, row.names = NULL)
head(boot_vals_m1)

# Extract AUC data
AUC_data <- boot_vals_m1 %>%
  dplyr::slice(2) %>%
  tidyr::pivot_longer(cols = 1:500,
                      names_to = "iteration",
                      values_to = "AUC_scores") %>%
  dplyr::mutate(iteration = 1:500)
head(AUC_data)

# Plot the predictions from the model
mod <- glm(established ~ clim_match,
           data = clim_data,
           family = binomial(link = "logit"))

# Get the link function for the binomial model
# - We need this to back-transform the predictions
ilink <- family(mod)$linkinv

# Create some data to predict at: 100 values over the range of clim_match
ndata <- expand.grid(# Range of female mass to predict over
  clim_match = seq(0, 1,
                   # How many points?
                   l = 100))

# Add the fitted values by predicting from the model for the new data
ndata <- add_column(ndata,
                    fit = predict(mod,
                                  newdata = ndata,
                                  type = 'response'))

# Add predictions from the model to a new df
ndata <- bind_cols(ndata, setNames(as_tibble(predict(mod,
                                                     ndata,
                                                     se.fit = TRUE)[1:2]),
                                   c("fit_link", "se_link")))
head(ndata)


# Create the confidence interval and back-transform
ndata <- ndata %>%
  dplyr::mutate(
    fit_resp = ilink(fit_link),
    right_upr = ilink(fit_link + (2 * se_link)),
    right_lwr = ilink(fit_link - (2 * se_link))
  )
head(ndata)

# Make the plot
logit_AUCdiff <- ggplot(data = ndata) +
  geom_ribbon(
    aes(x = clim_match,
        ymin = right_lwr,
        ymax = right_upr),
    fill = "gray40",
    alpha = 0.2
  ) +
  geom_line(aes(x = clim_match,
                y = fit_resp)) +
  geom_rug(data = subset(clim_data, established == 1),
           aes(x = clim_match),
           sides = "t") +
  geom_rug(data = subset(clim_data, established == 0),
           aes(x = clim_match),
           sides = "b") +
  coord_cartesian(clip = "off") +
  labs(x = "Maxent suitability score \n(climatic similarity to native range)",
       y = expression(paste(
         "Probability of ", italic("D. rubiformis "), establishment
       )),
       subtitle = "(c) AUC<sub>*diff*</sub>") +
  theme(plot.subtitle = element_markdown())
logit_AUCdiff

# Save
ggsave(
  "./models/figures/probability_of_establishment_AUCdiff.png",
  dpi = 600,
  height = 6,
  width = 7
)

# ----------------------------------------------------------------------------
# Evaluate performance of AUCdiff model:
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
print(mod_lrm_OR10)
summary(mod_lrm_OR10)

# Calculate beta for 0.01 unit increase in clim_match
beta = 2.01 # Input the effect from clim_match row here
exp(beta * 0.01)

# Validate model on out-of-sample AUC
d <- clim_data %>%
  dplyr::select(clim_match, established)
y <- clim_data %>%
  dplyr::select(established)
pred.logit <- predict(mod_lrm_AUCdiff, d)
phat <- 1 / (1 + exp(-pred.logit))
v <- val.prob(phat,
              y[1:77,])
print(v)

# Calculate uncertainty
set.seed(48572)
boot_vals_m1 <- boot_val(data = clim_data,
                         formula = as.formula("established ~ clim_match"))
boot_vals_m1 <- as.data.frame(boot_vals_m1, row.names = NULL)
head(boot_vals_m1)

# Extract AUC data
AUC_data <- boot_vals_m1 %>%
  dplyr::slice(2) %>%
  tidyr::pivot_longer(cols = 1:500,
                      names_to = "iteration",
                      values_to = "AUC_scores") %>%
  dplyr::mutate(iteration = 1:500)
head(AUC_data)

# Plot the predictions from the model
mod <- glm(established ~ clim_match,
           data = clim_data,
           family = binomial(link = "logit"))

# Get the link function for the binomial model
# - We need this to back-transform the predictions
ilink <- family(mod)$linkinv

# Create some data to predict at: 100 values over the range of clim_match
ndata <- expand.grid(# Range of female mass to predict over
  clim_match = seq(0, 1,
                   # How many points?
                   l = 100))

# Add the fitted values by predicting from the model for the new data
ndata <- add_column(ndata,
                    fit = predict(mod,
                                  newdata = ndata,
                                  type = 'response'))

# Add predictions from the model to a new df
ndata <- bind_cols(ndata, setNames(as_tibble(predict(mod,
                                                     ndata,
                                                     se.fit = TRUE)[1:2]),
                                   c("fit_link", "se_link")))
head(ndata)

# Create the confidence interval and back-transform
ndata <- ndata %>%
  dplyr::mutate(
    fit_resp = ilink(fit_link),
    right_upr = ilink(fit_link + (2 * se_link)),
    right_lwr = ilink(fit_link - (2 * se_link))
  )
head(ndata)

# Make the plot
logit_OR10 <- ggplot(data = ndata) +
  geom_ribbon(
    aes(x = clim_match,
        ymin = right_lwr,
        ymax = right_upr),
    fill = "gray40",
    alpha = 0.2
  ) +
  geom_line(aes(x = clim_match,
                y = fit_resp)) +
  geom_rug(data = subset(clim_data, established == 1),
           aes(x = clim_match),
           sides = "t") +
  geom_rug(data = subset(clim_data, established == 0),
           aes(x = clim_match),
           sides = "b") +
  coord_cartesian(clip = "off") +
  labs(x = "Maxent suitability score \n(climatic similarity to native range)",
       y = expression(paste(
         "Probability of ", italic("D. rubiformis "), establishment
       )),
       subtitle = "(d) OR<sub>10</sub>") +
  theme(plot.subtitle = element_markdown())
logit_OR10

# Save
ggsave(
  "./models/figures/probability_of_establishment_OR10.png",
  dpi = 600,
  height = 6,
  width = 7
)

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
print(mod_lrm_AICc)
summary(mod_lrm_AICc)

# Calculate beta for 0.01 unit increase in clim_match
beta = 2.88 # Input the effect from clim_match row here
exp(beta * 0.01)

# Calculate uncertainty
set.seed(48572)
boot_vals_m1 <- boot_val(data = clim_data,
                         formula = as.formula("established ~ clim_match"))
boot_vals_m1 <- as.data.frame(boot_vals_m1, row.names = NULL)
head(boot_vals_m1)

# Extract AUC data
AUC_data <- boot_vals_m1 %>%
  dplyr::slice(2) %>%
  tidyr::pivot_longer(cols = 1:500,
                      names_to = "iteration",
                      values_to = "AUC_scores") %>%
  dplyr::mutate(iteration = 1:500)
head(AUC_data)

# Plot the predictions from the model
mod <- glm(established ~ clim_match,
           data = clim_data,
           family = binomial(link = "logit"))

# Get the link function for the binomial model
# - We need this to back-transform the predictions
ilink <- family(mod)$linkinv

# Create some data to predict at: 100 values over the range of clim_match
ndata <- expand.grid(# Range of female mass to predict over
  clim_match = seq(0, 1,
                   # How many points?
                   l = 100))

# Add the fitted values by predicting from the model for the new data
ndata <- add_column(ndata,
                    fit = predict(mod,
                                  newdata = ndata,
                                  type = 'response'))

# Add predictions from the model to a new df
ndata <- bind_cols(ndata, setNames(as_tibble(predict(mod,
                                                     ndata,
                                                     se.fit = TRUE)[1:2]),
                                   c("fit_link", "se_link")))
head(ndata)

# Create the confidence interval and back-transform
ndata <- ndata %>%
  dplyr::mutate(
    fit_resp = ilink(fit_link),
    right_upr = ilink(fit_link + (2 * se_link)),
    right_lwr = ilink(fit_link - (2 * se_link))
  )
head(ndata)

# Make the plot
logit_AICc <- ggplot(data = ndata) +
  geom_ribbon(
    aes(x = clim_match,
        ymin = right_lwr,
        ymax = right_upr),
    fill = "gray40",
    alpha = 0.2
  ) +
  geom_line(aes(x = clim_match,
                y = fit_resp)) +
  geom_rug(data = subset(clim_data, established == 1),
           aes(x = clim_match),
           sides = "t") +
  geom_rug(data = subset(clim_data, established == 0),
           aes(x = clim_match),
           sides = "b") +
  coord_cartesian(clip = "off") +
  labs(x = "Maxent suitability score \n(climatic similarity to native range)",
       y = expression(paste(
         "Probability of ", italic("D. rubiformis "), establishment
       )),
       subtitle = "(e) *AIC<sub>c</sub>*") +
  theme(plot.subtitle = element_markdown())
logit_AICc

# Save
ggsave(
  "./models/figures/probability_of_establishment_AICc.png",
  dpi = 600,
  height = 6,
  width = 7
)


# Make figure 
combined <- cowplot::plot_grid(logit_default,
                               logit_AUCtest,
                               logit_AUCdiff,
                               logit_OR10,
                               logit_AICc)
combined

# Save to PC
ggsave("./ms_body/ms_figs/fig_x_logit.png",
       dpi = 600,
       height = 8, 
       width = 12)



