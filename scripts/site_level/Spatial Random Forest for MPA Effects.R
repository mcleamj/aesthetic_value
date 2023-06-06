

###########################################################
#'
#' SPATIAL RANDOM FOREST MODEL 
#' TO ASSESS THE SPATIAL SIGNAL OF MPA EFFECTS 
#' 
#' 
#########################################################

# LIBRARY PACKAGES
library(geodist)

# LOAD DATA
model_data <- readr::read_rds("data/model_data.rds")

# SELECT VARIABLES

rf_data <- model_data %>%
  select(aesthe_survey, MPA, HDI2017, fshD, gravtot2, SiteLongitude, SiteLatitude, SiteCode) %>%
  group_by(SiteCode) %>%
  summarise_all()

# PRE AVERAGE SURVEYS TO SITE LEVEL

rf_data <- 

rf_data <- log_aesth <- log(rf_data$aesthe_survey)

dependent.variable.name <- "log_aesth"
predictor.variable.names <- c("MPA", "HDI2017", "fshD", "gravtot2")

# CHECK FOR DATA REQUIREMENTS

sum(apply(rf_data, 2, is.na))

apply(rf_data, 2, var) == 0

sum(apply(scale(aesthe_survey), 2, is.nan))

# DEFINE SPATIAL COMPONENTS

xy <- rf_data %>%
  select(SiteLongitude, SiteLatitude) %>%
  rename(x = SiteLongitude, y = SiteLatitude)

dist.matrix <- geodist(xy)

distance.thresholds <- c()

## MAP 

world <- rnaturalearth::ne_countries(
  scale = "medium", 
  returnclass = "sf"
)

ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = world, 
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = RF_DATA,
    ggplot2::aes(
      x = x,
      y = y,
      color = log_aesth
    ),
    size = 2.5
  ) +
  ggplot2::scale_color_viridis_c(
    direction = -1, 
    option = "F"
  ) +
  ggplot2::theme_bw() +
  ggplot2::labs(color = "Log aesthetic") +
  ggplot2::scale_x_continuous(limits = c(-170, -30)) +
  ggplot2::scale_y_continuous(limits = c(-58, 80))  +
  ggplot2::ggtitle("Global Distribution of Aesthetic Value") + 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

# LOOK AT SPATIAL AUTOCORRELATION

spatialRF::plot_training_df_moran(
  data = rf_data,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix,
  distance.thresholds = distance.thresholds,
  fill.color = viridis::viridis(
    100,
    option = "F",
    direction = -1
  ),
  point.color = "gray40"
)

# NO NEED TO WORRY ABOUT MULTICOLLINEARTY IN THIS MODEL

# RUN THE SPATIAL RANDOM FOREST 

spatial.model <- spatialRF::rf_spatial(
  data = rf_data,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = distance.matrix
)


spatialRF::plot_residuals_diagnostics(
  spatial.model,
  verbose = FALSE
)


spatialRF::plot_importance(
  spatial.model,
  verbose = FALSE
)


# SITE LEVEL MODEL IMPORTANCE 

local.importance <- spatialRF::get_importance_local(model.non.spatial)

