

###########################################################
#'
#' SPATIAL RANDOM FOREST MODEL 
#' TO ASSESS THE SPATIAL SIGNAL OF MPA EFFECTS 
#' 
#' 
#########################################################

# LIBRARY PACKAGES

library(rnaturalearthdata)
library(spatialRF)
library(ranger)
library(geodist)
library(sjmisc )
library(dplyr)


# LOAD DATA

model_data <- readr::read_rds("data/model_data.rds")

# IMPORT AND ATTACH AESTHETIC DATA

aesthetic_survey_data <- read.csv("outputs/survey_aesth.csv")

model_data <- merge(model_data, aesthetic_survey_data, by="SurveyID")

# SELECT VARIABLES AND PRE-AVERAGE BY SITECODE

rf_data <- model_data %>%
  select(aesthe_survey, HDI2017, fshD, gravtot2, SiteLongitude, SiteLatitude, SiteCode) %>%
  group_by(SiteCode) %>%
  summarise_all(.funs=mean)

site_MPA <- model_data %>%
  select(SiteCode, MPA) %>%
  distinct(SiteCode, MPA)

MPA_dummy <- to_dummy(site_MPA$MPA) %>%
  select(MPA_2)

MPA_dummy <- data.frame(SiteCode=site_MPA$SiteCode, 
                        MPA_No_Take=MPA_dummy$MPA_2)

rf_data <- merge(rf_data, MPA_dummy, by="SiteCode")

# LOG TRANSFORM AESTHETIC VALUE

rf_data$log_aesth <- log(rf_data$aesthe_survey)

# CHECK FOR DATA REQUIREMENTS

sapply(rf_data, function(x) sum(is.na(x)))

rf_data <- na.omit(rf_data)

apply(rf_data, 2, var) == 0

# DEFINE SPATIAL COMPONENTS

xy <- rf_data %>%
  select(SiteLongitude, SiteLatitude) %>%
  rename(x = SiteLongitude, y = SiteLatitude)

dist.matrix <- geodist(xy, measure = "geodesic")
dist.matrix <- dist.matrix/1000 # CONVERT TO KM 
range(as.numeric(dist.matrix))

distance.thresholds <- c(0, 2000, 4000, 8000, 16000)

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
    data = rf_data,
    ggplot2::aes(
      x = SiteLongitude,
      y = SiteLatitude,
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
  ggplot2::scale_x_continuous(limits = c(-180, 180)) +
  ggplot2::scale_y_continuous(limits = c(-90, 90))  +
  ggplot2::ggtitle("Global Distribution of Aesthetic Value") + 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

# LOOK AT SPATIAL AUTOCORRELATION

dependent.variable.name <- "log_aesth"
predictor.variable.names <- c("MPA_No_Take", "HDI2017", "fshD", "gravtot2")

spatialRF::plot_training_df_moran(
  data = rf_data,
  dependent.variable.name = dependent.variable.name,
  predictor.variable.names = predictor.variable.names,
  distance.matrix = dist.matrix,
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
  distance.matrix = dist.matrix
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

local.importance <- spatialRF::get_importance_local(spatial.model)

local.importance <- cbind(xy, local.importance)

# MAP IMPORTANCE OF MPA 

#colors
color.low <- viridis::viridis(
  3,
  option = "F"
)[2]
color.high <- viridis::viridis(
  3,
  option = "F"
)[1]

#plot of MPA_No_Take
p1 <- ggplot2::ggplot() +
  ggplot2::geom_sf(
    data = world,
    fill = "white"
  ) +
  ggplot2::geom_point(
    data = local.importance,
    ggplot2::aes(
      x = x,
      y = y,
      color = gravtot2,
    ) 
  ) +
  ggplot2::scale_x_continuous(limits = c(-180, 180)) +
  ggplot2::scale_y_continuous(limits = c(-90, 90)) +
  ggplot2::scale_color_gradient2(
    low = color.low, 
    high = color.high
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom") + 
  ggplot2::ggtitle("x") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    legend.key.width = ggplot2::unit(1,"cm")
  ) + 
  ggplot2::labs(color = "Importance") + 
  ggplot2::xlab("Longitude") + 
  ggplot2::ylab("Latitude")

p1

# RESPONSE CURVES 

spatialRF::plot_response_curves(
  spatial.model,
  quantiles = 0.5,
  ncol = 4
)
