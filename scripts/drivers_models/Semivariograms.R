# Load required packages
library(gstat)
library(sp)
library(spdep)
library(sf)
library(Matrix)
library(glmmTMB)
library(brms)
library(mgcv)
library(dplyr)
library(gridExtra)
library(grid)
library(lattice)

# VARIOGRAM WITH RAW DATA

standardized_data <-  readRDS("outputs/standardized_data.rds")

raw_coords <- standardized_data %>%
  select(SiteCode, SiteLongitude, SiteLatitude, aesthe_survey_abund) %>%
  group_by(SiteCode) %>%
  summarise_all(.funs=mean) %>%
  ungroup() %>%
  mutate(log_aes = log(aesthe_survey_abund)) %>%
  na.omit() %>%
  distinct(SiteCode, SiteLongitude, SiteLatitude, log_aes)

# Convert to spatial object
raw_sf <- st_as_sf(raw_coords, coords = c("SiteLongitude", "SiteLatitude"), crs = 4326)

# Compute empirical semivariogram
raw_variogram <- variogram(log_aes ~ 1, data = raw_sf, cutoff=8000)

# Fit a theoretical semivariogram model
raw_var_model <- fit.variogram(raw_variogram, model = vgm("Sph"))

# Plot the empirical and fitted semivariogram
plot(raw_variogram, raw_var_model, main = "Semivariogram", 
     xlim=c(0,8000))

# Print model parameters
print(raw_var_model)





# VARIGRAOM FROM INTERCEPT ONLY MODEL

model_data <- standardized_data %>%
  select(SiteCode, Country, SiteLongitude, SiteLatitude, aesthe_survey_abund) %>%
  na.omit()

int_model <- brm(log(aesthe_survey_abund) ~ (1 | Country/SiteCode),
                 family=gaussian(),
                 chains=4, iter=4000, cores=ncores,
                 threads = threading(8),
                 family=gaussian(),
                 data=model_data)

int_mod_coords <- data.frame(SiteLongitude=model_data$SiteLongitude,
                        SiteLatitude=model_data$SiteLatitude,
                        residuals=resid(int_model))

# Convert to spatial object
int_model_sf <- st_as_sf(int_mod_coords, coords = c("SiteLongitude", "SiteLatitude"), crs = 4326)

# Compute empirical semivariogram
int_model_variogram <- variogram(residuals.Estimate ~ 1, data = int_model_sf, cutoff=8000)

# Fit a theoretical semivariogram model
int_model_var_model <- fit.variogram(int_model_variogram, model = vgm("Sph"))

# Plot the empirical and fitted semivariogram
plot(int_model_variogram, int_model_var_model, main = "Semivariogram", 
     xlim=c(0,8000))

# Print model parameters
print(int_model_var_model)





# VARIGRAOM FROM SPATIAL MODEL

model_data <- standardized_data %>%
  select(SiteCode, Country, SiteLongitude, SiteLatitude, aesthe_survey_abund) %>%
  na.omit()

spa_model <- brm(log(aesthe_survey_abund) ~ 
                       s(SiteLatitude, SiteLongitude, bs = "gp") +
                       (1 | Country/SiteCode),
                 control = list(max_treedepth = 15),
                 chains=4, iter=4000, cores=ncores,
                 family=gaussian(),
                 data=model_data)
saveRDS(spa_model, "outputs/BIG_FILES/spa_model.rds")

spa_mod_coords <- data.frame(SiteLongitude=model_data$SiteLongitude,
                         SiteLatitude=model_data$SiteLatitude,
                         residuals=resid(spa_model))

# Convert to spatial object
spa_model_sf <- st_as_sf(spa_mod_coords, coords = c("SiteLongitude", "SiteLatitude"), crs = 4326)

# Compute empirical semivariogram
spa_model_variogram <- variogram(residuals.Estimate ~ 1, data = spa_model_sf, cutoff=8000)

# Fit a theoretical semivariogram model
spa_model_var_model <- fit.variogram(spa_model_variogram,
                                     model = vgm("Sph"))

# Plot the empirical and fitted semivariogram
plot(spa_model_variogram, spa_model_var_model, main = "Semivariogram", xlim=c(0,8000))

# Print model parameters
print(spa_model_var_model)




# par(mfrow=c(1,3))
# 
# # Plot the raw semivariogram
# p1 <- plot(raw_variogram, pch=19,
#            #raw_var_model, 
#            xlab="Distance (km)", 
#            main = "Semivariogram of Raw Values", 
#            xlim=c(0,8000))
# mtext("A", side = 3, line = -1.5, adj = 0, cex = 1.5)
# 
# # Plot the model residual semivariogram
# p2 <- plot(int_model_variogram, pch=19,
#            #int_model_var_model,
#            xlab="Distance (km)", main = "Semivariogram of Model Residuals", xlim=c(0,8000))
# 
# # Plot the semivariogram
# p3 <- plot(spa_model_variogram, pch=19,
#            #spa_model_var_model,xlab="Distance (km)",
#            main = "Semivariogram of Model Residuals", xlim=c(0,8000))
# 
# gridExtra::grid.arrange(p1, p2, p3, ncol=3)




common_settings <- list(par.main.text = list(cex = 1.1, font = 1))  

p1 <- plot(raw_variogram, main = "Semivariogram of Raw Values\n ", 
           xlab = "Distance (km)", ylab="Semivariance",
           xlim = c(0,8000), pch=19,
           par.settings = common_settings)

p2 <- plot(int_model_variogram, main = "Semivariogram of Model Residuals\n(Hierarchical Random Effects Only)", 
           xlab = "Distance (km)", ylab="Semivariance",
           xlim = c(0,8000), pch=19,
           par.settings = common_settings)

p3 <- plot(spa_model_variogram, main = "Semivariogram of Model Residuals\n(Hierearchical + Spatial Random Effects)", 
           xlab = "Distance (km)", ylab="Semivariance",
           xlim = c(0,8000), pch=19,
           par.settings = common_settings)

grid.arrange(
  arrangeGrob(p1, top = textGrob("A", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))),
  arrangeGrob(p2, top = textGrob("B", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))),
  arrangeGrob(p3, top = textGrob("C", x = unit(0, "npc"), y = unit(1, "npc"), just = c("left", "top"), gp = gpar(fontsize = 16, fontface = "bold"))),
  ncol = 3
)
