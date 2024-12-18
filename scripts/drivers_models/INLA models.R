

####################################################################
## CODE TO RUN MODELS WITH GAUSSIAN PROCESS RANDOM EFFECTS
## USING STOCHASTIC PARTIAL DIFFERENTIAL EQUATIONS
## USING THE INLA PACKAGE
####################################################################

library(ggplot2)
library(INLA)
library(sp)
library(gstat)
library(sf)
library(tidyverse)

######################################
## LOAD DATA AND DEFINE COORDINATES ##
######################################

standardized_data <-  read_rds("outputs/standardized_data.rds")
coords <- cbind(standardized_data$SiteLongitude, standardized_data$SiteLatitude)

########################################################
## SET MAX EDGE BASED ON STANDARD PRACTICE 
## MAX.EDGE = 1/10 OF TOTAL RANGE OF STUDY AREA
## BASED ON DIAGONAL DISTANCE
## AND CUTOFF = 1/5 OF MAX.EDGE
########################################################

max.edge<- max(range(dist(coords))/10)
cutoff = max.edge/5

# Create the mesh
mesh <- inla.mesh.2d(
  loc = coords,
  cutoff = cutoff,              # Minimum distance between nodes
  max.edge = c(1,5)*max.edge,       # Maximum edge length
  offset = c(max.edge/2, max.edge/2)         # Extend mesh beyond data
)

# Plot the mesh
plot(mesh, main = "Spatial Mesh")
points(coords, col = "red", pch = 16)

spde <- inla.spde2.matern(
  mesh = mesh,
  alpha = 2 # Controls smoothness of the spatial random field
)

# Project spatial locations onto the mesh
A <- inla.spde.make.A(mesh, loc = coords)

# Create an index for the spatial random effect
spatial_index <- inla.spde.make.index(name = "spatial", n.spde = spde$n.spde)


###############
## SST MODEL ##
###############

sst_stack <- inla.stack(
  data = list(y = log(standardized_data$aesthe_survey_abund)),
  A = list(A, 1),  # Spatial projector matrix and design matrix
  effects = list(
    list(spatial_field = 1:mesh$n),
    data.frame(intercept = 1,
               sst_mean = standardized_data$sst_mean,
               abs_latitude = standardized_data$abs_latitude,
               Temperature_Zone = as.factor(standardized_data$Temperature_Zone),
               Country = standardized_data$Country,
               SiteCode = standardized_data$SiteCode)
  )
)



# Define the formula
sst_formula <- y ~ sst_mean + abs_latitude + as.factor(Temperature_Zone) + 
  f(Country, model = "iid") +
  f(SiteCode, model = "iid") +
  f(spatial_field, model = spde)

# Fit the model
sst_result <- inla(
  sst_formula,
  data = inla.stack.data(sst_stack),
  control.predictor = list(A = inla.stack.A(sst_stack), compute = TRUE),
  family = "gaussian" # Use Gaussian, Poisson, or another distribution if appropriate
)

sst_inla_summary <- sst_result$summary.fixed  %>%
  select(mean,  "0.025quant", "0.975quant") %>%
  rename(CI_2.5 = "0.025quant",
         CI_97.5 = "0.975quant",
         Estimate = mean) %>%
  mutate(model = rep("SPDE")) %>%
  rownames_to_column("variable") %>%
  filter( variable == "sst_mean")

###################
## GRAVITY MODEL ##
###################

# Combine data and spatial effects into a stack
gravtot2_stack <- inla.stack(
  data = list(y = log(standardized_data$aesthe_survey_abund)),
  A = list(A, 1),  # Spatial projector matrix and design matrix
  effects = list(
    list(spatial_field = 1:mesh$n),
    data.frame(intercept = 1,
               gravtot2 = standardized_data$gravtot2,
               HDI2017 = standardized_data$HDI2017,
               abs_latitude = standardized_data$abs_latitude,
               Temperature_Zone = as.factor(standardized_data$Temperature_Zone),
               Country = standardized_data$Country,
               SiteCode = standardized_data$SiteCode)
  )
)


# Define the formula
gravtot2_formula <- y ~ gravtot2 + HDI2017 + abs_latitude + as.factor(Temperature_Zone) + 
  f(Country, model = "iid") +
  f(SiteCode, model = "iid") +
  f(spatial_field, model = spde)

# Fit the model
gravtot2_result <- inla(
  gravtot2_formula,
  data = inla.stack.data(gravtot2_stack),
  control.predictor = list(A = inla.stack.A(gravtot2_stack), compute = TRUE),
  family = "gaussian" # Use Gaussian, Poisson, or another distribution if appropriate
)

gravtot2_inla_summary <- gravtot2_result$summary.fixed  %>%
  select(mean,  "0.025quant", "0.975quant") %>%
  rename(CI_2.5 = "0.025quant",
         CI_97.5 = "0.975quant",
         Estimate = mean) %>%
  mutate(model = rep("SPDE")) %>%
  rownames_to_column("variable") %>%
  filter( variable == "gravtot2")


###############
## NPP MODEL ##
###############

# Combine data and spatial effects into a stack
NPP_stack <- inla.stack(
  data = list(y = log(standardized_data$aesthe_survey_abund)),
  A = list(A, 1),  # Spatial projector matrix and design matrix
  effects = list(
    list(spatial_field = 1:mesh$n),
    data.frame(intercept = 1,
               NPP_mean = standardized_data$NPP_mean,
               abs_latitude = standardized_data$abs_latitude,
               BO_nitrate = standardized_data$BO_nitrate,
               BO_phosphate = standardized_data$BO_phosphate,
               sst_mean = standardized_data$sst_mean,
               wave_energy = standardized_data$wave_energy,
               Temperature_Zone = as.factor(standardized_data$Temperature_Zone),
               Country = standardized_data$Country,
               SiteCode = standardized_data$SiteCode)
  )
)


# Define the formula
NPP_formula <- y ~ NPP_mean + abs_latitude + 
  BO_nitrate + BO_phosphate +
  sst_mean + wave_energy +
  as.factor(Temperature_Zone) + 
  f(Country, model = "iid") +
  f(SiteCode, model = "iid") +
  f(spatial_field, model = spde)

# Fit the model
NPP_result <- inla(
  NPP_formula,
  data = inla.stack.data(NPP_stack),
  control.predictor = list(A = inla.stack.A(NPP_stack), compute = TRUE),
  family = "gaussian" # Use Gaussian, Poisson, or another distribution if appropriate
)

NPP_inla_summary <- NPP_result$summary.fixed  %>%
  select(mean,  "0.025quant", "0.975quant") %>%
  rename(CI_2.5 = "0.025quant",
         CI_97.5 = "0.975quant",
         Estimate = mean) %>%
  mutate(model = rep("SPDE")) %>%
  rownames_to_column("variable") %>%
  filter( variable == "NPP_mean")


#################
## DEPTH MODEL ##
#################

# Combine data and spatial effects into a stack
Depth_stack <- inla.stack(
  data = list(y = log(standardized_data$aesthe_survey_abund)),
  A = list(A, 1),  # Spatial projector matrix and design matrix
  effects = list(
    list(spatial_field = 1:mesh$n),
    data.frame(intercept = 1,
               Depth = standardized_data$Depth,
               Temperature_Zone = as.factor(standardized_data$Temperature_Zone),
               Country = standardized_data$Country,
               SiteCode = standardized_data$SiteCode)
  )
)


# Define the formula
Depth_formula <- y ~ Depth + as.factor(Temperature_Zone) + 
  f(Country, model = "iid") +
  f(SiteCode, model = "iid") +
  f(spatial_field, model = spde)

# Fit the model
Depth_result <- inla(
  Depth_formula,
  data = inla.stack.data(Depth_stack),
  control.predictor = list(A = inla.stack.A(Depth_stack), compute = TRUE),
  family = "gaussian" # Use Gaussian, Poisson, or another distribution if appropriate
)

Depth_inla_summary <- Depth_result$summary.fixed  %>%
  select(mean,  "0.025quant", "0.975quant") %>%
  rename(CI_2.5 = "0.025quant",
         CI_97.5 = "0.975quant",
         Estimate = mean) %>%
  mutate(model = rep("SPDE")) %>%
  rownames_to_column("variable") %>%
  filter( variable == "Depth")


################
## FSHD MODEL ##
################

# Combine data and spatial effects into a stack
fshD_stack <- inla.stack(
  data = list(y = log(standardized_data$aesthe_survey_abund)),
  A = list(A, 1),  # Spatial projector matrix and design matrix
  effects = list(
    list(spatial_field = 1:mesh$n),
    data.frame(intercept = 1,
               fshD = standardized_data$fshD,
               HDI2017 = standardized_data$HDI2017,
               Temperature_Zone = as.factor(standardized_data$Temperature_Zone),
               Country = standardized_data$Country,
               SiteCode = standardized_data$SiteCode)
  )
)


# Define the formula
fshD_formula <- y ~ fshD + HDI2017 + as.factor(Temperature_Zone) + 
  f(Country, model = "iid") +
  f(SiteCode, model = "iid") +
  f(spatial_field, model = spde)

# Fit the model
fshD_result <- inla(
  fshD_formula,
  data = inla.stack.data(fshD_stack),
  control.predictor = list(A = inla.stack.A(fshD_stack), compute = TRUE),
  family = "gaussian" # Use Gaussian, Poisson, or another distribution if appropriate
)

fshD_inla_summary <- fshD_result$summary.fixed  %>%
  select(mean,  "0.025quant", "0.975quant") %>%
  rename(CI_2.5 = "0.025quant",
         CI_97.5 = "0.975quant",
         Estimate = mean) %>%
  mutate(model = rep("SPDE")) %>%
  rownames_to_column("variable") %>%
  filter( variable == "fshD")

################
## HDI MODEL ##
################

# Combine data and spatial effects into a stack
HDI2017_stack <- inla.stack(
  data = list(y = log(standardized_data$aesthe_survey_abund)),
  A = list(A, 1),  # Spatial projector matrix and design matrix
  effects = list(
    list(spatial_field = 1:mesh$n),
    data.frame(intercept = 1,
               HDI2017 = standardized_data$HDI2017,
               abs_latitude = standardized_data$abs_latitude,
               sst_mean = standardized_data$sst_mean,
               Temperature_Zone = as.factor(standardized_data$Temperature_Zone),
               Country = standardized_data$Country,
               SiteCode = standardized_data$SiteCode)
  )
)


# Define the formula
HDI2017_formula <- y ~ HDI2017 + abs_latitude + sst_mean + as.factor(Temperature_Zone) + 
  f(Country, model = "iid") +
  f(SiteCode, model = "iid") +
  f(spatial_field, model = spde)

# Fit the model
HDI2017_result <- inla(
  HDI2017_formula,
  data = inla.stack.data(HDI2017_stack),
  control.predictor = list(A = inla.stack.A(HDI2017_stack), compute = TRUE),
  family = "gaussian" # Use Gaussian, Poisson, or another distribution if appropriate
)

HDI2017_inla_summary <- HDI2017_result$summary.fixed  %>%
  select(mean,  "0.025quant", "0.975quant") %>%
  rename(CI_2.5 = "0.025quant",
         CI_97.5 = "0.975quant",
         Estimate = mean) %>%
  mutate(model = rep("SPDE")) %>%
  rownames_to_column("variable") %>%
  filter( variable == "HDI2017")


################
## MPA MODEL ##
################

# Combine data and spatial effects into a stack
MPA_stack <- inla.stack(
  data = list(y = log(standardized_data$aesthe_survey_abund)),
  A = list(A, 1),  # Spatial projector matrix and design matrix
  effects = list(
    list(spatial_field = 1:mesh$n),
    data.frame(intercept = 1,
               MPA = standardized_data$MPA,
               HDI2017 = standardized_data$HDI2017,
               gravtot2 = standardized_data$gravtot2,
               Temperature_Zone = as.factor(standardized_data$Temperature_Zone),
               Ecoregion = standardized_data$Ecoregion,
               Country = standardized_data$Country,
               SiteCode = standardized_data$SiteCode)
  )
)


# Define the formula
MPA_formula <- y ~  MPA + HDI2017 + gravtot2 + as.factor(Temperature_Zone) + 
  f(Country, model = "iid") +
  f(SiteCode, model = "iid") +
  f(spatial_field, model = spde)

# Fit the model
MPA_result <- inla(
  MPA_formula,
  data = inla.stack.data(MPA_stack),
  control.predictor = list(A = inla.stack.A(MPA_stack), compute = TRUE),
  family = "gaussian" # Use Gaussian, Poisson, or another distribution if appropriate
)

MPA_inla_summary <- MPA_result$summary.fixed %>%
  select(mean,  "0.025quant", "0.975quant") %>%
  rename(CI_2.5 = "0.025quant",
         CI_97.5 = "0.975quant",
         Estimate = mean)  %>%
  mutate(model = rep("SPDE")) %>%
  rownames_to_column("variable") %>%
  filter( variable == "MPANo take" | variable == "MPARestricted take")


#######################
## BENTHIC PC1 MODEL ##
#######################

# Combine data and spatial effects into a stack
PC1_stack <- inla.stack(
  data = list(y = log(standardized_data$aesthe_survey_abund)),
  A = list(A, 1),  # Spatial projector matrix and design matrix
  effects = list(
    list(spatial_field = 1:mesh$n),
    data.frame(intercept = 1,
               PC1_imputed = standardized_data$PC1_imputed,
               Depth = standardized_data$Depth,
               gravtot2 = standardized_data$gravtot2,
               MPA = standardized_data$MPA,
               NPP_mean = standardized_data$NPP_mean,
               sst_mean = standardized_data$sst_mean,
               Temperature_Zone = as.factor(standardized_data$Temperature_Zone),
               Country = standardized_data$Country,
               SiteCode = standardized_data$SiteCode)
  )
)


# Define the formula
PC1_formula <- y ~ PC1_imputed + Depth +
  gravtot2 + MPA + NPP_mean +sst_mean +
  as.factor(Temperature_Zone) + 
  f(Country, model = "iid") +
  f(SiteCode, model = "iid") +
  f(spatial_field, model = spde)

# Fit the model
PC1_result <- inla(
  PC1_formula,
  data = inla.stack.data(PC1_stack),
  control.predictor = list(A = inla.stack.A(PC1_stack), compute = TRUE),
  family = "gaussian" # Use Gaussian, Poisson, or another distribution if appropriate
)

PC1_inla_summary <- PC1_result$summary.fixed  %>%
  select(mean,  "0.025quant", "0.975quant") %>%
  rename(CI_2.5 = "0.025quant",
         CI_97.5 = "0.975quant",
         Estimate = mean) %>%
  mutate(model = rep("SPDE")) %>%
  rownames_to_column("variable") %>%
  filter( variable == "PC1_imputed")


#######################
## BENTHIC PC2 MODEL ##
#######################

# Combine data and spatial effects into a stack
PC2_stack <- inla.stack(
  data = list(y = log(standardized_data$aesthe_survey_abund)),
  A = list(A, 1),  # Spatial projector matrix and design matrix
  effects = list(
    list(spatial_field = 1:mesh$n),
    data.frame(intercept = 1,
               PC2_imputed = standardized_data$PC2_imputed,
               Depth = standardized_data$Depth,
               gravtot2 = standardized_data$gravtot2,
               MPA = standardized_data$MPA,
               NPP_mean = standardized_data$NPP_mean,
               sst_mean = standardized_data$sst_mean,
               Temperature_Zone = as.factor(standardized_data$Temperature_Zone),
               Country = standardized_data$Country,
               SiteCode = standardized_data$SiteCode)
  )
)


# Define the formula
PC2_formula <- y ~ PC2_imputed + Depth +
  gravtot2 + MPA + NPP_mean +sst_mean +
  as.factor(Temperature_Zone) + 
  f(Country, model = "iid") +
  f(SiteCode, model = "iid") +
  f(spatial_field, model = spde)

# Fit the model
PC2_result <- inla(
  PC2_formula,
  data = inla.stack.data(PC2_stack),
  control.predictor = list(A = inla.stack.A(PC2_stack), compute = TRUE),
  family = "gaussian" # Use Gaussian, Poisson, or another distribution if appropriate
)

PC2_inla_summary <- PC2_result$summary.fixed  %>%
  select(mean,  "0.025quant", "0.975quant") %>%
  rename(CI_2.5 = "0.025quant",
         CI_97.5 = "0.975quant",
         Estimate = mean) %>%
  mutate(model = rep("SPDE")) %>%
  rownames_to_column("variable") %>%
  filter( variable == "PC2_imputed")


#######################
## DHW MODEL ##
#######################

# Combine data and spatial effects into a stack
DHW_stack <- inla.stack(
  data = list(y = log(standardized_data$aesthe_survey_abund)),
  A = list(A, 1),  # Spatial projector matrix and design matrix
  effects = list(
    list(spatial_field = 1:mesh$n),
    data.frame(intercept = 1,
               dhw_mean = standardized_data$dhw_mean,
               abs_latitude = standardized_data$abs_latitude,
               sst_mean = standardized_data$sst_mean,
               Temperature_Zone = as.factor(standardized_data$Temperature_Zone),
               Country = standardized_data$Country,
               SiteCode = standardized_data$SiteCode)
  )
)

# Define the formula
DHW_formula <- y ~ dhw_mean + abs_latitude + sst_mean + 
  as.factor(Temperature_Zone) + 
  f(Country, model = "iid") +
  f(SiteCode, model = "iid") +
  f(spatial_field, model = spde)

# Fit the model
DHW_result <- inla(
  DHW_formula,
  data = inla.stack.data(DHW_stack),
  control.predictor = list(A = inla.stack.A(DHW_stack), compute = TRUE),
  family = "gaussian" # Use Gaussian, Poisson, or another distribution if appropriate
)

DHW_inla_summary <- DHW_result$summary.fixed  %>%
  select(mean,  "0.025quant", "0.975quant") %>%
  rename(CI_2.5 = "0.025quant",
         CI_97.5 = "0.975quant",
         Estimate = mean) %>%
  mutate(model = rep("SPDE")) %>%
  rownames_to_column("variable") %>%
  filter( variable == "dhw_mean")


##################################
## REGROUP ALL EFFECTS TOGETHER ##
##################################

inla_output <- rbind(sst_inla_summary, gravtot2_inla_summary,
                                 NPP_inla_summary, Depth_inla_summary,
                                 fshD_inla_summary, HDI2017_inla_summary,
                                 MPA_inla_summary, PC1_inla_summary,
                     PC2_inla_summary, DHW_inla_summary)

inla_output <- inla_output %>%
  mutate(variable = case_when(
    variable == "gravtot2" ~ "Log Human Gravity",
    variable == "MPANo take" ~ "No Take MPA",
    variable == "MPARestricted take" ~ "Restricted Take MPA",
    variable == "sst_mean" ~ "Sea Surface Temperature",
    variable == "NPP_mean" ~ "Net Primary Productivity",
    variable == "dhw_mean" ~ "Degree Heating Weeks",
    variable ==  "Depth" ~  "Depth",
    variable == "fshD" ~ "Fisheries Dependency",
    variable ==  "HDI2017" ~ "Human Dev. Index",
    variable == "PC1_imputed" ~ "Benthic Composition (PC1)",
    variable == "PC2_imputed" ~ "Benthic Composition (PC2)"))

saveRDS(inla_output, "outputs/inla_output.rds")



