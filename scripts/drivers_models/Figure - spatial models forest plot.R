

####################################################################
## CODE TO RUN MODELS WITH GUASSIAN PROCESS SPATIAL RANDOM EFFECTS
## USING THE GLMBB TMB PACKAGE
####################################################################

library(ggplot2)
library(brms)
library(tidyverse)

# IMPORT ALL MODEL OUTPUTS
original_output <- readRDS("outputs/dag_output.rds")
inla_output <- readRDS("outputs/inla_output.rds")

original_output <- as.data.frame(posterior_summary(original_output)) %>%
  dplyr::select(Estimate, Q2.5, Q97.5) %>%
  rename("CI_2.5" = Q2.5,
         "CI_97.5" = Q97.5) %>%
  rownames_to_column("variable") %>%
  mutate(model = rep("original")) %>%
  dplyr::select(variable, Estimate, CI_2.5, CI_97.5, model)

#############################
## COMBINE DATA FOR GGPLOT ##
#############################

plot_data <- rbind(original_output, inla_output)

plot_data$model <- factor(plot_data$model, 
                          levels=rev(c("original","SPDE")))

######################
## MAKE FOREST PLOT ##
######################

forest_plot <- ggplot(plot_data, aes(x = Estimate, y = variable, color=model)) +
  geom_vline(xintercept = 0, color = "grey", linetype = "solid", size = 0.8) +
  geom_pointrange(aes(xmin=CI_2.5, xmax=CI_97.5), size=0.75, linewidth=1,
                  position = position_dodge(width=0.5)) + # Add points for the estimates
  theme(
    strip.background = element_rect(color = "black", fill = "lightgrey", size = 1), # Add border to facet titles
    axis.text.y = element_text(size = 12), # Style predictor labels
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.text = element_text(size=12),
    legend.title = element_text(size=12),
    panel.border = element_rect(color = "black", fill=NA)
  ) +
  labs(
    x = "Coefficient Estimate",
    y = NULL,
    fill = "Variable Type",
    title = NULL
  ) 

forest_plot

ggsave("forest_plot_spatial.png", forest_plot,
                      device = 'png',
                      dpi = 300, 
                      units = 'in', 
                      width = 8, 
                      height = 6, 
                      bg = 'white')





