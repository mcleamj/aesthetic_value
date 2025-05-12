

############################################################################
## CODE TO PLOT MODEL COEFFICIENTS FROM ORIGINAL AND WITH GAUSSIAN PROCESS
############################################################################

library(ggplot2)
library(brms)
library(tidyverse)

# IMPORT  MODEL OUTPUTS
original_output <- readRDS("outputs/dag_output.rds")

original_output <- as.data.frame(posterior_summary(original_output)) %>%
  dplyr::select(Estimate, Q2.5, Q97.5) %>%
  rename("CI_2.5" = Q2.5,
         "CI_97.5" = Q97.5) %>%
  rownames_to_column("variable") %>%
  mutate(model = rep("original")) 


spatial_outputs <- readRDS("outputs/BIG_FILES/spatial_tests/full_model_table.rds") %>%
  rename("CI_2.5" = conf.low,
         "CI_97.5" = conf.high,
         "Estimate" = estimate,
         "variable" = term) %>%
  filter(m_val == 2) %>%
  select(-m_val) %>%
  mutate(k_val = ifelse(k_val == "-1", "default", k_val)) %>%
  mutate(model = paste(rep("spatial (k="),k_val,")",sep="")) %>%
  select(-k_val)

#############################
## COMBINE DATA FOR GGPLOT ##
#############################

plot_data <- rbind(original_output, spatial_outputs)

unique(plot_data$model)

plot_data$model <- factor(plot_data$model, 
                          levels=rev(c("original" ,
                                       "spatial (k=10)" ,
                                       "spatial (k=default)",
                                       "spatial (k=50)" ,
                                       "spatial (k=100)" )))  

plot_data$model

######################
## MAKE FOREST PLOT ##
######################

# forest_plot <- ggplot(plot_data, aes(x = Estimate, y = variable, 
#                                      color=interaction(model, k_val),
#                       group = interaction(model, k_val))) +
#   geom_vline(xintercept = 0, color = "grey", linetype = "solid", size = 0.8) +
#   geom_pointrange(aes(xmin=CI_2.5, xmax=CI_97.5), size=0.75, linewidth=1,
#                   position = position_dodge(width=0.5)) + # Add points for the estimates
#   theme(
#     strip.background = element_rect(color = "black", fill = "lightgrey", size = 1), # Add border to facet titles
#     axis.text.y = element_text(size = 12), # Style predictor labels
#     axis.text.x = element_text(size = 12),
#     axis.title.x = element_text(size = 12),
#     legend.text = element_text(size=12),
#     legend.title = element_text(size=12),
#     panel.border = element_rect(color = "black", fill=NA)
#   ) +
#   theme_bw()+
#   labs(
#     x = "Coefficient Estimate",
#     y = NULL,
#     fill = "Variable Type",
#     title = NULL
#   ) 
# 
# forest_plot



library(ggplot2)

forest_plot <- ggplot(plot_data, 
                      aes(x = Estimate, y = variable, 
                          color = model))  +
  geom_vline(xintercept = 0, color = "grey40", linetype = "dashed", linewidth = 1) +
  geom_pointrange(
    aes(xmin = CI_2.5, xmax = CI_97.5),
    position = position_dodge(width = -0.75),  # Original model on top
    size = 0.8, linewidth = 1, fatten = 2
  ) +
  # Fixed color scale: Explicitly define values and labels
  scale_color_manual(
    values = c(
      "original" = "#1f77b4",  # Keep original model blue
      "spatial (k=10)"      = "#ffcccc",  # Very light red
      "spatial (k=default)" = "#ff9999",  # Light red
      "spatial (k=50)"      = "#ff6666",  # Medium red
      "spatial (k=100)"     = "#cc0000"   # Dark red
    ),
    labels = c(
      "original"           = "Original",
      "spatial (k=10)"     = "Spatial (k=10)",
      "spatial (k=default)" = "Spatial (k=Default)",
      "spatial (k=50)"     = "Spatial (k=50)",
      "spatial (k=100)"    = "Spatial (k=100)"
    )) +
  theme_bw(base_size = 14) +
  labs(x = "Coefficient Estimate", y = NULL) +
  theme(legend.position = "right")

forest_plot


ggsave(file=here::here("figures_tables/revised_submission/supplementary_figures",
                       "spatial_test_forest.png"), 
       forest_plot,
       device = 'png',
       dpi = 300, 
       units = 'in', 
       width = 8, 
       height = 6, 
       bg = 'white')





