# Running LMMs on SPRING-FISH profiler data and covariates
# AGS
# December 7, 2024

library(tidyverse)
library(glmmTMB)

# figure out how to use a gamma distribution for the family in glmmTMB()

# Model without any random effects
mod1 <- glmmTMB(chlorophyll.a ~ water.temperature.z + oxygen.concentration.z +
                     fish_treatment + timing +
                     TP_ugL.z + TN_ugL.z + total_zoop_density_L.z + avg_zoop_length_um.z,
                   data = data_final,
                family = Gamma(link = "log"))
summary(mod1)

# Mixed effects model with enclosure as a random effect
mod2 <- glmmTMB(chlorophyll.a ~ water.temperature.z + oxygen.concentration.z +
                  fish_treatment + timing +
                  TP_ugL.z + TN_ugL.z + total_zoop_density_L.z + avg_zoop_length_um.z +
                  (1|enclosure),
                REML = T,
                data = data_final,
                family = Gamma(link = "log"))
summary(mod2)

# Use AIC to compare which model is a better fit
AIC(mod1, mod2)

# Mixed effects model with enclosure as a random effect and week and depth as AR(1) terms
