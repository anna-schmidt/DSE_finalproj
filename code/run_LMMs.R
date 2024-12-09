# Running LMMs on SPRING-FISH profiler data and covariates
# AGS
# December 7, 2024

library(tidyverse)
library(glmmTMB)
library(DHARMa)

# Model without any random effects
mod1 <- glmmTMB(chlorophyll.a ~ water.temperature.z + oxygen.concentration.z +
                     fish_treatment + timing +
                     TP_ugL.z + TN_ugL.z + total_zoop_density_L.z + avg_zoop_length_um.z,
                   data = data_final,
                family = Gamma(link = "log"))
summary(mod1)

# Mixed effects model with enclosure as a random effect but ignoring time and space
mod2 <- glmmTMB(chlorophyll.a ~ water.temperature.z + oxygen.concentration.z +
                  fish_treatment + timing +
                  TP_ugL.z + TN_ugL.z + total_zoop_density_L.z + avg_zoop_length_um.z +
                  (1|enclosure), # constant slope, varying intercepts
                REML = T,
                data = data_final,
                family = Gamma(link = "log"))
summary(mod2)

# Use AIC to compare which model is a better fit
AIC(mod1, mod2)
# mod2 is better

# DHARMa for mod2
simOut <- simulateResiduals(mod1, n = 250)
plot(simOut)
# pretty bad

# Mixed effects model with enclosure as a random effect but and ar1() terms for both time and space

# First, make week and specified.depth factors
data_final$week <- as.factor(data_final$week)
data_final$specified.depth <- as.factor(data_final$specified.depth)

mod3 <- glmmTMB(chlorophyll.a ~ oxygen.concentration.z + # removed temperature
                  fish_treatment + timing +
                  fish_treatment*timing +
                  TP_ugL.z + TN_ugL.z + total_zoop_density_L.z + avg_zoop_length_um.z +
                  ar1(week + 0 | enclosure) +           # temporal AR1 within each enclosure
                  ar1(specified.depth + 0 | enclosure), # spatial AR1 within each enclosure
                REML = F, # F when you're doing fixed effects model selection
                data = data_final,
                family = Gamma(link = "log"))
summary(mod3)

# Use AIC to compare which model is a better fit
AIC(mod1, mod2, mod3)

# DHARMa for mod3
simOut <- simulateResiduals(mod3, n = 250)
plot(simOut)
# not very good

# if specified.depth is a factor, how does it know which depths are next to each other? and same with week?

# Notes:
# The + 0 ensures that no intercept is included in the AR1 term (this is required syntax for the AR1 structure in glmmTMB)
#  Specifying two independent AR1 correlation structures, one for week (temporal) and one for specified.depth (spatial)
# The AR1 correlation structure in glmmTMB requires a grouping variable to define independent clusters where autocorrelation is assumed
# 1 as the grouping variable specifies that AR1 correlations are applied globally (i.e., across the whole dataset)
# The ar1 function in glmmTMB uses the levels of the factor to determine the temporal order. This avoids errors from non-sequential or improperly sorted numeric values