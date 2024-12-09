# EDA of SPRING-FISH profiler data and covariates
# AGS
# December 7, 2024

library(tidyverse)
library(GGally)

summary(data_final)
str(data_final)

# check for NAs
is.na(data_final)

# Look at the distribution of our response variable, chlorophylla.a
hist(data_final$chlorophyll.a)
# all values are above 0, continuous, slightly skewed
# will use the gamma(?) distribution

# Preliminary look at the the effects of each covariate on chlorophyll-a
plot(data_final$water.temperature, data_final$specified.depth)
plot(data_final$oxygen.concentration, data_final$chlorophyll.a)
plot(data_final$TP_ugL, data_final$chlorophyll.a)
plot(data_final$TN_ugL, data_final$chlorophyll.a)
plot(data_final$total_zoop_density_L, data_final$chlorophyll.a)
plot(data_final$avg_zoop_length_um, data_final$chlorophyll.a)

plot(data_final$fish_treatment, data_final$chlorophyll.a)
plot(data_final$timing, data_final$chlorophyll.a)

plot(data_final$enclosure, data_final$chlorophyll.a)
plot(data_final$specified.depth, data_final$chlorophyll.a)
plot(data_final$week, data_final$chlorophyll.a)

# Assess collinearity between covariates
ggpairs(data_final)
# water temperature and depth are highly correlated (-0.8), is that okay? to be expected obviously

plot(data_final$total_zoop_density_L, data_final$avg_zoop_length_um)

# should I z-transform (scale) all the continuous data?
data_final <- data_final %>% mutate(water.temperature.z = scale(water.temperature)[,1],
                                    oxygen.concentration.z = scale(oxygen.concentration)[,1],
                                    TP_ugL.z = scale(TP_ugL)[,1],
                                    TN_ugL.z = scale(TN_ugL)[,1],
                                    total_zoop_density_L.z = scale(total_zoop_density_L)[,1],
                                    avg_zoop_length_um.z = scale(avg_zoop_length_um)[,1])
