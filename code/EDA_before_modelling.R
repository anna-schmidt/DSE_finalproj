# EDA of SPRING-FISH profiler data and covariates
# AGS
# December 7, 2024

library(tidyverse)
library(GGally)

summary(data_final)
str(data_final)

# Look at the distribution of our response variable, chlorophylla.a
hist(data_final$chlorophyll.a)
# all values are above 0, continuous, slightly right-skewed
# will use the gamma distribution

# Preliminary look at the the effects of each covariate on chlorophyll-a
plot(data_final$water.temperature, data_final$chlorophyll.a)
plot(data_final$oxygen.concentration, data_final$chlorophyll.a)
plot(data_final$TP_ugL, data_final$chlorophyll.a)
plot(data_final$TN_ugL, data_final$chlorophyll.a)
plot(data_final$total_zoop_biomass_ug, data_final$chlorophyll.a)
plot(data_final$avg_zoop_length_um, data_final$chlorophyll.a)

plot(data_final$fish_treatment, data_final$chlorophyll.a)
plot(data_final$timing, data_final$chlorophyll.a)

plot(data_final$enclosure, data_final$chlorophyll.a)
plot(data_final$specified.depth, data_final$chlorophyll.a)
plot(data_final$day, data_final$chlorophyll.a)

# Assess collinearity between covariates
ggpairs(data_final)
plot(data_final$water.temperature, data_final$specified.depth)
# water temperature and depth are highly correlated (-0.8) (to be expected) - decision to remove water.temperature as a covariate because it will be dealt with in the specified.depth AR term

# Models ran better if I removed the deep depths that were only in some of the mesocosms. Set the cutoff to <16 m to deal with this.
#data_final <- data_final %>% subset(specified.depth < 10) # this makes DHARMa better but is removing a of of the data....

# Z-transform (scale) all the numeric covariates
# and remove columns I won't be using anymore
data_final <- data_final %>% 
  mutate(oxygen.concentration.z = scale(oxygen.concentration)[,1],
         TP_ugL.z = scale(TP_ugL)[,1],
         TN_ugL.z = scale(TN_ugL)[,1],
         total_zoop_biomass_ugL.z = scale(total_zoop_biomass_ugL)[,1],
         avg_zoop_length_um.z = scale(avg_zoop_length_um)[,1]) %>%
  select(-c(oxygen.concentration,
            TP_ugL,
            TN_ugL,
            total_zoop_biomass_ugL,
            avg_zoop_length_um,
            water.temperature))



