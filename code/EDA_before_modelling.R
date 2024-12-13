# EDA of SPRING-FISH profiler data and covariates
# AGS
# December 7, 2024

library(tidyverse)
library(GGally)

summary(data_final)
str(data_final)

# Look at the distribution of our response variable, chlorophylla.a
p6 <- ggplot(data_final, aes(x = chlorophyll.a)) +
  geom_histogram(binwidth = 0.5, color = "white", fill = "#5D69B1") +
  theme_bw(base_size = 15) +
  labs(x = "Chlorophyll-a concentration (\u00b5g/L)", y = "Frequency")
ggsave(path = "figures", filename = "chlorophyll_histogram.png", plot = p6, device = "png", width = 17, height = 19, units = "cm", dpi = 300)
# all values are above 0, continuous, right-skewed
# will use the gamma distribution

# Preliminary look at the the effects of each covariate on chlorophyll-a
plot(data_final$water.temperature, data_final$chlorophyll.a)
plot(data_final$oxygen.concentration, data_final$chlorophyll.a)
plot(data_final$TP_ugL, data_final$chlorophyll.a)
plot(data_final$TN_ugL, data_final$chlorophyll.a)
plot(data_final$total_zoop_biomass_ug, data_final$chlorophyll.a)

plot(data_final$fish_treatment, data_final$chlorophyll.a)
plot(data_final$timing, data_final$chlorophyll.a)

plot(data_final$enclosure, data_final$chlorophyll.a)
plot(data_final$specified.depth, data_final$chlorophyll.a)
plot(data_final$day, data_final$chlorophyll.a)

# Assess collinearity between covariates
# want to do this on the data frame without NAs
data_final_noNAs <- data_final %>% subset(chlorophyll.a > 0)
p7 <- ggpairs(data_final_noNAs)
ggsave(path = "figures", filename = "collinearity.png", plot = p7, device = "png", width = 30, height = 30, units = "cm", dpi = 300)

# Plot of water temperature and depth correlation
p8 <- ggplot(data_final, aes(x = water.temperature, y = specified.depth)) +
  scale_y_reverse() +
  geom_point(alpha = 1) +
  theme_bw(base_size = 15) +
  labs(x = expression(Temperature~"("*degree*C*")"), y = "Depth (m)")
ggsave(path = "figures", filename = "temp_vs_depth.png", plot = p8, device = "png", width = 17, height = 17, units = "cm", dpi = 300)
# water temperature and depth are highly correlated (-0.8) (to be expected) - decision to remove water.temperature as a covariate because it will be dealt with in the specified.depth AR term

# Z-transform (scale) all the numeric covariates
# and remove columns I won't be using anymore
data_final <- data_final %>% 
  mutate(oxygen.concentration.z = scale(oxygen.concentration)[,1],
         TP_ugL.z = scale(TP_ugL)[,1],
         TN_ugL.z = scale(TN_ugL)[,1],
         total_zoop_biomass_ugL.z = scale(total_zoop_biomass_ugL)[,1]) %>%
  select(-c(oxygen.concentration, TP_ugL, TN_ugL, total_zoop_biomass_ugL, water.temperature))