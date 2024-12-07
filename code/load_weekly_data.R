# Load in SPRING-FISH weekly zooplankton and nutrient data and combine with profiler data
# AGS
# December 7, 2024

library(tidyverse)

# ------------------ Nutrient data ------------------

# Read in "nutrient_data" RDS object that was generated after some cleaning steps in my "SPRING-FISH" R project (nutrients_plotting.R)
nutrients <- readRDS(file = "data/nutrient_data.rds")

# Remove columns that I won't be using in this analysis
# and remove L01 because I won't be using the lake station in this analysis
nutrients_clean <- nutrients %>%
  select(c("enclosure", "depth_layer", "TP_ugL", "TN_ugL", "week")) %>% # (there are other nutrient data columns, but I am choosing just to use TN and TP)
  subset(enclosure != "L01")

# Visualize change in nutrients across the weeks, to assess collinearity and determine if I should keep nutrients in the model
# all mesocosms together - should I be doing this for each mesocosm separately? won't have error bars though
ggplot(nutrients_clean, aes(x = week, y = TP_ugL)) +
  geom_boxplot() +
  facet_wrap(~depth_layer) +
  theme_bw()
ggplot(nutrients_clean, aes(x = week, y = TN_ugL)) +
  geom_boxplot() +
  facet_wrap(~depth_layer) +
  theme_bw()
# sort of a relationship over time with meta TN...enough to justify removing it from the analysis??

# Join with profiler data
data1 <- left_join(profiler, nutrients_clean, by = c("enclosure", "depth_layer", "week"))

# ------------------ Zooplankton density data ------------------

# Read in "zooplankton_experiment_summary.csv" that was generated in my "SPF_analysis" R project
zoop_density <- read.csv("data/zooplankton_experiment_summary.csv")

# Remove columns that I won't be using in this analysis,
# remove L01 because I won't be using the lake station in this analysis,
# sum up the densities for all taxa within each week and enclosure,
# and make week a character
zoop_density_clean <- zoop_density %>%
  select(c("week", "enclosure", "taxon", "counts.L")) %>% 
  subset(enclosure != "L01") %>%
  group_by(week, enclosure) %>%
  summarize(total_zoop_density_L = sum(counts.L)) %>%
  ungroup() %>%
  mutate_at("week", as.character)

# Visualize change in zoop density across the weeks, to assess collinearity and determine if I should keep zoop density in the model
# all mesocosms together - should I be doing this for each mesocosm separately? won't have error bars though
ggplot(zoop_density_clean, aes(x = week, y = total_zoop_density_L)) +
  geom_boxplot() +
  theme_bw()

# Join with profiler data + nutrient data
data2 <- left_join(data1, zoop_density_clean, by = c("enclosure", "week"))

# ------------------ Zooplankton length data ------------------

# Read in "zooplankton_experiment_lengths.csv" that was generated in my "SPF_analysis" R project
zoop_lengths <- read.csv("data/zooplankton_experiment_lengths.csv")

# Remove columns that I won't be using in this analysis,
# remove L01 because I won't be using the lake station in this analysis,
# remove week 0,
# average the lengths for all taxa within each week and enclosure,
# and make week a character
zoop_lengths_clean <- zoop_lengths %>%
  select(c("week", "enclosure", "length_um")) %>% 
  subset(enclosure != "L01") %>%
  subset(week != 0) %>%
  group_by(week, enclosure) %>%
  summarize(avg_zoop_length_um = mean(length_um)) %>%
  ungroup() %>%
  mutate_at("week", as.character)

# Visualize change in zoop lengths across the weeks, to assess collinearity and determine if I should keep zoop density in the model
# all mesocosms together - should I be doing this for each mesocosm separately? won't have error bars though
ggplot(zoop_lengths_clean, aes(x = week, y = avg_zoop_length_um)) +
  geom_boxplot() +
  theme_bw()
# sort of a relationship over time with zoop lengths...enough to justify removing it from the analysis??

# Join with profiler data + nutrient data + zooplankton density data
data3 <- left_join(data2, zoop_lengths_clean, by = c("enclosure", "week"))

# ------------------ Create final dataframe for running EDA and models  ------------------

# Remove columns that I won't be using in this analysis
data_final <- data3 %>%
  select(-c(hour_profile.datetime, date, depth_layer, photosynthetically.active.radiation.up)) # assuming I'll use week as my AR(1) term

# Check the structure
str(data_final)

# Make enclosure, fish_treatment, week, and timing factors instead of characters
data_final$enclosure <- as.factor(data_final$enclosure)
data_final$fish_treatment <- as.factor(data_final$fish_treatment)
data_final$week <- as.factor(data_final$week)
data_final$timing <- as.factor(data_final$timing)

# Check the structure again
str(data_final)