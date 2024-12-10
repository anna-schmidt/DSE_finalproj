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
  theme_bw()
ggplot(nutrients_clean, aes(x = week, y = TN_ugL)) +
  geom_boxplot() +
  theme_bw()

# Join with profiler data
data1 <- left_join(profiler, nutrients_clean, by = c("enclosure", "depth_layer", "week"))

# ------------------ Zooplankton density data ------------------

# Read in "zooplankton_experiment_summary.csv" that was generated in my "SPF_analysis" R project
zoop_density <- read.csv("data/zooplankton_experiment_summary.csv")
# Read in LM-table that I made for calculating biomass from lengths for each taxon
zoop_LMtable <- read.csv("data/zooplankton_LMtable.csv")

# Calculate biomass for each taxon in each sample
zoop_density_biomass <- zoop_density %>%
  # first, fix taxa that are mis-spelled in some samples
  mutate(taxon = str_replace(taxon, "Harpacticoida", "Harpacticoid"),
         taxon = str_replace(taxon, "Alona Spp.", "Alona spp."),
         taxon = str_replace(taxon, "D. longispinpa", "D. longispina"),
         taxon = str_replace(taxon, "Daphnia longispina", "D. longispina")) %>%
  left_join(zoop_LMtable, by = "taxon") %>%
  mutate(biomass_ug_ind = exp(b*log(avg_length_mm) + lnA)) %>%
  mutate(biomass_ug_L = biomass_ug_ind*counts.L)
  

# Remove columns that I won't be using in this analysis,
# remove L01 because I won't be using the lake station in this analysis,
# sum up the densities for all taxa within each week and enclosure,
# and make week a character
zoop_biomass_clean <- zoop_density_biomass %>%
  select(c("week", "enclosure", "taxon", "biomass_ug_L")) %>% 
  subset(enclosure != "L01") %>%
  group_by(week, enclosure) %>%
  summarize(total_zoop_biomass_ugL = sum(biomass_ug_L)) %>%
  ungroup() %>%
  mutate_at("week", as.character)

# Visualize change in zoop density across the weeks, to assess collinearity and determine if I should keep zoop density in the model
# all mesocosms together - should I be doing this for each mesocosm separately? won't have error bars though
ggplot(zoop_biomass_clean, aes(x = week, y = total_zoop_biomass_ugL)) +
  geom_boxplot() +
  facet_wrap(~enclosure) +
  theme_bw()

# Join with profiler data + nutrient data
data2 <- left_join(data1, zoop_biomass_clean, by = c("enclosure", "week"))

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

# Make a new column for "day" to use day as my AR(1) term instead of week (based on Mark's advice)
# and remove columns that I won't be using in this analysis
data_final <- data3 %>%
  mutate(day = case_when(date == "2023-04-25" ~ 1,
                         date == "2023-04-26" ~ 2,
                         date == "2023-04-29" ~ 5,
                         date == "2023-04-30" ~ 6,
                         date == "2023-05-01" ~ 7,
                         date == "2023-05-02" ~ 8,
                         date == "2023-05-03" ~ 9,
                         date == "2023-05-04" ~ 10,
                         date == "2023-05-05" ~ 11,
                         date == "2023-05-06" ~ 12,
                         date == "2023-05-07" ~ 13,
                         date == "2023-05-08" ~ 14,
                         date == "2023-05-09" ~ 15,
                         date == "2023-05-10" ~ 16,
                         date == "2023-05-11" ~ 17,
                         date == "2023-05-12" ~ 18,
                         date == "2023-05-13" ~ 19,
                         date == "2023-05-14" ~ 20,
                         date == "2023-05-15" ~ 21,
                         date == "2023-05-16" ~ 22,
                         date == "2023-05-17" ~ 23,
                         date == "2023-05-18" ~ 24,
                         date == "2023-05-19" ~ 25,
                         date == "2023-05-20" ~ 26,
                         date == "2023-05-21" ~ 27,
                         date == "2023-05-22" ~ 28,
                         date == "2023-05-23" ~ 29,
                         date == "2023-05-24" ~ 30,
                         date == "2023-05-25" ~ 31,
                         date == "2023-05-26" ~ 32,
                         date == "2023-05-27" ~ 33,
                         date == "2023-05-28" ~ 34,
                         date == "2023-05-29" ~ 35,
                         date == "2023-05-30" ~ 36,
                         date == "2023-05-31" ~ 37,
                         date == "2023-06-01" ~ 38,
                         date == "2023-06-02" ~ 39)) %>%
  select(-c(hour_profile.datetime, date, week, depth_layer))

# For AR1 terms, need to use expand_grid to fill in NAs for dates and specified.depths that are missing

# Crate a dataframe of every combination of enclosure, day (filling in gaps), and specified.depth
empty_names <- data_final %>%
  expand(enclosure, full_seq(day, 1), specified.depth) %>%
  rename(day = `full_seq(day, 1)`)

# Join on those columns
data_final <- data_final %>%
  full_join(empty_names)

# Now, there are NAs everywhere there is missing data, but all combinations of specified.depth and day are here for each enclosure.

# Check the structure
str(data_final)

# Make enclosure, fish_treatment, and timing factors instead of characters
data_final$enclosure <- as.factor(data_final$enclosure)
data_final$fish_treatment <- as.factor(data_final$fish_treatment)
data_final$timing <- as.factor(data_final$timing)

# Check the structure again
str(data_final)