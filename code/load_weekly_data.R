# Load in SPRING-FISH weekly zooplankton and nutrient data and combine with profiler data
# AGS
# December 7, 2024

library(tidyverse)
library(patchwork)

# ------------------ Nutrient data ------------------

# Read in "nutrient_data" RDS object that was generated after some cleaning steps in my "SPRING-FISH" R project (nutrients_plotting.R)
nutrients <- readRDS(file = "data/nutrient_data.rds")

# Remove columns that I won't be using in this analysis
# and remove L01 because I won't be using the lake station in this analysis
nutrients_clean <- nutrients %>%
  select(c("enclosure", "depth_layer", "TP_ugL", "TN_ugL", "week")) %>% # (there are other nutrient data columns, but I am choosing just to use TN and TP)
  subset(enclosure != "L01")

# Join with profiler data
data1 <- left_join(profiler, nutrients_clean, by = c("enclosure", "depth_layer", "week"))

# ------------------ Zooplankton biomass data ------------------

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

# Join with profiler data + nutrient data
data2 <- left_join(data1, zoop_biomass_clean, by = c("enclosure", "week"))

# ------------------ Assess changes in numeric covariates across time to justify inclusion in the model ------------------

# Visualize change in these parameters across the weeks, to assess collinearity and determine if I should keep them in the model

my_colors <- c("#E58606", "#99C945", "#5D69B1", "#CC61B0", "#A3D5FF","#24796C", "#DAA51B", "#2F8AC4", "#764E9F", "#dcbeff",  "#ED645A", "#a9a9a9", "#8D021F", "#000075")

# Oxygen concentration
p2 <- ggplot(profiler, aes(x = week, y = oxygen.concentration, fill = enclosure)) +
  geom_boxplot() +
  theme_bw(base_size = 15) +
  labs(x = "Week", y = "Oxygen concentration (mg/L)", color = "Enclosure") +
  scale_fill_manual(values = my_colors) + theme(legend.position = "top")

# Total phosphorus
p3 <- ggplot(nutrients_clean, aes(x = week, y = TP_ugL, fill = enclosure)) +
  geom_boxplot() +
  theme_bw(base_size = 15) +
  labs(x = "Week", y = "Total phosphorus concentration (\u00b5g/L)", color = "Enclosure") +
  scale_fill_manual(values = my_colors) + theme(legend.position = "none")

# Total nitrogen
p4 <- ggplot(nutrients_clean, aes(x = week, y = TN_ugL, fill = enclosure)) +
  geom_boxplot() +
  theme_bw(base_size = 15) +
  labs(x = "Week", y = "Total nitrogen concentration (\u00b5g/L)", color = "Enclosure") +
  scale_fill_manual(values = my_colors) + theme(legend.position = "none")

# Zooplankton biomass
p5 <- ggplot(zoop_biomass_clean, aes(x = week, y = total_zoop_biomass_ugL, color = enclosure)) +
  geom_boxplot() +
  theme_bw(base_size = 15) +
  labs(x = "Week", y = "Total zooplankton biomass (\u00b5g/L)", color = "Enclosure") +
  scale_color_manual(values = my_colors) + theme(legend.position = "none")

p6 <- (p2+p4) / (p3+p5)
ggsave(path = "figures", filename = "num_byweek_plots.png", plot = p6, device = "png", width = 48, height = 29, units = "cm", dpi = 300)

# ------------------ Create final dataframe for running EDA and models  ------------------

# Make a new column for "day" to use day as my AR(1) term instead of week (based on Mark's advice) (probably an easier way to do this but I couldn't figure it out, so just doing the long way)
# and remove columns that I won't be using in this analysis
data_final <- data2 %>%
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
  tidyr::expand(enclosure, full_seq(day, 1), specified.depth) %>%
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