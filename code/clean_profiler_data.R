# Data cleaning of SPRING-FISH profiler data
# AGS
# December 7, 2024

library(tidyverse)

# Read in "profiler_cleaned_final" RDS object that was generated after some cleaning steps in my "SPF_profiler" R project
profiler_cleaned_final <- readRDS(file = "data/profiler_cleaned_final.rds")

# Make a column with just the date,
# remove columns that I won't be using in this analysis, 
# subset to only use data from midnight to deal with non-photochemical quenching of fluorescence sensors,
# remove L01 because I won't be using the lake station in this analysis,
# add a column for week of the experiment,
# add a column for "before" or "after" the fish were added,
# add a column for depth_layer to allow for later joining with the nutrient data,
# and make date and week columns characters for potting
profiler <- profiler_cleaned_final %>%
  mutate(date = date(profile.datetime)) %>%
  select(-c(id, profile.datetime, conductivity, pH.value, phycocyanin, photosynthetically.active.radiation.up, probe.serial.number, day_profiledatetime, month_profiledatetime, year_profiledatetime)) %>%
  subset(hour_profile.datetime == 0) %>%
  subset(enclosure != "L01") %>%
  mutate(week = case_when(date %in% c("2023-04-25", "2023-04-26", "2023-04-27", "2023-04-28", "2023-04-29") ~ 1,
                          date %in% c("2023-04-30", "2023-05-01", "2023-05-02", "2023-05-03", "2023-05-04", "2023-05-05", "2023-05-06") ~ 2,
                          date %in% c("2023-05-07", "2023-05-08", "2023-05-09", "2023-05-10", "2023-05-11", "2023-05-12", "2023-05-13") ~ 3,
                          date %in% c("2023-05-14", "2023-05-15", "2023-05-16", "2023-05-17", "2023-05-18", "2023-05-19", "2023-05-20") ~ 4,
                          date %in% c("2023-05-21", "2023-05-22", "2023-05-23", "2023-05-24", "2023-05-25", "2023-05-26", "2023-05-27") ~ 5,
                          date %in% c("2023-05-28", "2023-05-29", "2023-05-30", "2023-05-31", "2023-06-01", "2023-06-02") ~ 6)) %>%
  mutate(timing = case_when(date %in% c("2023-04-25", "2023-04-26", "2023-04-29", "2023-04-30", "2023-05-01", "2023-05-02", "2023-05-03", "2023-05-04") ~ "before",
                            .default = "after")) %>%
  mutate(depth_layer = case_when(specified.depth <= 6 ~ "epi",
                                 .default = "meta"))  %>%
  mutate_at(c("week", "date"), as.character)

# Check str
str(profiler)

# Visualize vertical depth plots of all mesocosms, all dates, all depths to show how many data points we have
# could also make heatmaps of chlorophyll for each mesocosm, but would be harder to include as a single figure in the report, so maybe use this one
ggplot(profiler, aes(x = chlorophyll.a, y = specified.depth, color = date)) +
  geom_point(alpha = 0.4) +
  geom_line(orientation = "y", alpha = 0.4, aes(group = date)) +
  scale_y_reverse() +
  facet_wrap(~enclosure, ncol = 7) +
  theme_bw()
