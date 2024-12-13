# Running LMMs on SPRING-FISH profiler data and covariates
# AGS
# December 7, 2024

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(sjPlot)
library(viridis)
library(patchwork)
library(MuMIn)
set.seed(123)

# ------------------ Running models ------------------

# Model without any random effects
mod1 <- glmmTMB(chlorophyll.a ~ oxygen.concentration.z +
                     fish_treatment + timing +
                     fish_treatment*timing +
                     TP_ugL.z + TN_ugL.z + total_zoop_biomass_ugL.z,
                   data = data_final,
                family = Gamma(link = "log"))
summary(mod1)
tab_model(mod1, transform = NULL)

# Mixed effects model with enclosure as a random effect but ignoring time and space
mod2 <- glmmTMB(chlorophyll.a ~ oxygen.concentration.z +
                  fish_treatment + timing +
                  fish_treatment*timing +
                  TP_ugL.z + TN_ugL.z + total_zoop_biomass_ugL.z +
                  (1|enclosure), # constant slope, varying intercepts
                REML = T,
                data = data_final,
                family = Gamma(link = "log"))
summary(mod2)
tab_model(mod2, transform = NULL)

# Use AIC to compare which model is a better fit
AIC(mod1, mod2)
# mod2 is better

# DHARMa for mod2
simOut <- simulateResiduals(mod2, n = 250)
plot(simOut)

# Mixed effects model with enclosure as a random effect and ar1() terms for both time and space

# first, make day and specified.depth factors
data_final$day <- as.factor(data_final$day)
data_final$specified.depth <- as.factor(data_final$specified.depth)

mod3 <- glmmTMB(chlorophyll.a ~ oxygen.concentration.z +
                  fish_treatment + timing +
                  fish_treatment*timing +
                  TP_ugL.z + TN_ugL.z + total_zoop_biomass_ugL.z +
                  ar1(specified.depth + 0 | enclosure) + # spatial AR1 within each enclosure
                  ar1(day + 0 | enclosure),              # temporal AR1 within each enclosure
                REML = T, # T when you're still deciding on random effects structure
                data = data_final,
                family = Gamma(link = "log"))
summary(mod3)

# Use AIC to compare which model is a better fit
AIC(mod1, mod2, mod3)
# mod3 is the best -> move forward with our AR(1) structures

# ------------------ Model diagnostics ------------------

# Re-run mod3 with REML = F
mod3 <- glmmTMB(chlorophyll.a ~ oxygen.concentration.z +
                  fish_treatment + timing +
                  fish_treatment*timing +
                  TP_ugL.z + TN_ugL.z + total_zoop_biomass_ugL.z +
                  ar1(specified.depth + 0 | enclosure) + # spatial AR1 within each enclosure
                  ar1(day + 0 | enclosure),              # temporal AR1 within each enclosure
                REML = F, # F when you're doing fixed effects model selection
                data = data_final,
                family = Gamma(link = "log"))
summary(mod3)
tab_model(mod3, transform = NULL)

# Plot the model results
# oxygen concentration vs. chlorophyll-a
p9 <- plot_model(mod3, type = "pred", terms = "oxygen.concentration.z", title = NULL,
                 axis.title = c("Oxygen concentration", "Chlorophyll-a")) + ggtitle("") + theme_bw(base_size = 15)
# TN vs. chlorophyll-a
p10 <- plot_model(mod3, type = "pred", terms = "TN_ugL.z", title = NULL,
                  axis.title = c("Total nitrogen", "Chlorophyll-a")) + ggtitle("") + theme_bw(base_size = 15)
# TP vs. chlorophyll-a
p11 <- plot_model(mod3, type = "pred", terms = "TP_ugL.z", title = NULL,
                  axis.title = c("Total phosphorus", "Chlorophyll-a")) + ggtitle("") + theme_bw(base_size = 15)
p12 <- p9+p10+p11 + plot_annotation(tag_levels = 'A')
ggsave(path = "figures", filename = "model3_output.png", plot = p12, device = "png", width = 35, height = 18, units = "cm", dpi = 300)

# DHARMa for mod3
simOut <- simulateResiduals(mod3, n = 250)
plot(simOut)
# not very good

# Plot observed vs. predicted values for chlorophyll-a
# day and specified.depth need to be characters for this to work (see https://github.com/glmmTMB/glmmTMB/issues/602)
data_final$day <- as.character(data_final$day)
data_final$specified.depth <- as.character(data_final$specified.depth)
data_final$chlorophyll.a_predicted <- predict(mod3, newdata = data_final, type = "response")
data_final$day <- as.numeric(data_final$day)
data_final$specified.depth <- as.numeric(data_final$specified.depth)
p13 <- ggplot(data_final, aes(x = chlorophyll.a, y = chlorophyll.a_predicted, color = enclosure)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Observed chlorophyll-a", y = "Predicted chlorophyll-a", color = "Enclosure") +
  theme_bw(base_size = 15) +
  scale_color_manual(values = c("#E58606", "#99C945", "#5D69B1", "#CC61B0", "#A3D5FF","#24796C", "#DAA51B", "#2F8AC4", "#764E9F", "#dcbeff",  "#ED645A", "#a9a9a9", "#8D021F", "#000075"))
p14 <- ggplot(data_final, aes(x = chlorophyll.a, y = chlorophyll.a_predicted, color = day)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Observed chlorophyll-a", y = "Predicted chlorophyll-a", color = "Day") +
  theme_bw(base_size = 15) +
  scale_color_viridis()
p15 <- ggplot(data_final, aes(x = chlorophyll.a, y = chlorophyll.a_predicted, color = specified.depth)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Observed chlorophyll-a", y = "Predicted chlorophyll-a", color = "Depth") +
  theme_bw(base_size = 15) +
  scale_color_viridis()
p16 <- p13+p14+p15 + plot_annotation(tag_levels = 'A')
ggsave(path = "figures", filename = "obs_vs_pred.png", plot = p16, device = "png", width = 48, height = 16, units = "cm", dpi = 300)

# make day and specified.depth factors again
data_final$day <- as.factor(data_final$day)
data_final$specified.depth <- as.factor(data_final$specified.depth)

# Look for patterns in residuals across days and depths
# NAs need to be removed first
data_final_noNAs <- data_final %>% subset(chlorophyll.a > 0)
# extract the residuals
data_final_noNAs$residuals <- resid(mod3, type = "pearson")
p17 <- ggplot(data_final_noNAs, aes(x = day, y = residuals)) +
  geom_boxplot() +
  labs(x = "Day", y = "Residuals") +
  theme_bw(base_size = 15)
p18 <- ggplot(data_final_noNAs, aes(x = specified.depth, y = residuals)) +
  geom_boxplot() +
  labs(x = "Depth", y = "Residuals") +
  theme_bw(base_size = 15)
p19 <- p17 + p18 + plot_annotation(tag_levels = 'A')
ggsave(path = "figures", filename = "resid_patterns.png", plot = p19, device = "png", width = 46, height = 20, units = "cm", dpi = 300)

# Calculate goodness-of-fit metrics (marginal and conditional r-squared) 
r.squaredGLMM(mod3)

# Look at VCVs
vcov(mod3)
VarCorr(mod3)

# Calculate fit statistics (MSE and RMSE)
mse <- mean((obs_vs_pred$observed - obs_vs_pred$predicted)^2) # mean squared error
rmse <- sqrt(mse) # root mean squared error

# -----------------------------------------
# Other things that I tried - ignore here and below
# -----------------------------------------

# Plot autocorrelation in residuals: 
# Use the autocorrelation function (ACF) plot to check if the AR1 terms have sufficiently accounted for autocorrelation (sig peaks in the ACF plot suggest remaining temporal or spatial autocorrelation)
acf(residuals)
acf(residuals(mod3, type = "pearson"))

# Could simplify the AR1 structure: include AR1 for time only and treat depth as a random intercept or slope

# if I take out the ar1 for depth then the dharma looks good...interesting
data_final$day <- as.factor(data_final$day)
data_final$specified.depth <- as.numeric(data_final$specified.depth)
mod4 <- glmmTMB(chlorophyll.a ~ oxygen.concentration.z +
                  fish_treatment + timing +
                  fish_treatment*timing +
                  (1|specified.depth) +
                  # specified.depth
                  TP_ugL.z + TN_ugL.z + total_zoop_biomass_ugL.z +
                  ar1(day + 0 | enclosure), # spatial AR1 within each enclosure
                REML = F, # F when you're doing fixed effects model selection
                data = data_final,
                family = Gamma(link = "log"))
summary(mod4)
plot(data_final$specified.depth, data_final$chlorophyll.a)

# Resources
# https://cran.r-project.org/web/packages/glmmTMB/vignettes/covstruct.html#the-ar1-covariance-structure
# https://cran.r-project.org/web/packages/DHARMa/DHARMa.pdf
# https://ladal.edu.au/regression.html#Mixed-Effects_Ordinal_Regression
# https://ourcodingclub.github.io/tutorials/mixed-models/
# https://github.com/florianhartig/DHARMa/issues/368
# https://github.com/florianhartig/DHARMa/issues/364
# https://bbolker.github.io/mixedmodels-misc/notes/corr_braindump.html
# https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html
# https://bbolker.github.io/mixedmodels-misc/ecostats_chap.html