# Running LMMs on SPRING-FISH profiler data and covariates
# AGS
# December 7, 2024

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(sjPlot)

# ------------------ Running models ------------------

# Model without any random effects
mod1 <- glmmTMB(chlorophyll.a ~ oxygen.concentration.z +
                     fish_treatment + timing +
                     fish_treatment*timing +
                     TP_ugL.z + TN_ugL.z + total_zoop_biomass_ugL.z,
                   data = data_final,
                family = Gamma(link = "log"))
summary(mod1)

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

# Use AIC to compare which model is a better fit
AIC(mod1, mod2)
# mod2 is better

# DHARMa for mod2
simOut <- simulateResiduals(mod1, n = 250)
plot(simOut)
# okay

# Mixed effects model with enclosure as a random effect and ar1() terms for both time and space

# In order to fit the model with glmmTMB we must first specify a time variable as a factor. The factor levels correspond to unit spaced time points. It is a common mistake to forget some factor levels due to missing data or to order the levels incorrectly.

# first, make day and specified.depth factors
data_final$day <- as.factor(data_final$day)
data_final$specified.depth <- as.factor(data_final$specified.depth)

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

# Use AIC to compare which model is a better fit
AIC(mod1, mod2, mod3)

# ------------------ Model diagnostics ------------------

# Look at VCVs
vcov(mod3)
VarCorr(mod3)

# DHARMa for mod3
simOut <- simulateResiduals(mod3, n = 250)
plot(simOut)
# not very good

# Calculate fit statistics
mse <- mean((obs_vs_pred$observed - obs_vs_pred$predicted)^2) # mean squared error
rmse <- sqrt(mse) # root mean squared error

# Plot observed vs. predicted values for chlorophyll-a
# day and specified.depth need to be characters for this to work (see https://github.com/glmmTMB/glmmTMB/issues/602)
data_final$day <- as.character(data_final$day)
data_final$specified.depth <- as.character(data_final$specified.depth)
data_final$chlorophyll.a_predicted <- predict(mod3, newdata = data_final, type = "response")
data_final$day <- as.numeric(data_final$day)
ggplot(data_final, aes(x = chlorophyll.a, y = chlorophyll.a_predicted, 
               color = day)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  labs(x = "Observed chlorophyll-a", y = "Predicted chlorophyll-a") +
  theme_bw()
# color by enclosure, day, and depth to show the residuals are random

# Plot the residuals by day and specified.depth and look for patterns
#  Ideally, residuals should show no clear pattern if the AR1 structure has effectively captured the temporal and spatial autocorrelation.
# NAs need to be removed first
data_final_noNAs <- data_final %>% subset(chlorophyll.a > 0)
residuals <- resid(mod3)

plot(data_final_noNAs$day, residuals, main = "Residuals vs Time Index")
plot(data_final_noNAs$specified.depth, residuals, main = "Residuals vs Depth Index")

# Plot the model results
# oxygen concentration vs. chlorophyll-a
plot_model(mod3, type = "pred", terms = "oxygen.concentration.z",
           axis.title = c("Oxygen concentration", "Chlorophyll-a"))




# Plot autocorrelation in Residuals: Use the autocorrelation function (ACF) plot to check if the AR1 terms have sufficiently accounted for autocorrelation. Significant peaks in the ACF plot suggest remaining temporal or spatial autocorrelation.
acf(residuals)
acf(residuals(mod3, type = "pearson"))

# Goodness-of-Fit Metrics: Marginal and Conditional R²: Use the MuMIn package to calculate marginal R-squared (variance explained by fixed effects) and conditional R-squared (variance explained by both fixed and random effects).
library(MuMIn)
r.squaredGLMM(mod3)

# Evaluate Random Effect Variance: Examine the estimated variances for AR1 random effects for time and depth. Large variances suggest these structures are capturing meaningful patterns, while near-zero variances suggest overfitting or unnecessary complexity.

# Simplifying the AR1 structure: include AR1 for time only and treat depth as a random intercept or slope, or vice versa.

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


#----------------------------------------------------
# Below here is just trying out plotting model output

plot_model(mod3, type = 'pred', terms = 'oxygen.concentration.z')

library(ggeffects)
plot(ggpredict(mod3, terms = "oxygen.concentration.z"))
plot(ggpredict(mod3, terms = "TN_ugL.z"))


# if specified.depth is a factor, how does it know which depths are next to each other? and same with week?

# Notes:
# The + 0 ensures that no intercept is included in the AR1 term (this is required syntax for the AR1 structure in glmmTMB)
#  Specifying two independent AR1 correlation structures, one for week (temporal) and one for specified.depth (spatial)
# The AR1 correlation structure in glmmTMB requires a grouping variable to define independent clusters where autocorrelation is assumed
# 1 as the grouping variable specifies that AR1 correlations are applied globally (i.e., across the whole dataset)
# The ar1 function in glmmTMB uses the levels of the factor to determine the temporal order. This avoids errors from non-sequential or improperly sorted numeric values
# When using an AR1 term in a glmmTMB model, you are assuming that the underlying correlation structure is not of primary interest, but rather a nuisance factor that needs to be accounted for to accurately estimate the effects of your fixed variables.

# DHARMa generates simulated residuals under the assumption of independence. If the AR1 terms successfully capture the autocorrelation in your data, DHARMa might still flag residuals as problematic because its default diagnostics do not fully account for the modeled correlation structure.

# Trying altenative 

time_variable <- as.numeric(as.character(data_final$day))
# Run DHARMa with Autocorrelation Testing
testTemporalAutocorrelation(simulationOutput = simulateResiduals(mod3),
                            time = time_variable)

# Our conclusions from the DHARMa simulated residuals are as follows: The QQ plot shows a pretty linear trend of points that tend to fall close to the red diagonal line, though some points do deviate from the line near the middle. We detect some overall deviations from the expected distribution with the dispersion test, but the KS test and outlier test are not significant. The residuals plot shows that some of the residual vs. model predictions lines are curved (red), which refers to quantile deviations, and we can also see four red stars which represent simulation outliers (data points that are outside the range of simulated values), which should be carefully interpreted. So overall, we have some reason to be hesitant about our model because of the items mentioned above that are pointed out in red on the figure. However, residual patterns don't indicate that the model is necessarily unusable, as we discussed in class.
