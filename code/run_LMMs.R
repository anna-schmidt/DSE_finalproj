# Running LMMs on SPRING-FISH profiler data and covariates
# AGS
# December 7, 2024

library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(sjPlot)

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
                  TP_ugL.z + TN_ugL.z + total_zoop_biomass_ugL.z,
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

# In order to fit the model with glmmTMB we must first specify a time variable as a factor. The factor levels correspond to unit spaced time points. It is a common mistake to forget some factor levels due to missing data or to order the levels incorrectly. We therefore recommend to construct factors with explicit levels, using the levels argument to the factor function

# first, make day and specified.depth factors
data_final$day <- as.factor(data_final$day)
data_final$specified.depth <- as.factor(data_final$specified.depth)

mod3 <- glmmTMB(chlorophyll.a ~ oxygen.concentration.z +
                  fish_treatment + timing +
                  fish_treatment*timing +
                  TP_ugL.z + TN_ugL.z + total_zoop_biomass_ugL.z +
                  ar1(specified.depth + 0 | enclosure) + # spatial AR1 within each enclosure
                  ar1(day + 0 | enclosure),          # temporal AR1 within each enclosure
                REML = F, # F when you're doing fixed effects model selection
                data = data_final,
                family = Gamma(link = "log"))
summary(mod3)

# Use AIC to compare which model is a better fit
AIC(mod1,  mod3)

# Look at VCVs
vcov(mod3)
VarCorr(mod3)

# DHARMa for mod3
simOut <- simulateResiduals(mod3, n = 250)
plot(simOut)
# not very good
# if I take out the ar1 for depth then the dharma looks good...interesting

# Plot predicted vs. observed values
predicted_values <- predict(mod3, type = "response")
observed_values <- mod3$frame$chlorophyll.a
results <- data.frame(
  Observed = observed_values,
  Predicted = predicted_values)
ggplot(results, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(x = "Observed Values", y = "Predicted Values") +
  theme_minimal() +
  ggtitle("Predicted vs. Observed Values")

# Calclate fit statistics
# Mean squared error
mse <- mean((results$Observed - results$Predicted)^2)
# Root mean squared error
rmse <- sqrt(mse)
# Print metrics
cat("MSE:", mse, "\nRMSE:", rmse)

#----------------------------------------------------
# Trying out ideas from online

# Check Residual Patterns by Time and Depth: Plot the residuals against time and depth indices. Ideally, residuals should show no clear pattern if the AR1 structure has effectively captured the temporal and spatial autocorrelation.
data_final2 <- data_final %>% subset(chlorophyll.a > 0)
residuals <- resid(mod3)
# NAs need to be removed first
plot(data_final2$day, residuals, main = "Residuals vs Time Index")
plot(data_final2$specified.depth, residuals, main = "Residuals vs Depth Index")

# Autocorrelation in Residuals: Use the autocorrelation function (ACF) plot to check if the AR1 terms have sufficiently accounted for autocorrelation. Significant peaks in the ACF plot suggest remaining temporal or spatial autocorrelation.
acf(residuals, main = "ACF of Residuals")

# Goodness-of-Fit Metrics: Marginal and Conditional R²: Use the MuMIn package to calculate marginal R-squared (variance explained by fixed effects) and conditional R-squared (variance explained by both fixed and random effects).
library(MuMIn)
r.squaredGLMM(mod3)

# Compare predicted vs observed values (A good fit will show points closely clustered around the 1:1 line)
predicted <- predict(mod3, re.form = NULL)  # Include random effects
data_final2 <- data_final %>% subset(chlorophyll.a > 0)
plot(data_final2$chlorophyll.a, predicted, main = "Observed vs Predicted", 
     xlab = "Observed", ylab = "Predicted")
abline(0, 1, col = "red")

# Evaluate Random Effect Variance: Examine the estimated variances for AR1 random effects for time and depth. Large variances suggest these structures are capturing meaningful patterns, while near-zero variances suggest overfitting or unnecessary complexity.

#Consider whether the Gamma family with a log link is the most appropriate for your data. If your residual diagnostics continue to look poor:
# Check for overdispersion (e.g., dispersion parameter much larger than 1).
# Try alternative distributions (e.g., negative binomial for overdispersion or log-normal for skewed continuous data).
residual_deviance <- sum(resid(mod3, type = "pearson")^2)
residual_df <- df.residual(mod3)

# Calculate overdispersion ratio
overdispersion_ratio <- residual_deviance / residual_df
print(overdispersion_ratio)

# If diagnostics consistently suggest poor fit, consider: Simplifying the AR1 structure: For instance, include AR1 for time only and treat depth as a random intercept or slope, or vice versa.


#----------------------------------------------------
# Below here is just trying out plotting model output

plot_model(mod3, type = 'pred', terms = 'oxygen.concentration.z')

library(ggeffects)
plot(ggpredict(mod3, terms = "oxygen.concentration.z"))
plot(ggpredict(mod3, terms = "TN_ugL.z"))


## population-level prediction
nd_pop <- data.frame(Days=unique(data_final$day),
                     enclosure=NA)
predict(mod3, newdata = nd_pop)



# if specified.depth is a factor, how does it know which depths are next to each other? and same with week?

# Notes:
# The + 0 ensures that no intercept is included in the AR1 term (this is required syntax for the AR1 structure in glmmTMB)
#  Specifying two independent AR1 correlation structures, one for week (temporal) and one for specified.depth (spatial)
# The AR1 correlation structure in glmmTMB requires a grouping variable to define independent clusters where autocorrelation is assumed
# 1 as the grouping variable specifies that AR1 correlations are applied globally (i.e., across the whole dataset)
# The ar1 function in glmmTMB uses the levels of the factor to determine the temporal order. This avoids errors from non-sequential or improperly sorted numeric values
# When using an AR1 term in a glmmTMB model, you are assuming that the underlying correlation structure is not of primary interest, but rather a nuisance factor that needs to be accounted for to accurately estimate the effects of your fixed variables.

# DHARMa generates simulated residuals under the assumption of independence. If the AR1 terms successfully capture the autocorrelation in your data, DHARMa might still flag residuals as problematic because its default diagnostics do not fully account for the modeled correlation structure.







# Plot Autocorrelation in Residuals
acf(residuals(mod3, type = "pearson"))

time_variable <- as.numeric(as.character(data_final$day))
# Run DHARMa with Autocorrelation Testing
testTemporalAutocorrelation(simulationOutput = simulateResiduals(mod3),
                            time = time_variable)

############### Trying another model...just the end of the experiment
# So no time component
data_final_lastday <- data_final %>% subset(day == 39)

mod4 <- glmmTMB(chlorophyll.a ~ oxygen.concentration.z +
                  fish_treatment +
                  TP_ugL.z + TN_ugL.z + total_zoop_biomass_ugL.z +
                  ar1(specified.depth + 0 | enclosure), # spatial AR1 within each enclosure
                REML = F, # F when you're doing fixed effects model selection
                data = data_final_lastday,
                family = Gamma(link = "log"))
summary(mod4)
# pretty bad
