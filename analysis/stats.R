#statistics

library(tidyverse)
library(lme4)
library(car)
library(MuMIn)
#library(whatever it is for AICc)


setwd("~/GitHub/EB_Ch1/analysis")


df <- read_csv("eb_results_summary_wide.csv")

pops <- readRDS("C:/Users/jmsmi/OneDrive/Documents/GitHub/EB_Ch1/data/pops.rds")

by <- join_by(species==spp, site_orig==site, sex)
df <- left_join(df, pops, by)

#data preparation
# 1. Average across sexes for main analyses
df_avg <- df %>%
  group_by(species, year, site_orig, site_clim, year_period) %>%
  summarize(
    energy_gain = mean(partial_shade_low_veg/mass, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Keep sex-specific data for robustness checks
df_male <- df %>% filter(sex == "M") %>% mutate(energy_gain=partial_shade_low_veg/mass)
df_female <- df %>% filter(sex == "F") %>% mutate(energy_gain=partial_shade_low_veg/mass)


#analysis 1: contemporary reciprocal transplant (fig 3)
# Filter to contemporary only
df_contemporary <- df_avg %>% 
  filter(year_period == "contemporary")

# Model 1a: Full reciprocal transplant model
mod_rt <- lmer(energy_gain ~ 
                 site_orig * site_clim +  # Main effects + interaction
                 (1|year),                # Random year effect
               data = df_contemporary %>% filter(species == "MB"))

# Model 1a: Full reciprocal transplant model
mod_rt_alt <- lmer(energy_gain ~ 
                 site_orig + site_clim +  # Main effects 
                 (1|year),                # Random year effect
               data = df_contemporary %>% filter(species == "MB"))

#mod_rt (more complicated model does better)
#AIC(mod_rt, mod_rt_alt)
AICc(mod_rt, mod_rt_alt)
print("select mod_rt, but let's look at both")
summary(mod_rt)
summary(mod_rt_alt)

# Model 1b: Same for MS
mod_rt_MS <- lmer(energy_gain ~ 
                    site_orig * site_clim + 
                    (1|year),
                  data = df_contemporary %>% filter(species == "MS"))

mod_rt_MS_alt <- lmer(energy_gain ~ 
                    site_orig + site_clim + 
                    (1|year),
                  data = df_contemporary %>% filter(species == "MS"))

#mod_rt (more complicated model does better)
AIC(mod_rt_MS, mod_rt_MS_alt)


#specific hypotheses to test
library(emmeans)
library(tidyverse)

# 1. Test for local adaptation (home vs away)
df_contemporary <- df_contemporary %>%
  mutate(is_home = site_orig == site_clim)

mod_local <- lmer(energy_gain ~ 
                    site_clim + site_orig + is_home +
                    (1|year),
                  data = df_contemporary %>% filter(species == "MB"))

# Extract coefficient for is_home
summary(mod_local)

# 2. Decompose variance: What matters more, climate or physiology?
# Compare nested models
mod_climate_only <- lmer(energy_gain ~ site_clim + (1|year), 
                         data = df_contemporary %>% filter(species == "MB"))
mod_physiol_only <- lmer(energy_gain ~ site_orig + (1|year), 
                         data = df_contemporary %>% filter(species == "MB"))

# Compare model fit
AIC(mod_climate_only, mod_physiol_only, mod_rt)

# Variance partitioning
library(performance)
r2(mod_climate_only)  # How much variance does climate explain?
r2(mod_physiol_only)  # How much variance does physiology explain?
r2(mod_rt)           # How much does the full model explain?





#THIS
#a more sophisticated approach
library(glmm.hp)

# For MB contemporary
mb_data <- df_contemporary %>% filter(species == "MB")
# mod_rt <- lmer(log_energy_gain ~ site_orig + site_clim + (1|year), data = mb_data)
# 
# # Partition variance
# hp_result <- glmm.hp(mod_rt)
# print(hp_result)
# plot(hp_result)



# For MB contemporary
ms_data <- df_contemporary %>% filter(species == "MS")
# mod_rt <- lmer(log_energy_gain ~ site_orig + site_clim  + (1|year), data = ms_data)
# 
# # Partition variance
# hp_result <- glmm.hp(mod_rt)
# print(hp_result)
# plot(hp_result)

#including the interaction term might actually be bad? (for hierarchical partitioning)

#other ideas

# Test if interaction is significant
mod_with_int <- lmer(energy_gain ~ site_orig * site_clim + (1|year), data = mb_data)
mod_without_int <- lmer(energy_gain ~ site_orig + site_clim + (1|year), data = mb_data)

anova(mod_with_int, mod_without_int)  # Is interaction significant?

AIC(mod_with_int, mod_without_int)

#MB
# # Log-transform for proportional (multiplicative) effects
# df_contemporary <- df_contemporary %>%
#   mutate(log_energy_gain = log(energy_gain))
# 
# mb_data <- df_contemporary %>% filter(species == "MB")
# 
# mod_no_int_log <- lmer(log_energy_gain ~ site_orig + site_clim + (1|year), 
#                    data = mb_data, REML=FALSE)
# 
# mod_int_log <- lmer(log_energy_gain ~ site_orig * site_clim + (1|year), 
#                        data = mb_data, REML=FALSE)
# 
# 
# anova(mod_int_log, mod_no_int_log)  # Is interaction significant?
# 
# AICc(mod_int_log, mod_no_int_log)
# 
# #refitting with REML (since we've already done the interaction vs. no interaction comparison)
# mod_no_int_log <- lmer(log_energy_gain ~ site_orig + site_clim + (1|year), 
#                        data = mb_data)
# 
# mod_int_log <- lmer(log_energy_gain ~ site_orig * site_clim + (1|year), 
#                     data = mb_data)
# 
# summary(mod_no_int_log)
# 
# summary(mod_int_log)
# 
# # Partition variance
# hp_result <- glmm.hp(mod_no_int_log)
# print(hp_result)

#and no log
mod_no_int <- lmer(energy_gain ~ site_orig + site_clim + (1|year), 
                       data = mb_data, REML=FALSE)

mod_int <- lmer(energy_gain ~ site_orig * site_clim + (1|year), 
                    data = mb_data, REML=FALSE)


anova(mod_int, mod_no_int)  # Is interaction significant?

AICc(mod_int, mod_no_int)

#refitting with REML (since we've already done the interaction vs. no interaction comparison)
mod_no_int <- lmer(energy_gain ~ site_orig + site_clim + (1|year), 
                       data = mb_data)

mod_int <- lmer(energy_gain ~ site_orig * site_clim + (1|year), 
                    data = mb_data)

summary(mod_no_int)

summary(mod_int)

Anova(mod_int, type="III")

# Partition variance
hp_result <- glmm.hp(mod_no_int)
print(hp_result)



#MS
# # Log-transform for proportional (multiplicative) effects
# df_contemporary <- df_contemporary %>%
#   mutate(log_energy_gain = log(energy_gain))
# 
# ms_data <- df_contemporary %>% filter(species == "MS")
# 
# mod_no_int_log <- lmer(log_energy_gain ~ site_orig + site_clim + (1|year), 
#                        data = ms_data, REML=FALSE)
# 
# mod_int_log <- lmer(log_energy_gain ~ site_orig * site_clim + (1|year), 
#                     data = ms_data, REML=FALSE)
# 
# 
# anova(mod_int_log, mod_no_int_log)  # Is interaction significant?
# 
# AICc(mod_int_log, mod_no_int_log)
# 
# #refitting with REML (since we've already done the interaction vs. no interaction comparison)
# mod_no_int_log <- lmer(log_energy_gain ~ site_orig + site_clim + (1|year), 
#                        data = ms_data)
# 
# mod_int_log <- lmer(log_energy_gain ~ site_orig * site_clim + (1|year), 
#                     data = ms_data)
# 
# summary(mod_no_int_log)
# 
# summary(mod_int_log)
# 
# # Partition variance
# hp_result <- glmm.hp(mod_no_int_log)
# print(hp_result)

#and no log
mod_no_int <- lmer(energy_gain ~ site_orig + site_clim + (1|year), 
                   data = ms_data, REML=FALSE)

mod_int <- lmer(energy_gain ~ site_orig * site_clim + (1|year), 
                data = ms_data, REML=FALSE)


anova(mod_int, mod_no_int)  # Is interaction significant?

AICc(mod_int, mod_no_int) #76

#refitting with REML (since we've already done the interaction vs. no interaction comparison)
mod_no_int <- lmer(energy_gain ~ site_orig + site_clim + (1|year), 
                   data = ms_data)

mod_int <- lmer(energy_gain ~ site_orig * site_clim + (1|year), 
                data = ms_data)

summary(mod_no_int)

summary(mod_int)

Anova(mod_int, type="III")

# Partition variance
hp_result <- glmm.hp(mod_no_int)
print(hp_result)


#moving this up for ease, but first the shorter version

# #MB log
df_temporal <- df_avg %>%
  filter(site_orig == site_clim) %>%
  mutate(log_energy_gain = log(energy_gain))
# 
# mod_temp_log_MB <- lmer(log_energy_gain ~ 
#                       year_period * site_orig +  # Period × site interaction
#                       (1|year),                   # Random year effect
#                     data = df_temporal %>% filter(species == "MB"))
# summary(mod_temp_log_MB)
# 
# #specific hypotheses to test
# # 1. Overall temporal change (across all sites)
# emm_period <- emmeans(mod_temp_log_MB, ~ year_period)
# pairs(emm_period)
# 
# 
# # 2. Site-specific temporal changes
# emm_site_period <- emmeans(mod_temp_log_MB, ~ year_period | site_orig)
# contrasts_temporal <- pairs(emm_site_period, adjust = "none")
# summary(contrasts_temporal)


#MB no log
mod_temp_MB <- lmer(energy_gain ~ 
                      year_period * site_orig +  # Period × site interaction
                      (1|year),                   # Random year effect
                    data = df_temporal %>% filter(species == "MB"))
summary(mod_temp_MB)

#specific hypotheses to test
# 1. Overall temporal change (across all sites)
emm_period <- emmeans(mod_temp_MB, ~ year_period)
pairs(emm_period)


# 2. Site-specific temporal changes
emm_site_period <- emmeans(mod_temp_MB, ~ year_period | site_orig)
contrasts_temporal <- pairs(emm_site_period, adjust = "none")
summary(contrasts_temporal)

Anova(mod_temp_MB, type="III")



#MS no log
mod_temp_MS <- lmer(energy_gain ~ 
                      year_period * site_orig +  # Period × site interaction
                      (1|year),                   # Random year effect
                    data = df_temporal %>% filter(species == "MS"))
summary(mod_temp_MS)

Anova(mod_temp_MS, type="III")

#specific hypotheses to test
# 1. Overall temporal change (across all sites)
emm_period <- emmeans(mod_temp_MS, ~ year_period)
pairs(emm_period)


# 2. Site-specific temporal changes
emm_site_period <- emmeans(mod_temp_MS, ~ year_period | site_orig)
contrasts_temporal <- pairs(emm_site_period, adjust = "none")
summary(contrasts_temporal)

#End of THIS




#the longer version (AICc comparisons to select interaction)

#Analysis 2: time periods (fig 4)
# Use all data, but only "home" sites (no transplant)
df_temporal <- df_avg %>%
  filter(site_orig == site_clim) %>%
  mutate(log_energy_gain = log(energy_gain))

# Model 2: Temporal change across sites (MB)
mod_temp_MB <- lmer(energy_gain ~ 
                   year_period * site_orig +  # Period × site interaction
                   (1|year),                   # Random year effect
                   REML=FALSE,
                 data = df_temporal %>% filter(species == "MB"))

mod_temp_no_int_MB <- lmer(energy_gain ~ 
                      year_period + site_orig +  # Period × site interaction
                      (1|year),                   # Random year effect
                      REML=FALSE,
                    data = df_temporal %>% filter(species == "MB"))

AICc(mod_temp_MB, mod_temp_no_int_MB)

#after comparison, change back to REML:
mod_temp_MB <- lmer(energy_gain ~ 
                      year_period * site_orig +  # Period × site interaction
                      (1|year),                   # Random year effect
                    data = df_temporal %>% filter(species == "MB"))


summary(mod_temp_MB)
Anova(mod_temp_MB, type="III")

#specific hypotheses to test
# 1. Overall temporal change (across all sites)
emm_period <- emmeans(mod_temp_MB, ~ year_period)
pairs(emm_period)


# 2. Site-specific temporal changes
emm_site_period <- emmeans(mod_temp_MB, ~ year_period | site_orig)
contrasts_temporal <- pairs(emm_site_period, adjust = "none")
summary(contrasts_temporal)



#now log
log_mod_temp_MB <- lmer(log_energy_gain ~ 
                      year_period * site_orig +  # Period × site interaction
                      (1|year),                   # Random year effect
                      REML=FALSE,
                    data = df_temporal %>% filter(species == "MB"))

log_mod_temp_no_int_MB <- lmer(log_energy_gain ~ 
                             year_period + site_orig +  # Period × site interaction
                             (1|year),                   # Random year effect
                             REML=FALSE,
                           data = df_temporal %>% filter(species == "MB"))

AICc(log_mod_temp_MB, log_mod_temp_no_int_MB)

#back to REML
log_mod_temp_MB <- lmer(log_energy_gain ~ 
                          year_period * site_orig +  # Period × site interaction
                          (1|year),                   # Random year effect
                        data = df_temporal %>% filter(species == "MB"))

log_mod_temp_no_int_MB <- lmer(log_energy_gain ~ 
                                 year_period + site_orig +  # Period × site interaction
                                 (1|year),                   # Random year effect
                               data = df_temporal %>% filter(species == "MB"))

summary(log_mod_temp_MB)
summary(log_mod_temp_no_int_MB)


#specific hypotheses to test
# 1. Overall temporal change (across all sites)
emm_period <- emmeans(log_mod_temp_no_int_MB, ~ year_period)
pairs(emm_period)


# 2. Site-specific temporal changes
emm_site_period <- emmeans(log_mod_temp_MB, ~ year_period | site_orig)
contrasts_temporal <- pairs(emm_site_period, adjust = "none")
summary(contrasts_temporal)



#MS
#Analysis 2: time periods (fig 4)
# Use all data, but only "home" sites (no transplant)
df_temporal <- df_avg %>%
  filter(site_orig == site_clim) %>%
  mutate(log_energy_gain = log(energy_gain))

# Model 2: Temporal change across sites (MS)
mod_temp_MS <- lmer(energy_gain ~ 
                      year_period * site_orig +  # Period × site interaction
                      (1|year),                   # Random year effect
                    REML=FALSE,
                    data = df_temporal %>% filter(species == "MS"))

mod_temp_no_int_MS <- lmer(energy_gain ~ 
                             year_period + site_orig +  # Period × site interaction
                             (1|year),                   # Random year effect
                           REML=FALSE,
                           data = df_temporal %>% filter(species == "MS"))

AICc(mod_temp_MS, mod_temp_no_int_MS)

#after comparison, change back to REML
mod_temp_MS <- lmer(energy_gain ~ 
                      year_period * site_orig +  # Period × site interaction
                      (1|year),                   # Random year effect
                    data = df_temporal %>% filter(species == "MS"))

mod_temp_no_int_MS <- lmer(energy_gain ~ 
                             year_period + site_orig +  # Period × site interaction
                             (1|year),                   # Random year effect
                           data = df_temporal %>% filter(species == "MS"))

summary(mod_temp_MS)
summary(mod_temp_no_int_MS)

Anova(mod_temp_MS, type="III")

#specific hypotheses to test
# 1. Overall temporal change (across all sites)
emm_period <- emmeans(mod_temp, ~ year_period)
pairs(emm_period)


# 2. Site-specific temporal changes
emm_site_period <- emmeans(mod_temp, ~ year_period | site_orig)
contrasts_temporal <- pairs(emm_site_period, adjust = "none")
summary(contrasts_temporal)



#now log
log_mod_temp_MS <- lmer(log_energy_gain ~ 
                          year_period * site_orig +  # Period × site interaction
                          (1|year),                   # Random year effect
                        REML=FALSE,
                        data = df_temporal %>% filter(species == "MS"))

log_mod_temp_no_int_MS <- lmer(log_energy_gain ~ 
                                 year_period + site_orig +  # Period × site interaction
                                 (1|year),                   # Random year effect
                               REML=FALSE,
                               data = df_temporal %>% filter(species == "MS"))

AICc(log_mod_temp_MS, log_mod_temp_no_int_MS)

summary(log_mod_temp_MS)
summary(log_mod_temp_no_int_MS)


#specific hypotheses to test
# 1. Overall temporal change (across all sites)
emm_period <- emmeans(log_mod_temp_no_int_MS, ~ year_period)
pairs(emm_period)


# 2. Site-specific temporal changes
emm_site_period <- emmeans(log_mod_temp_MS, ~ year_period | site_orig)
contrasts_temporal <- pairs(emm_site_period, adjust = "none")
summary(contrasts_temporal)






#THE SAME BUT FOR MS
mod_local <- lmer(energy_gain ~ 
                    site_clim + site_orig + is_home +
                    (1|year),
                  data = df_contemporary %>% filter(species == "MS"))

# Extract coefficient for is_home
summary(mod_local)

# 2. Decompose variance: What matters more, climate or physiology?
# Compare nested models
mod_climate_only <- lmer(energy_gain ~ site_clim + (1|year), 
                         data = df_contemporary %>% filter(species == "MS"))
mod_physiol_only <- lmer(energy_gain ~ site_orig + (1|year), 
                         data = df_contemporary %>% filter(species == "MS"))

# Compare model fit
AIC(mod_climate_only, mod_physiol_only, mod_rt)

# Variance partitioning
library(performance)
r2(mod_climate_only)  # How much variance does climate explain?
r2(mod_physiol_only)  # How much variance does physiology explain?
r2(mod_rt)           # How much does the full model explain?



# 3. Pairwise comparisons at each climate site
emm_rt <- emmeans(mod_rt, ~ site_orig | site_clim)
pairs(emm_rt, adjust = "tukey")

# 4. Test for elevational patterns
# Add elevation as continuous variable
elev_data <- data.frame(
  site = c("Eldo", "A1", "B1", "C1", "D1"),
  elevation = c(1740, 2195, 2591, 3014, 3515)  # Fill in D1 elevation
)

df_contemporary <- df_contemporary %>%
  left_join(elev_data, by = c("site_orig" = "site")) %>%
  rename(elev_orig = elevation) %>%
  left_join(elev_data, by = c("site_clim" = "site")) %>%
  rename(elev_clim = elevation)

# Test if energy gain decreases with elevation
mod_elev <- lmer(energy_gain ~ 
                   scale(elev_clim) * scale(elev_orig) +
                   (1|year),
                 data = df_contemporary %>% filter(species == "MB"))


#Analysis 2: time periods (fig 4)
# Use all data, but only "home" sites (no transplant)
df_temporal <- df_avg %>%
  filter(site_orig == site_clim)

# Model 2a: Temporal change across sites for MB
mod_temp <- lmer(energy_gain ~ 
                   year_period * site_orig +  # Period × site interaction
                   (1|year),                   # Random year effect
                 data = df_temporal %>% filter(species == "MB"))

# Model 2b: Same for MS
mod_temp_MS <- lmer(energy_gain ~ 
                      year_period * site_orig +
                      (1|year),
                    data = df_temporal %>% filter(species == "MS"))

#specific hypotheses to test
# 1. Overall temporal change (across all sites)
emm_period <- emmeans(mod_temp, ~ year_period)
pairs(emm_period)

# Calculate percent change
period_means <- summary(emm_period)
pct_change <- (period_means$emmean[2] - period_means$emmean[1]) / 
  period_means$emmean[1] * 100

# 2. Site-specific temporal changes
emm_site_period <- emmeans(mod_temp, ~ year_period | site_orig)
contrasts_temporal <- pairs(emm_site_period, adjust = "none")
summary(contrasts_temporal)

# 3. Does temporal change vary with elevation?
mod_temp_elev <- lmer(energy_gain ~ 
                        year_period * scale(elev_orig) +
                        (1|site_orig) +
                        (1|year),
                      data = df_temporal %>% 
                        left_join(elev_data, by = c("site_orig" = "site")) %>%
                        rename(elev_orig = elevation) %>%
                        filter(species == "MB"))

# Test the interaction - does climate change impact vary by elevation?
anova(mod_temp_elev)

#sex robustness check
# Run the same models for each sex separately
# For reciprocal transplant
mod_rt_male <- lmer(energy_gain ~ site_orig * site_clim + (1|year),
                    data = df_male %>% 
                      filter(year_period == "contemporary", species == "MB"))

mod_rt_female <- lmer(energy_gain ~ site_orig * site_clim + (1|year),
                      data = df_female %>% 
                        filter(year_period == "contemporary", species == "MB"))

# Compare coefficients
coef_comparison <- data.frame(
  effect = rownames(summary(mod_rt)$coefficients),
  avg_model = summary(mod_rt)$coefficients[,"Estimate"],
  male_model = summary(mod_rt_male)$coefficients[,"Estimate"],
  female_model = summary(mod_rt_female)$coefficients[,"Estimate"]
)

# If patterns are similar, report: "Results were qualitatively similar when 
# analyzed separately by sex (see supplementary materials)"



#summary statistics table

# Create a results summary table
results_summary <- df_avg %>%
  group_by(species, year_period, site_orig, site_clim) %>%
  summarize(
    mean_energy = mean(energy_gain),
    se_energy = sd(energy_gain)/sqrt(n()),
    n_years = n(),
    .groups = "drop"
  )

# For the paper, create a clean table of key contrasts
key_results <- data.frame(
  Analysis = c("Local adaptation (MB)", 
               "Local adaptation (MS)",
               "Climate effect (MB)",
               "Physiology effect (MB)",
               "Temporal change - overall (MB)",
               "Temporal change - low elev",
               "Temporal change - high elev"),
  Estimate = NA,  # Fill from models
  SE = NA,
  P_value = NA
)






#trying to address the R2 thing (for MS)
# Check for singularity
isSingular(mod_climate_only)  # Likely TRUE
isSingular(mod_physiol_only)  # Likely TRUE
isSingular(mod_rt_MS)         # Likely FALSE

# Look at the variance components
VarCorr(mod_climate_only)
VarCorr(mod_physiol_only)
VarCorr(mod_rt_MS)

# Solutions:

# Option 1: If year variance is truly negligible, use fixed effects model
mod_climate_only_fixed <- lm(energy_gain ~ site_clim, 
                             data = df_contemporary %>% filter(species == "MS"))
mod_physiol_only_fixed <- lm(energy_gain ~ site_orig, 
                             data = df_contemporary %>% filter(species == "MS"))

# Compare R-squared (now it's regular R2)
summary(mod_climate_only_fixed)$r.squared  # Climate explains X%
summary(mod_physiol_only_fixed)$r.squared  # Physiology explains X%

# Option 2: Force the calculation with lower tolerance
r2(mod_climate_only, tolerance = 1e-8)

# Option 3: Use a different metric - compare AIC weights
library(MuMIn)
AICc_climate <- AICc(mod_climate_only)
AICc_physiol <- AICc(mod_physiol_only)
AICc_full <- AICc(mod_rt_MS)

# Calculate Akaike weights
models <- list(climate = mod_climate_only, 
               physiol = mod_physiol_only, 
               full = mod_rt_MS)
aictab <- AICc(climate = mod_climate_only, 
               physiol = mod_physiol_only, 
               full = mod_rt_MS)