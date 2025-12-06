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


df_contemporary <- df_avg %>%
  filter(year_period == "contemporary")

library(glmm.hp)

# For MB contemporary
mb_data <- df_contemporary %>% filter(species == "MB")




# For MS contemporary
ms_data <- df_contemporary %>% filter(species == "MS")


#MB

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












library(lme4)
library(lmerTest)
library(car)
library(emmeans)
library(kableExtra)
library(dplyr)

setwd("~/GitHub/EB_Ch1/analysis")

df <- read_csv("eb_results_summary_wide.csv")
pops <- readRDS("C:/Users/jmsmi/OneDrive/Documents/GitHub/EB_Ch1/data/pops.rds")

by <- join_by(species==spp, site_orig==site, sex)
df <- left_join(df, pops, by)

# Site elevation lookup
site_coords <- data.frame(
  name = c("Eldo", "A1", "B1", "C1", "D1"),
  elev = c(1740, 2195, 2591, 3048, 3515)
)

# Function to replace site names with elevations
replace_sites <- function(term_names) {
  term_names <- gsub("site_orig", "pop", term_names)
  term_names <- gsub("site_clim", "clim", term_names)
  term_names <- gsub("year_period", "time_period", term_names)
  for(i in 1:nrow(site_coords)) {
    term_names <- gsub(site_coords$name[i], 
                       paste0(site_coords$elev[i], "m"), 
                       term_names)
  }
  return(term_names)
}

# Function to order rows by elevation
order_by_elevation <- function(df, term_col = "Term") {
  # Extract elevation numbers from terms
  df$elev_order <- NA
  for(i in 1:nrow(site_coords)) {
    pattern <- paste0(site_coords$elev[i], "m")
    matches <- grepl(pattern, df[[term_col]], fixed = TRUE)
    df$elev_order[matches] <- site_coords$elev[i]
  }
  
  # Put rows without elevation at top, then order by elevation
  df <- df %>%
    arrange(is.na(elev_order), elev_order) %>%
    select(-elev_order)
  
  return(df)
}

# 1. Average across sexes for main analyses
df_avg <- df %>%
  group_by(species, year, site_orig, site_clim, year_period) %>%
  summarize(
    energy_gain = mean(partial_shade_low_veg/mass, na.rm = TRUE),
    .groups = "drop"
  )

# ========================================
# CONTEMPORARY TRANSPLANT EXPERIMENT
# ========================================

df_contemporary <- df_avg %>%
  filter(year_period == "contemporary")

# Fit models for both species
mb_data <- df_contemporary %>% filter(species == "MB")
mod_int_mb <- lmer(energy_gain ~ site_orig * site_clim + (1|year), 
                   data = mb_data)

ms_data <- df_contemporary %>% filter(species == "MS")
mod_int_ms <- lmer(energy_gain ~ site_orig * site_clim + (1|year), 
                   data = ms_data)

# ===== TABLE 1: ANOVA RESULTS (CONTEMPORARY) =====
anova_mb <- Anova(mod_int_mb, type="III")
anova_ms <- Anova(mod_int_ms, type="III")

anova_table <- data.frame(
  Term = replace_sites(rownames(anova_mb)),
  MB_Chisq = anova_mb$Chisq,
  MB_df = anova_mb$Df,
  MB_p = anova_mb$`Pr(>Chisq)`,
  MS_Chisq = anova_ms$Chisq,
  MS_df = anova_ms$Df,
  MS_p = anova_ms$`Pr(>Chisq)`
)

anova_table$MB_p <- ifelse(anova_table$MB_p < 0.001, "< 0.001", 
                           sprintf("%.3f", anova_table$MB_p))
anova_table$MS_p <- ifelse(anova_table$MS_p < 0.001, "< 0.001", 
                           sprintf("%.3f", anova_table$MS_p))

anova_table <- order_by_elevation(anova_table)

anova_kable <- anova_table %>%
  kbl(
    col.names = c("", "χ²", "df", "p", "χ²", "df", "p"),
    digits = 1,
    align = c("l", "r", "r", "r", "r", "r", "r"),
    escape = FALSE,
    row.names = FALSE
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    position = "left"
  ) %>%
  add_header_above(c(" " = 1, "M. boulderensis" = 3, "M. sanguinipes" = 3)) %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1, width = "3cm") %>%
  column_spec(2:7, width = "1.5cm")

print(anova_kable)

# ===== TABLE 2: FIXED EFFECTS (CONTEMPORARY) =====
coef_mb <- summary(mod_int_mb)$coefficients
coef_ms <- summary(mod_int_ms)$coefficients

# Get all unique terms from both models
all_terms <- unique(c(rownames(coef_mb), rownames(coef_ms)))

# Create data frame with NA for missing terms
coef_table <- data.frame(
  Term = replace_sites(all_terms),
  MB_Estimate = NA,
  MB_SE = NA,
  MB_t = NA,
  MB_p = NA,
  MS_Estimate = NA,
  MS_SE = NA,
  MS_t = NA,
  MS_p = NA
)

# Fill in MB data
mb_match <- match(all_terms, rownames(coef_mb))
for(i in 1:length(all_terms)) {
  if(!is.na(mb_match[i])) {
    coef_table$MB_Estimate[i] <- coef_mb[mb_match[i], "Estimate"]
    coef_table$MB_SE[i] <- coef_mb[mb_match[i], "Std. Error"]
    coef_table$MB_t[i] <- coef_mb[mb_match[i], "t value"]
    coef_table$MB_p[i] <- coef_mb[mb_match[i], "Pr(>|t|)"]
  }
}

# Fill in MS data
ms_match <- match(all_terms, rownames(coef_ms))
for(i in 1:length(all_terms)) {
  if(!is.na(ms_match[i])) {
    coef_table$MS_Estimate[i] <- coef_ms[ms_match[i], "Estimate"]
    coef_table$MS_SE[i] <- coef_ms[ms_match[i], "Std. Error"]
    coef_table$MS_t[i] <- coef_ms[ms_match[i], "t value"]
    coef_table$MS_p[i] <- coef_ms[ms_match[i], "Pr(>|t|)"]
  }
}

# Format p-values
coef_table$MB_p <- ifelse(is.na(coef_table$MB_p), "—",
                          ifelse(coef_table$MB_p < 0.001, "< 0.001", 
                                 sprintf("%.3f", coef_table$MB_p)))
coef_table$MS_p <- ifelse(is.na(coef_table$MS_p), "—",
                          ifelse(coef_table$MS_p < 0.001, "< 0.001", 
                                 sprintf("%.3f", coef_table$MS_p)))

coef_table <- order_by_elevation(coef_table)

coef_kable <- coef_table %>%
  kbl(
    col.names = c("", "Coefficient", "SE", "t", "p", 
                  "Coefficient", "SE", "t", "p"),
    digits = 3,
    align = c("l", "r", "r", "r", "r", "r", "r", "r", "r"),
    escape = FALSE,
    row.names = FALSE
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    position = "left"
  ) %>%
  add_header_above(c(" " = 1, "M. boulderensis" = 4, "M. sanguinipes" = 4)) %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1, width = "4cm") %>%
  column_spec(2:9, width = "1.5cm")

print(coef_kable)

# ========================================
# TEMPORAL COMPARISON
# ========================================

# Prepare temporal data
df_temporal <- df_avg %>%
  filter(site_orig == site_clim)

# Fit models for both species
mod_temp_MB <- lmer(energy_gain ~ year_period * site_orig + (1|year),
                    data = df_temporal %>% filter(species == "MB"))

mod_temp_MS <- lmer(energy_gain ~ year_period * site_orig + (1|year),
                    data = df_temporal %>% filter(species == "MS"))

# ===== TABLE 3: ANOVA RESULTS (TEMPORAL) =====
anova_mb <- Anova(mod_temp_MB, type="III")
anova_ms <- Anova(mod_temp_MS, type="III")

anova_table <- data.frame(
  Term = replace_sites(rownames(anova_mb)),
  MB_Chisq = anova_mb$Chisq,
  MB_df = anova_mb$Df,
  MB_p = anova_mb$`Pr(>Chisq)`,
  MS_Chisq = anova_ms$Chisq,
  MS_df = anova_ms$Df,
  MS_p = anova_ms$`Pr(>Chisq)`
)

anova_table$MB_p <- ifelse(anova_table$MB_p < 0.001, "< 0.001", 
                           sprintf("%.3f", anova_table$MB_p))
anova_table$MS_p <- ifelse(anova_table$MS_p < 0.001, "< 0.001", 
                           sprintf("%.3f", anova_table$MS_p))

anova_table <- order_by_elevation(anova_table)

anova_kable_temporal <- anova_table %>%
  kbl(
    col.names = c("", "χ²", "df", "p", "χ²", "df", "p"),
    digits = 1,
    align = c("l", "r", "r", "r", "r", "r", "r"),
    escape = FALSE,
    row.names = FALSE
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    position = "left"
  ) %>%
  add_header_above(c(" " = 1, "M. boulderensis" = 3, "M. sanguinipes" = 3)) %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1, width = "3cm") %>%
  column_spec(2:7, width = "1.5cm")

print(anova_kable_temporal)

# ===== TABLE 4: FIXED EFFECTS (TEMPORAL) =====
coef_mb <- summary(mod_temp_MB)$coefficients
coef_ms <- summary(mod_temp_MS)$coefficients

# Get all unique terms from both models
all_terms <- unique(c(rownames(coef_mb), rownames(coef_ms)))

# Create data frame with NA for missing terms
coef_table <- data.frame(
  Term = replace_sites(all_terms),
  MB_Estimate = NA,
  MB_SE = NA,
  MB_t = NA,
  MB_p = NA,
  MS_Estimate = NA,
  MS_SE = NA,
  MS_t = NA,
  MS_p = NA
)

# Fill in MB data
mb_match <- match(all_terms, rownames(coef_mb))
for(i in 1:length(all_terms)) {
  if(!is.na(mb_match[i])) {
    coef_table$MB_Estimate[i] <- coef_mb[mb_match[i], "Estimate"]
    coef_table$MB_SE[i] <- coef_mb[mb_match[i], "Std. Error"]
    coef_table$MB_t[i] <- coef_mb[mb_match[i], "t value"]
    coef_table$MB_p[i] <- coef_mb[mb_match[i], "Pr(>|t|)"]
  }
}

# Fill in MS data
ms_match <- match(all_terms, rownames(coef_ms))
for(i in 1:length(all_terms)) {
  if(!is.na(ms_match[i])) {
    coef_table$MS_Estimate[i] <- coef_ms[ms_match[i], "Estimate"]
    coef_table$MS_SE[i] <- coef_ms[ms_match[i], "Std. Error"]
    coef_table$MS_t[i] <- coef_ms[ms_match[i], "t value"]
    coef_table$MS_p[i] <- coef_ms[ms_match[i], "Pr(>|t|)"]
  }
}

# Format p-values
coef_table$MB_p <- ifelse(is.na(coef_table$MB_p), "—",
                          ifelse(coef_table$MB_p < 0.001, "< 0.001", 
                                 sprintf("%.3f", coef_table$MB_p)))
coef_table$MS_p <- ifelse(is.na(coef_table$MS_p), "—",
                          ifelse(coef_table$MS_p < 0.001, "< 0.001", 
                                 sprintf("%.3f", coef_table$MS_p)))

coef_table <- order_by_elevation(coef_table)

coef_kable_temporal <- coef_table %>%
  kbl(
    col.names = c("", "Coefficient", "SE", "t", "p", 
                  "Coefficient", "SE", "t", "p"),
    digits = 3,
    align = c("l", "r", "r", "r", "r", "r", "r", "r", "r"),
    escape = FALSE,
    row.names = FALSE
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    position = "left"
  ) %>%
  add_header_above(c(" " = 1, "M. boulderensis" = 4, "M. sanguinipes" = 4)) %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1, width = "4cm") %>%
  column_spec(2:9, width = "1.5cm")

print(coef_kable_temporal)

# ===== TABLE 5: SITE-SPECIFIC CONTRASTS =====
# Get contrasts for MB
emm_site_period_mb <- emmeans(mod_temp_MB, ~ year_period | site_orig)
contrasts_mb <- pairs(emm_site_period_mb, adjust = "none")
contrast_summary_mb <- summary(contrasts_mb)

# Get contrasts for MS
emm_site_period_ms <- emmeans(mod_temp_MS, ~ year_period | site_orig)
contrasts_ms <- pairs(emm_site_period_ms, adjust = "none")
contrast_summary_ms <- summary(contrasts_ms)

# Get all unique sites from both species
all_sites <- unique(c(as.character(contrast_summary_mb$site_orig), 
                      as.character(contrast_summary_ms$site_orig)))

# Create data frame with NA for missing sites
contrast_table <- data.frame(
  Site = replace_sites(all_sites),
  MB_Estimate = NA,
  MB_SE = NA,
  MB_t = NA,
  MB_p = NA,
  MS_Estimate = NA,
  MS_SE = NA,
  MS_t = NA,
  MS_p = NA
)

# Fill in MB data
mb_site_match <- match(all_sites, as.character(contrast_summary_mb$site_orig))
for(i in 1:length(all_sites)) {
  if(!is.na(mb_site_match[i])) {
    contrast_table$MB_Estimate[i] <- contrast_summary_mb$estimate[mb_site_match[i]]
    contrast_table$MB_SE[i] <- contrast_summary_mb$SE[mb_site_match[i]]
    contrast_table$MB_t[i] <- contrast_summary_mb$t.ratio[mb_site_match[i]]
    contrast_table$MB_p[i] <- contrast_summary_mb$p.value[mb_site_match[i]]
  }
}

# Fill in MS data
ms_site_match <- match(all_sites, as.character(contrast_summary_ms$site_orig))
for(i in 1:length(all_sites)) {
  if(!is.na(ms_site_match[i])) {
    contrast_table$MS_Estimate[i] <- contrast_summary_ms$estimate[ms_site_match[i]]
    contrast_table$MS_SE[i] <- contrast_summary_ms$SE[ms_site_match[i]]
    contrast_table$MS_t[i] <- contrast_summary_ms$t.ratio[ms_site_match[i]]
    contrast_table$MS_p[i] <- contrast_summary_ms$p.value[ms_site_match[i]]
  }
}

# Format p-values
contrast_table$MB_p <- ifelse(is.na(contrast_table$MB_p), "—",
                              ifelse(contrast_table$MB_p < 0.001, "< 0.001", 
                                     sprintf("%.4f", contrast_table$MB_p)))
contrast_table$MS_p <- ifelse(is.na(contrast_table$MS_p), "—",
                              ifelse(contrast_table$MS_p < 0.001, "< 0.001", 
                                     sprintf("%.4f", contrast_table$MS_p)))

contrast_table <- order_by_elevation(contrast_table, term_col = "Site")

contrast_kable <- contrast_table %>%
  kbl(
    col.names = c("Site", "Estimate", "SE", "t", "p", 
                  "Estimate", "SE", "t", "p"),
    digits = 3,
    align = c("l", "r", "r", "r", "r", "r", "r", "r", "r"),
    escape = FALSE,
    row.names = FALSE
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    position = "left"
  ) %>%
  add_header_above(c(" " = 1, "M. boulderensis" = 4, "M. sanguinipes" = 4)) %>%
  row_spec(0, bold = TRUE) %>%
  column_spec(1, width = "2.5cm") %>%
  column_spec(2:9, width = "1.5cm")

print(contrast_kable)