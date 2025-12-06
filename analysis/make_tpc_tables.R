# create_tpc_tables.R
# Script to create formatted tables of TPC parameters with MAP + 90% HDI
library(tidyverse)
library(bayestestR)
library(coda)
library(kableExtra)

# Load the model fits
ms_fit <- readRDS("~/GitHub/thermal_perf/ms_improved7_betterprior.rds")
mb_fit <- readRDS("~/GitHub/thermal_perf/mb_both_years_with_yeareffect3betterpriors.rds")

# Set method
method <- "MAP"
CI <- 0.9

# Site elevations mapping
site_elevations <- c(
  "Eldo" = "1740m",
  "A1" = "2195m",
  "B1" = "2591m",
  "C1" = "3048m",
  "D1" = "3515m"
)

# ===== HELPER FUNCTION TO FORMAT VALUES =====
format_estimate <- function(point, lower, upper, digits = 2) {
  sprintf("%.*f (%.*f, %.*f)", 
          digits, point, 
          digits, lower, 
          digits, upper)
}

# ===== EXTRACT MB PARAMETERS =====

n_chains <- mb_fit$fit@sim$chains
mb_posteriors <- c()
for(i in 1:n_chains){
  b <- as.data.frame(as.mcmc(mb_fit)[[i]])
  b$chain <- i
  mb_posteriors <- rbind(b, mb_posteriors)
}

# Extract population-level parameters for MB
pops_mb <- c("A1", "B1", "C1", "D1")
mb_table_data <- data.frame()

for(x in pops_mb){
  # Extract parameters
  Tmin_col <- parse(text = paste0('mb_posteriors$`b_Tmin_pop', x, '`'), keep.source = FALSE)[[1]]
  Topt_col <- parse(text = paste0('mb_posteriors$`b_Topt_pop', x, '`'), keep.source = FALSE)[[1]]
  Tmax_col <- parse(text = paste0('mb_posteriors$`b_Above_pop', x, '` + mb_posteriors$`b_Topt_pop', x, '`'), keep.source = FALSE)[[1]]
  Ropt_col <- parse(text = paste0('mb_posteriors$`b_Ropt_pop', x, '`'), keep.source = FALSE)[[1]]
  
  row_data <- data.frame(
    Population = site_elevations[x],
    Tmin = format_estimate(
      as.numeric(point_estimate(eval(Tmin_col), centrality = method)),
      as.numeric(ci(eval(Tmin_col), ci = CI, method = "HDI")[2]),
      as.numeric(ci(eval(Tmin_col), ci = CI, method = "HDI")[3])
    ),
    Topt = format_estimate(
      as.numeric(point_estimate(eval(Topt_col), centrality = method)),
      as.numeric(ci(eval(Topt_col), ci = CI, method = "HDI")[2]),
      as.numeric(ci(eval(Topt_col), ci = CI, method = "HDI")[3])
    ),
    Tmax = format_estimate(
      as.numeric(point_estimate(eval(Tmax_col), centrality = method)),
      as.numeric(ci(eval(Tmax_col), ci = CI, method = "HDI")[2]),
      as.numeric(ci(eval(Tmax_col), ci = CI, method = "HDI")[3])
    ),
    Ropt = format_estimate(
      as.numeric(point_estimate(eval(Ropt_col), centrality = method)),
      as.numeric(ci(eval(Ropt_col), ci = CI, method = "HDI")[2]),
      as.numeric(ci(eval(Ropt_col), ci = CI, method = "HDI")[3]),
      digits = 3
    ),
    stringsAsFactors = FALSE
  )
  
  mb_table_data <- rbind(mb_table_data, row_data)
}

# Extract fixed effects for MB
mb_sex_effect <- format_estimate(
  as.numeric(point_estimate(mb_posteriors$`b_sexeffect_Intercept`, centrality = method)),
  as.numeric(ci(mb_posteriors$`b_sexeffect_Intercept`, ci = CI, method = "HDI")[2]),
  as.numeric(ci(mb_posteriors$`b_sexeffect_Intercept`, ci = CI, method = "HDI")[3]),
  digits = 3
)

mb_year_effect <- format_estimate(
  as.numeric(point_estimate(mb_posteriors$`b_yeareffectheight_Intercept`, centrality = method)),
  as.numeric(ci(mb_posteriors$`b_yeareffectheight_Intercept`, ci = CI, method = "HDI")[2]),
  as.numeric(ci(mb_posteriors$`b_yeareffectheight_Intercept`, ci = CI, method = "HDI")[3]),
  digits = 3
)

# Extract random effect SD for MB
mb_ropt_sd <- format_estimate(
  as.numeric(point_estimate(mb_posteriors$`sd_full_ID__Roptraneff_Intercept`, centrality = method)),
  as.numeric(ci(mb_posteriors$`sd_full_ID__Roptraneff_Intercept`, ci = CI, method = "HDI")[2]),
  as.numeric(ci(mb_posteriors$`sd_full_ID__Roptraneff_Intercept`, ci = CI, method = "HDI")[3]),
  digits = 3
)

# Create effect rows - headers first, then indented effects
mb_fixed_effects_header <- data.frame(
  Population = "Fixed Effects",
  Tmin = "",
  Topt = "",
  Tmax = "",
  Ropt = "",
  stringsAsFactors = FALSE
)

mb_sex_row <- data.frame(
  Population = "  Sex effect",
  Tmin = "",
  Topt = "",
  Tmax = "",
  Ropt = mb_sex_effect,
  stringsAsFactors = FALSE
)

mb_year_row <- data.frame(
  Population = "  Year effect",
  Tmin = "",
  Topt = "",
  Tmax = "",
  Ropt = mb_year_effect,
  stringsAsFactors = FALSE
)

mb_random_effects_header <- data.frame(
  Population = "Random Effects",
  Tmin = "",
  Topt = "",
  Tmax = "",
  Ropt = "",
  stringsAsFactors = FALSE
)

mb_individual_sd_row <- data.frame(
  Population = "  Individual SD",
  Tmin = "",
  Topt = "",
  Tmax = "",
  Ropt = mb_ropt_sd,
  stringsAsFactors = FALSE
)

# Combine all MB data with headers before effects
mb_table_complete <- rbind(
  mb_table_data,
  mb_fixed_effects_header,
  mb_sex_row,
  mb_year_row,
  mb_random_effects_header,
  mb_individual_sd_row
)

# ===== EXTRACT MS PARAMETERS =====

n_chains <- ms_fit$fit@sim$chains
ms_posteriors <- c()
for(i in 1:n_chains){
  b <- as.data.frame(as.mcmc(ms_fit)[[i]])
  b$chain <- i
  ms_posteriors <- rbind(b, ms_posteriors)
}

# Extract population-level parameters for MS
pops_ms <- c("Eldo", "A1", "B1")
ms_table_data <- data.frame()

for(x in pops_ms){
  # Extract parameters
  Tmin_col <- parse(text = paste0('ms_posteriors$`b_Tmin_pop', x, '`'), keep.source = FALSE)[[1]]
  Topt_col <- parse(text = paste0('ms_posteriors$`b_Topt_pop', x, '`'), keep.source = FALSE)[[1]]
  Tmax_col <- parse(text = paste0('ms_posteriors$`b_Above_pop', x, '` + ms_posteriors$`b_Topt_pop', x, '`'), keep.source = FALSE)[[1]]
  Ropt_col <- parse(text = paste0('ms_posteriors$`b_Ropt_pop', x, '`'), keep.source = FALSE)[[1]]
  
  row_data <- data.frame(
    Population = site_elevations[x],
    Tmin = format_estimate(
      as.numeric(point_estimate(eval(Tmin_col), centrality = method)),
      as.numeric(ci(eval(Tmin_col), ci = CI, method = "HDI")[2]),
      as.numeric(ci(eval(Tmin_col), ci = CI, method = "HDI")[3])
    ),
    Topt = format_estimate(
      as.numeric(point_estimate(eval(Topt_col), centrality = method)),
      as.numeric(ci(eval(Topt_col), ci = CI, method = "HDI")[2]),
      as.numeric(ci(eval(Topt_col), ci = CI, method = "HDI")[3])
    ),
    Tmax = format_estimate(
      as.numeric(point_estimate(eval(Tmax_col), centrality = method)),
      as.numeric(ci(eval(Tmax_col), ci = CI, method = "HDI")[2]),
      as.numeric(ci(eval(Tmax_col), ci = CI, method = "HDI")[3])
    ),
    Ropt = format_estimate(
      as.numeric(point_estimate(eval(Ropt_col), centrality = method)),
      as.numeric(ci(eval(Ropt_col), ci = CI, method = "HDI")[2]),
      as.numeric(ci(eval(Ropt_col), ci = CI, method = "HDI")[3]),
      digits = 3
    ),
    stringsAsFactors = FALSE
  )
  
  ms_table_data <- rbind(ms_table_data, row_data)
}

# Extract fixed effects for MS
ms_sex_effect <- format_estimate(
  as.numeric(point_estimate(ms_posteriors$`b_sexeffect_Intercept`, centrality = method)),
  as.numeric(ci(ms_posteriors$`b_sexeffect_Intercept`, ci = CI, method = "HDI")[2]),
  as.numeric(ci(ms_posteriors$`b_sexeffect_Intercept`, ci = CI, method = "HDI")[3]),
  digits = 3
)

# Extract random effect SD for MS
ms_ropt_sd <- format_estimate(
  as.numeric(point_estimate(ms_posteriors$`sd_full_ID__Roptraneff_Intercept`, centrality = method)),
  as.numeric(ci(ms_posteriors$`sd_full_ID__Roptraneff_Intercept`, ci = CI, method = "HDI")[2]),
  as.numeric(ci(ms_posteriors$`sd_full_ID__Roptraneff_Intercept`, ci = CI, method = "HDI")[3]),
  digits = 3
)

# Create effect rows
ms_fixed_effects_header <- data.frame(
  Population = "Fixed Effects",
  Tmin = "",
  Topt = "",
  Tmax = "",
  Ropt = "",
  stringsAsFactors = FALSE
)

ms_sex_row <- data.frame(
  Population = "  Sex effect",
  Tmin = "",
  Topt = "",
  Tmax = "",
  Ropt = ms_sex_effect,
  stringsAsFactors = FALSE
)

ms_random_effects_header <- data.frame(
  Population = "Random Effects",
  Tmin = "",
  Topt = "",
  Tmax = "",
  Ropt = "",
  stringsAsFactors = FALSE
)

ms_individual_sd_row <- data.frame(
  Population = "  Individual SD",
  Tmin = "",
  Topt = "",
  Tmax = "",
  Ropt = ms_ropt_sd,
  stringsAsFactors = FALSE
)

# Combine all MS data with headers before effects
ms_table_complete <- rbind(
  ms_table_data,
  ms_fixed_effects_header,
  ms_sex_row,
  ms_random_effects_header,
  ms_individual_sd_row
)

# ===== CREATE FORMATTED TABLES =====

# MB Table
mb_kable <- mb_table_complete %>%
  kbl(
    col.names = c("Population", "T<sub>min</sub> (°C)", "T<sub>opt</sub> (°C)", "T<sub>max</sub> (°C)", "R<sub>opt</sub> (mg/g/hr)"),
    align = c("l", "c", "c", "c", "c"),
    escape = FALSE,
    row.names = FALSE  # Remove row numbers
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    position = "left"
  ) %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(c(5, 8), bold = TRUE, italic = TRUE)  # Headers are bold and italic

# MS Table
ms_kable <- ms_table_complete %>%
  kbl(
    col.names = c("Population", "T<sub>min</sub> (°C)", "T<sub>opt</sub> (°C)", "T<sub>max</sub> (°C)", "R<sub>opt</sub> (mg/g/hr)"),
    align = c("l", "c", "c", "c", "c"),
    escape = FALSE,
    row.names = FALSE  # Remove row numbers
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    position = "left"
  ) %>%
  row_spec(0, bold = TRUE) %>%
  row_spec(c(4, 6), bold = TRUE, italic = TRUE)  # Headers are bold and italic

# Print tables
cat("\n=== M. boulderensis TPC Parameters ===\n")
print(mb_kable)

cat("\n\n=== M. sanguinipes TPC Parameters ===\n")
print(ms_kable)

# Create tables directory if it doesn't exist
dir.create("tables", showWarnings = FALSE)

# Save tables as HTML
save_kable(mb_kable, "tables/mb_tpc_parameters.html")
save_kable(ms_kable, "tables/ms_tpc_parameters.html")

# Also save as plain text for easy copying
sink("tables/mb_tpc_parameters.txt")
cat("M. boulderensis TPC Parameters (MAP with 90% HDI)\n\n")
print(mb_table_complete, row.names = FALSE)
sink()

sink("tables/ms_tpc_parameters.txt")
cat("M. sanguinipes TPC Parameters (MAP with 90% HDI)\n\n")
print(ms_table_complete, row.names = FALSE)
sink()

cat("\nTables created and saved to tables/ directory\n")