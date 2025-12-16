library(tidyverse)
library(ggplot2)
library(grid)
library(gtable)

setwd("~/GitHub/EB_Ch1/analysis")

thing0 <- read_csv("eb_results_summary_wide.csv") 
pops <- readRDS("C:/Users/jmsmi/OneDrive/Documents/GitHub/EB_Ch1/data/pops.rds")

by <- join_by(species==spp, site_orig==site, sex)
thing <- left_join(thing0, pops, by)

# Define color schemes
site_colors <- c(
  "Eldo" = "#FF4500",
  "A1" = "#FF8C00",
  "B1" = "#FFD700",
  "C1" = "#4682B4",
  "D1" = "#0000CD"
)

site_elevations <- c(
  "Eldo" = "1740m",
  "A1" = "2195m",
  "B1" = "2591m",
  "C1" = "3048m",
  "D1" = "3515m"
)

#######################################################
# FIGURE 3: Sex-Averaged Reciprocal Transplant - Contemporary Only
# Colored by population, solid fill for home, dots for transplant
#######################################################

create_sex_averaged_reciprocal_transplant_boxplot <- function(data, use_log_scale = FALSE) {
  site_order <- c("Eldo", "A1", "B1", "C1", "D1")
  
  averaged_data_combined <- data %>%
    filter(year_period == "contemporary") %>%
    group_by(species, year, site_orig, site_clim) %>%
    summarize(
      avg_energy_per_mass = mean(shade_NA_h_NA / mass, na.rm = TRUE),
      n_sex = n(),
      .groups = "drop"
    ) %>%
    mutate(
      site_orig = factor(site_orig, levels = site_order),
      site_clim = factor(site_clim, levels = site_order),
      climate_type = ifelse(site_orig == site_clim, "Home climate", "Transplanted climate"),
      population_label = site_elevations[as.character(site_orig)],
      population_label = factor(population_label, levels = site_elevations),
      site_clim_label = site_elevations[as.character(site_clim)],
      site_clim_label = factor(site_clim_label, levels = site_elevations)
    )
  
  if (use_log_scale) {
    averaged_data_combined <- averaged_data_combined %>%
      mutate(avg_energy_per_mass = log(avg_energy_per_mass))
  }
  
  # Split data for home vs transplanted
  home_data <- averaged_data_combined %>% filter(climate_type == "Home climate")
  transplant_data <- averaged_data_combined %>% filter(climate_type == "Transplanted climate")
  
  # Create elevation-based color palette
  elev_colors <- setNames(site_colors, site_elevations)
  
  p <- ggplot(averaged_data_combined, 
              aes(x = population_label, 
                  y = avg_energy_per_mass)) +
    facet_grid(species ~ site_clim_label, scales = "free_y",
               labeller = labeller(
                 species = as_labeller(c("MB" = "bolditalic('M. boulderensis')", 
                                         "MS" = "bolditalic('M. sanguinipes')"), 
                                       label_parsed))) +
    # Boxplots with consistent alpha
    geom_boxplot(aes(fill = population_label),
                 alpha = 0.4, 
                 color = "gray30",
                 outlier.shape = NA) +
    # Jittered points for transplanted (X shape)
    geom_jitter(data = transplant_data,
                aes(shape = "Transplanted"),
                color = "gray30", 
                width = 0.2, 
                alpha = 0.5, 
                size = 1.5,
                stroke = 0.7) +
    # Jittered points for home (filled dots)
    geom_jitter(data = home_data,
                aes(shape = "Home"),
                color = "gray30",
                width = 0.2, 
                alpha = 0.5, 
                size = 1.5) +
    scale_fill_manual(values = elev_colors, name = "Population") +
    scale_shape_manual(values = c("Home" = 16, "Transplanted" = 4),
                       name = "Climate type") +
    scale_x_discrete(drop = FALSE, limits = site_elevations) +
    labs(
      x = "Population",
      y = if(use_log_scale) "log(Total Net Energy Gained (kJ/g))" else "Total Net Energy Gained (kJ/g)"
    ) +
    theme_bw() +
    theme(
      legend.position = c(0.95, 0.08),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.margin = margin(t = 4, r = 4, b = 4, l = 4),
      legend.key.size = unit(0.8, "lines"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9, face = "bold"),
      legend.spacing.y = unit(0.1, "cm"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "#D2B48C"),
      strip.text = element_text(face = "bold", color = "black"),
      panel.grid.minor = element_blank()
    ) +
    guides(
      fill = guide_legend(order = 1, keyheight = unit(0.7, "lines"), keywidth = unit(0.7, "lines")),
      shape = guide_legend(order = 2, keyheight = unit(0.7, "lines"), keywidth = unit(0.7, "lines"))
    )
  
  # Convert to grob for manipulation
  g <- ggplotGrob(p)
  
  # Add "Climate" label at the top
  climate_label <- textGrob("Climate", gp = gpar(fontsize = 12, fontface = "bold", col = "#8B4513"))
  g <- gtable_add_rows(g, heights = unit(0.6, "cm"), pos = 0)
  panel_cols <- which(grepl("panel", g$layout$name))
  g <- gtable_add_grob(g, climate_label, t = 1, 
                       l = min(g$layout$l[panel_cols]),
                       r = max(g$layout$r[panel_cols]))
  
  # Add "Species" label on the right
  species_label <- textGrob("Species", rot = 270, gp = gpar(fontsize = 12, fontface = "bold", col = "#8B4513"))
  g <- gtable_add_cols(g, widths = unit(0.6, "cm"), pos = -1)
  panel_rows <- range(which(grepl("panel", g$layout$name)))
  g <- gtable_add_grob(g, species_label, 
                       t = min(g$layout$t[panel_rows]),
                       b = max(g$layout$b[panel_rows]),
                       l = ncol(g), r = ncol(g))
  
  return(g)
}

#######################################################
# FIGURE 4: Historical vs Contemporary Comparison
# Points aligned with boxplots using position_jitterdodge
#######################################################

create_sex_averaged_historical_contemporary_boxplot <- function(data, use_log_scale = FALSE) {
  site_order <- c("Eldo", "A1", "B1", "C1", "D1")
  
  averaged_data_combined <- data %>%
    filter(site_orig == site_clim) %>%
    group_by(species, year, year_period, site_orig) %>%
    summarize(
      avg_energy_per_mass = mean(shade_NA_h_NA / mass, na.rm = TRUE),
      n_sex = n(),
      .groups = "drop"
    ) %>%
    mutate(
      site = site_orig,
      site = factor(site, levels = site_order),
      year_period = factor(year_period, levels = c("historical", "contemporary")),
      site_label = site_elevations[as.character(site)],
      site_label = factor(site_label, levels = site_elevations)
    )
  
  if (use_log_scale) {
    averaged_data_combined <- averaged_data_combined %>%
      mutate(avg_energy_per_mass = log(avg_energy_per_mass))
  }
  
  p <- ggplot(averaged_data_combined, 
              aes(x = site_label, 
                  y = avg_energy_per_mass, 
                  fill = year_period)) +
    facet_grid(species ~ ., scales = "free_y",
               labeller = labeller(
                 species = as_labeller(c("MB" = "bolditalic('M. boulderensis')", 
                                         "MS" = "bolditalic('M. sanguinipes')"), 
                                       label_parsed))) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 0.75)) +
    geom_point(color = "gray30",
               position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
               alpha = 0.5, 
               size = 1) +
    scale_fill_manual(values = c("historical" = "#8DA0CB", "contemporary" = "#FC8D62"),
                      labels = c("Historical (1950-1959)", "Contemporary (2014-2024)"),
                      name = "Time Period") +
    scale_x_discrete(limits = site_elevations) +
    labs(
      x = "Site",
      y = if(use_log_scale) "log(Total Net Energy Gained (kJ/g))" else "Total Net Energy Gained (kJ/g)"
    ) +
    theme_bw() +
    theme(
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      legend.box.margin = margin(6, 6, 6, 6),
      legend.key.size = unit(0.8, "lines"),
      legend.text = element_text(size = 9),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#D2B48C"),
      strip.text = element_text(face = "bold", color = "black")
    ) +
    guides(fill = guide_legend(keyheight = unit(0.8, "lines"),
                               keywidth = unit(1, "lines")))
  
  # Convert to grob for manipulation
  g <- ggplotGrob(p)
  
  # Add "Species" label on the right
  species_label <- textGrob("Species", rot = 270, 
                            gp = gpar(fontsize = 12, fontface = "bold", col = "#8B4513"))
  g <- gtable_add_cols(g, widths = unit(0.6, "cm"), pos = -1)
  panel_rows <- which(grepl("panel", g$layout$name))
  g <- gtable_add_grob(g, species_label, 
                       t = min(g$layout$t[panel_rows]),
                       b = max(g$layout$b[panel_rows]),
                       l = ncol(g), r = ncol(g))
  
  return(g)
}

# Generate and display plots
plot1 <- create_sex_averaged_reciprocal_transplant_boxplot(thing, use_log_scale = FALSE)
grid.newpage()
grid.draw(plot1)

plot2 <- create_sex_averaged_historical_contemporary_boxplot(thing, use_log_scale = FALSE)
grid.newpage()
grid.draw(plot2)

# Save plots
ggsave("fig3_reciprocal_transplant.png", plot1, width = 12, height = 7, dpi = 300)
ggsave("fig3_reciprocal_transplant.pdf", plot1, width = 12, height = 7)

ggsave("fig4_historical_contemporary.png", plot2, width = 8, height = 6, dpi = 300)
ggsave("fig4_historical_contemporary.pdf", plot2, width = 8, height = 6)