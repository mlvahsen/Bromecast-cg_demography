# This code brings together all components of Figure 6. Subcomponents were
# generated from 'Fig6_dens_dep_climate.R' and 'Fig6_densxlocaladapt.R'.

# Load libraries
library(tidyverse)

# Read in rds files of subplots
survival_temp <- read_rds("outputs/survival_temp.rds")
fecundity_temp <- read_rds("outputs/fecundity_temp.rds")
fitness_temp <- read_rds("outputs/fitness_temp.rds")

survival_vwc <- read_rds("outputs/survival_vwc.rds")
fecundity_vwc<- read_rds("outputs/fecundity_vwc.rds")
fitness_vwc <- read_rds("outputs/fitness_vwc.rds")

survival_sub <- read_rds("outputs/survival_sub.rds")
fecundity_sub<- read_rds("outputs/fecundity_sub.rds")
fitness_sub <- read_rds("outputs/fitness_sub.rds")

# Set up different themes for plots (gridded for local adaptation)
eco_theme <- theme_classic(base_size = 16)
evo_theme <- theme_bw(base_size = 16)

# Helpers to remove duplicate labels
no_y_title <- theme(axis.title.y = element_blank())
no_y_ticks <- theme(axis.text.y  = element_blank())
no_x_title <- theme(axis.title.x = element_blank())

# Row-specific y scales
y_surv <- scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05)))
y_fec  <- scale_y_continuous(limits = c(2, 8), expand = expansion(mult = c(0, 0.05)))
y_fit  <- scale_y_continuous(limits = c(1, 6.2), expand = expansion(mult = c(0, 0.05)))

dens_combined_plot <-
  
  # Row 1 — Survival
  (survival_temp + eco_theme + y_surv + no_x_title) +
  (survival_vwc  + eco_theme + y_surv + no_y_title + no_y_ticks + no_x_title) +
  (survival_sub  + evo_theme + y_surv + no_y_title + no_y_ticks + no_x_title) +
  
  # Row 2 — Fecundity
  (fecundity_temp + eco_theme + y_fec + no_x_title) +
  (fecundity_vwc  + eco_theme + y_fec + no_y_title + no_y_ticks + no_x_title) +
  (fecundity_sub  + evo_theme + y_fec + no_y_title + no_y_ticks + no_x_title) +
  
  # Row 3 — Fitness
  (fitness_temp  + eco_theme + y_fit) +
  (fitness_vwc   + eco_theme + y_fit + no_y_title + no_y_ticks) +
  (fitness_sub   + evo_theme + y_fit + no_y_title + no_y_ticks) +
  
  plot_layout(
    nrow = 3,
    ncol = 3,
    guides = "collect",
    axis_titles = "collect"
  ) +
  
  plot_annotation(tag_levels = "a") &
  
  theme(
    legend.position = "right",
    legend.background = element_blank()
  )

# Save
png("figs/Fig6_all.png", height = 10, width = 12, res = 300, units = "in")
dens_combined_plot
dev.off()
