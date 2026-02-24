
# Load libraries
library(rjags);library(tidyverse);library(cmdstanr);library(posterior);
library(bayesplot); library(janitor); library(patchwork); library(lubridate); 
library(loo); library(maps); library(ggrepel)


# Source data for modeling 
#source("supp_code/data_prep.R")
cg_model <- read_csv("data/cg_model_data.csv")# Updated 3 April 2025

# Make data set of just plants that survived and reproduced
cg_model$survived <- ifelse((cg_model$inflor_mass) > 0, 1, 0)

# Make subset of data that survived
cg_model %>% 
  filter(survived == 1) -> cg_model_fecun

# Make site year gravel index
cg_model_fecun$site_year_gravel <- paste(cg_model_fecun$site, cg_model_fecun$year, cg_model_fecun$albedo, sep = "_")

cg_model_fecun %>% 
  dplyr::select(site_year_gravel,temp_fecun, vwc_avg) %>%
  distinct() %>% 
  mutate(site = substr(site_year_gravel, 1, 2),
         year = paste(20, substr(site_year_gravel, 6, 7), sep = ""),
         site_year = paste(site, year, sep = "\n"),
         gravel = substr(site_year_gravel, 9, 13)) -> plot_data

plot_data %>% 
  # Change site names to include temp/moisture indicators
  mutate(site = case_when(site == "SS" ~ "cold (SS)",
                             site == "BA" ~ "cool & dry (BA)",
                             site == "CH" ~ "cool & wet (CH)",
                             site == "WI" ~ "hot (WI)")) %>% 
  mutate(site = factor(site, levels = c("cold (SS)", "cool & wet (CH)",
                                        "cool & dry (BA)", "hot (WI)"))) %>% 
  ggplot(aes(x = temp_fecun, y = vwc_avg)) +
  geom_point(size = 6, aes(color = site, fill = gravel, shape = year), stroke = 3) +
  scale_fill_manual(values = c("black", "white")) +
  scale_color_manual(values = c("#2b83ba","#abdda4","#fdae61", "#d7191c")) +
  scale_shape_manual(values = c(21,22)) +
  theme_classic(base_size = 16) +
  labs(y = "average daily soil moisture\n(volumetric water content)",
       x = "average daily soil temperature (°C)") +
  guides(colour = guide_legend(override.aes = list(size=5, stroke = 2)),
         fill = guide_legend(override.aes = list(size=5, stroke = 2, shape = 21)),
         shape = guide_legend(override.aes = list(size=5, stroke = 2)))-> xy_plot

# Sample GPS coordinates in Wyoming and Idaho
gps_data <- data.frame(
  site = c("Baltzor, ID", "Wildcat, ID", "Sheep Station, ID", "Cheyenne, WY"),
  lat = c(43.208482, 43.47437, 44.24559, 41.212078),
  lon = c(-116.995198, -116.90177, -112.214337, -104.852543)
)

# Get US map data
states_map <- map_data("state")

# Plot the map
ggplot() +
  geom_polygon(data = states_map, aes(x = long, y = lat, group = group),
               fill = "gray90", color = "white") +
  geom_point(data = gps_data, aes(x = lon, y = lat, color = site), size = 6,
             shape = 16, stroke = 2) +
  ggrepel::geom_text_repel(data = gps_data, aes(x = lon, y = lat, label = site, color = site),
                  size = 6, box.padding = 1.5) +  # Increased point_padding
  coord_quickmap(xlim = c(-117, -104), ylim = c(41, 46)) +
  theme_minimal(base_size = 16) +
  labs(x = "longitude", y = "latitude") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#fdae61", "#abdda4", "#2b83ba", "#d7191c")) +
  annotate(geom = "text", x = -114.4, y = 41.7, label = "(1573 m)", color = "#fdae61", size = 5)+ 
  annotate(geom = "text", x = -109, y = 44.8, label = "(1647 m)", color = "#2b83ba", size = 5)+ 
  annotate(geom = "text", x = -107.7, y = 41.6, label = "(1915 m)", color = "#abdda4", size = 5)+ 
  annotate(geom = "text", x = -114.6, y = 43.8, label = "(823 m)", color = "#d7191c", size = 5)-> map


design <- "BBBAAA
           BBBAAA
           BBBAAA
           CCCAAA
           CCCAAA"

fig2 <- xy_plot + map + plot_spacer() + plot_layout(design = design)
ggsave(filename = "figs/Fig2_setup.svg", plot = fig2, width = 12.88, height = 6.66)

