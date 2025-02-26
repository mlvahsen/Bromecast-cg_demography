
# Load libraries
library(rjags);library(tidyverse);library(cmdstanr);library(posterior);
library(bayesplot); library(janitor); library(patchwork); library(lubridate); 
library(loo)
# Source data for modeling 
source("supp_code/data_prep.R")

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
         year = substr(site_year_gravel, 6, 7),
         site_year = paste(site, year, sep = "\n"),
         gravel = substr(site_year_gravel, 9, 13)) %>% 
  ggplot(aes(x = temp_fecun, y = 0)) +
  annotate("segment",x=3,xend=15, y=0, yend=0, linewidth=2) +
  annotate("segment",x=3,xend=3, y=-0.1,yend=0.1, linewidth=2) +
  annotate("segment",x=15,xend=15, y=-0.1,yend=0.1, linewidth=2) +
  annotate("text", x = 3, y = 0, label = "3", vjust = 2, size = 5) +
  annotate("text", x = 15, y = 0, label = "15", vjust = 2, size = 5) +
  geom_point(size = 15, shape = 21, aes(color = site, fill = gravel), stroke = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = site_year), col="black", nudge_y = 1, fontface = "bold") +
  scale_x_continuous(limits = c(3,15)) +
  scale_y_continuous(limits = c(-1,1)) +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("black", "white")) +
  scale_color_manual(values = c("#44AA99", "#AA4499", "#332288", "#6699CC")) +
  ggtitle("Temperature (°C)")-> a


cg_model_fecun %>% 
  dplyr::select(site_year_gravel,temp_fecun, vwc_avg) %>%
  distinct() %>% 
  mutate(site = substr(site_year_gravel, 1, 2),
         year = substr(site_year_gravel, 6, 7),
         site_year = paste(site, year, sep = "\n"),
         gravel = substr(site_year_gravel, 9, 13)) %>% 
  ggplot(aes(x = vwc_avg, y = 0, color = color)) +
  annotate("segment",x=0.10,xend=0.25, y=0, yend=0, size=2) +
  annotate("segment",x=0.10,xend=0.10, y=-0.1,yend=0.1, size=2) +
  annotate("segment",x=0.25,xend=0.25, y=-0.1,yend=0.1, size=2) +
  annotate("text", x = 0.10, y = 0, label = "10", vjust = 2, size = 5) +
  annotate("text", x = 0.25, y = 0, label = "25", vjust = 2, size = 5) +
  geom_point(size = 15, shape = 21, aes(color = site, fill = gravel), stroke = 3, alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = site_year), col="black", nudge_y = 1, fontface = "bold") +
  scale_x_continuous(limits = c(0.10,0.25)) +
  scale_y_continuous(limits = c(-1,1)) +
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = c("black", "white")) +
  scale_color_manual(values = c("#44AA99", "#AA4499", "#332288", "#6699CC")) +
  ggtitle("Soil moisture (% vwc)") -> b

png("figs/Fig1_microclimate.png", height = 3.9, width = 9.2, res = 300, units = "in")
a/b 
dev.off()
