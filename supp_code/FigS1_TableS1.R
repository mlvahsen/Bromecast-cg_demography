# Create a map and table to source locations

## Preliminaries ####
library(tidyverse); library(usmap); library(geosphere);
library(ggnewscale); library(factoextra); library(patchwork)
# Read in common garden data
cg_data <- read_csv("data/cg_model_data.csv")
# Read in genotype and source location conversions
site2genotype <- read_csv("supp_data/sitecode2genotypenumber.csv")
# Read in GPS locations
gps <- read_csv("supp_data/gps_sites.csv")

# Get unique genotype information
cg_data %>% 
  select(genotype, PC1) %>% 
  distinct() -> cg_data

site2genotype %>% 
  mutate(genotype = parse_number(genotypeID)) %>% 
  select(genotype, site = Site.code) -> site2genotype

merge(cg_data, site2genotype) %>% 
  merge(gps) %>% 
  distinct() -> plot_data

## Make plot of source locations ####
# Get US map with states
state <- map_data("state")
# Get Canada map
canada <- map_data("world") %>% filter(region == "Canada")
# Join into single region
region <- rbind(state, canada)

# Plot of PC1 by site for genotypes
ggplot(data=region, aes(x=long, y=lat, group = group)) +
  geom_polygon(color = "gray", fill = "white")+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  coord_cartesian(ylim = c(34, 51), xlim = c(-125, -102)) +
  geom_point(data = plot_data, aes(x = lon, y = lat, group = NA, fill = PC1),
             color = "black", shape = 21, size = 4) +
  theme_classic(base_size = 18) +
  scale_fill_distiller(palette = "RdYlBu", direction = 1) +
  labs(fill = "PC 1",
       x = "Longitude",
       y = "Latitude",
       shape = "Common garden") +
  scale_shape_manual(values = c(2,5,0,6)) -> plot

png("figs/FigS1_SourceLocations.png", height = 8, width = 6.2, units = "in", res = 300)
plot & theme(legend.position = "bottom",
           legend.key.width=unit(1.5,"cm")) &
  labs(fill = "PC 1 (hot & dry → cool & wet)") &
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))
dev.off()
