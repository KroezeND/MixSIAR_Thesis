library(tidyverse)
library(here)
here::here()

source("src/fit_lm.R")
bgb_biomass <- read_rds("output/bgb_biomass.rds")

# Get unique seed_years
unique_years <- as.character(sort(unique(bgb_biomass$seed_year)))

# Generate data for each year
year_data <- list()
for (year in unique_years) {
  # Update year_data with predictions for each year
  year_data <- predict_and_store(bgb_biomass, year_data, year)
  rm(year)
}

combined_data <- do.call(rbind, year_data)
# Define a color palette 
colors <- c("#c51b7d", "#de77ae", "#f1b6da", "#e6f5d0", "#b8e186", "#7fbc41", "#4d9221","darkgray")
# Corn Ancestral:
# 1872, 1920, 1938, 
# #Blackwater Ancestral:
# # 1979, 1985
# Blackwater Modern
# 2011, 2014, 2017
# Corn Modern 
# 2017

color_map <- setNames(colors, unique_years)

ggplot2::ggplot() +
  geom_line(data = combined_data, aes(x = x, y = y, color = year), linewidth = 0.5) +  # Adjust line size, color by interaction
  geom_jitter(data = bgb_biomass, aes(x = midpoint, y = log(weight_density), color = factor(seed_year)), alpha = 0.5,width = 0.3) +  # Adjust point transparency
  # facet_wrap(~ seed_year) +  # Facet by seed_year
  labs(title = "Depth Layer Midpoint vs. Log Transformed Weight Density Total",
       x = "Depth Layer Midpoint",
       y = "Log(Weight Density Total)") +
  scale_color_manual(values = color_map) +  # Set color based on year
  theme_bw()  # Apply black and white theme


unique_treatment <- unique(bgb_biomass$env_treatment)
env_list <- list()
for(env in unique_treatment){
  unique_years <- as.character(sort(unique(bgb_biomass$seed_year[bgb_biomass$env_treatment==env])))
  env_list[[env]] <- list()
  for(year in unique_years){
    # Update year_data with predictions for each year
    env_list[[env]] <- predict_and_store(bgb_biomass[bgb_biomass$env_treatment==env,], env_list[[env]], year)
  }
}

g <- list()
for(i in names(env_list)){
  treatment <- env_list[[i]]
  treatment <- do.call(rbind, treatment)
  g[[i]] <- ggplot2::ggplot() +
    geom_line(data = treatment, aes(x = x, y = y, color = year), linewidth = 0.5) +  # Adjust line size, color by interaction
    geom_jitter(data = bgb_biomass[bgb_biomass$env_treatment==i,], 
                aes(x = midpoint, y = log(weight_density), 
                    color = factor(seed_year)),alpha = 0.5,width = 0.3) +  # Adjust point transparency
         labs(x = "Depth Layer Midpoint",
         y = "Log(Weight Density Total)") +
    scale_color_manual(values = color_map) +  # Set color based on year
    ylim(-11,-3) +
    theme_bw()  # Apply black and white theme
}

g$`high_sal/high_elev` + labs(title = "High Salinity/High Elevation")
g$`low_sal/high_elev` + labs(title = "Low Salinity/High Elevation")
g$`high_sal/low_elev` + labs(title = "High Salinity/Low Elevation")
g$`low_sal/low_elev` + labs(title = "Low Salinity/Low Elevation")
