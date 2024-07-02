library(tidyverse)
library(here)
here::here()

bgb_biomass <- read_rds("output/bgb_biomass.rds")
traits <- read.csv("output/data/Compiled_Traits_appended.csv")
traits$env_treatment <- paste0(ifelse(traits$salinity == 0, "low_sal", "high_sal"),
                                         "/",
                                         ifelse(traits$elev_bi == 0, "low_elev", "high_elev"))
traits <- traits[-which(traits$species == "empty"), ]  # Modify the original data frame

# Filtered data (assuming missing values in tot_bgb and rsr are NA)
traits_filtered <- traits %>% 
  filter(!is.na(tot_bgb) & !is.na(rsr))

# Summarize by env_treatment and species
summary_data <- traits_filtered %>%
  group_by(env_treatment, species) %>%
  summarise(mean_tot_bgb = mean(tot_bgb),
            sd_tot_bgb = sd(tot_bgb),
            mean_rsr = mean(rsr),
            sd_rsr = sd(rsr))


ggplot(traits_filtered, aes(y = tot_bgb, x = factor(env_treatment), color = factor(species))) +
  geom_jitter(width = 0.2, height = 0.1) +  # Adjust jitter amount
  theme_bw() +
  ylim(0, 26) +
  geom_point(aes(x = env_treatment, y = mean_tot_bgb, size = 1.5), data = summary_data) +
  # geom_errorbar(aes(x = env_treatment, 
  #                   ymin = summary_data$mean_tot_bgb - summary_data$sd_tot_bgb, 
  #                   ymax = summary_data$mean_tot_bgb + summary_data$sd_tot_bgb), 
  #               linetype = "solid", size = 1, 
  #               color = c("#EE964B","#19647E", "#ED254E","#EE964B","#19647E", "#ED254E",
  #                         "#EE964B","#19647E", "#ED254E","#EE964B","#19647E", "#ED254E"), 
  #               data = summary_data,inherit.aes = FALSE) +
  labs(x = "Environmental Conditions", y = "Total Belowground Biomass (g)", color = "Competition") + # Custom labels
  scale_color_manual(values = c("#EE964B","#19647E", "#ED254E"),labels = c("comp", "scam","sppa"))+ # Set custom colors
  guides(size = guide_none())

ggplot(traits_filtered, aes(y= rsr, x = factor(env_treatment), color = factor(species))) +
  geom_jitter(width = 0.2, height = 0.1) +  # Adjust jitter amount
  theme_bw() +
  ylim(0,3.6) +
  labs(x = "Environmental Conditions", y = "Root-to-Shoot Ratio", color = "Competition") + # Custom labels
  scale_color_manual(values = c("#EE964B","#19647E", "#ED254E"),labels = c("Polyculture", "SCAM","SPPA")) +# Set custom colors
  geom_point(aes(x = env_treatment, y = mean_rsr, size = 1.5), data = summary_data) +
  guides(size = guide_none())
  
  
# Salinity
ggplot(traits_filtered, aes(y= tot_bgb, x = factor(salinity), color = factor(species))) +
  geom_jitter(width = 0.2, height = 0.1) +  # Adjust jitter amount
  theme_bw() +
  ylim(0,26) +
  scale_x_discrete(labels = c("Low", "High"))  +  # Custom labels
  labs(x = "Salinity", y = "Total Belowground Biomass (g)", color = "Competition") + # Custom labels
  scale_color_manual(values = c("#EE964B","#19647E", "#ED254E"),labels = c("Polyculture", "SCAM","SPPA")) # Set custom colors

ggplot(traits_filtered, aes(y= rsr, x = factor(salinity), color = factor(species))) +
  geom_jitter(width = 0.2, height = 0.1) +  # Adjust jitter amount
  theme_bw() +
  ylim(0,3.6) +
  scale_x_discrete(labels = c("Low", "High"))  +  # Custom labels
  labs(x = "Salinity", y = "Root-to-Shoot Ratio", color = "Competition") + # Custom labels
  scale_color_manual(values = c("#EE964B","#19647E", "#ED254E"),labels = c("Polyculture", "SCAM","SPPA")) # Set custom colors

# Inundation
ggplot(traits_filtered, aes(y= tot_bgb, x = factor(elev_bi), color = factor(species))) +
  geom_jitter(width = 0.2, height = 0.1) +  # Adjust jitter amount
  theme_bw() +
  ylim(0,26) +
  scale_x_discrete(labels = c("Low", "High"))  +  # Custom labels
  labs(x = "Elevation", y = "Total Belowground Biomass", color = "Competition") + # Custom labels
  scale_color_manual(values = c("#EE964B","#19647E", "#ED254E"),labels = c("Polyculture", "SCAM","SPPA")) # Set custom colors

ggplot(traits_filtered, aes(y= rsr, x = factor(salinity), color = factor(species))) +
  geom_jitter(width = 0.2, height = 0.1) +  # Adjust jitter amount
  theme_bw() +
  ylim(0,3.6) +
  scale_x_discrete(labels = c("Low", "High"))  +  # Custom labels
  labs(x = "Elevation", y = "Root-to-Shoot Ratio", color = "Competition") + # Custom labels
  scale_color_manual(values = c("#EE964B","#19647E", "#ED254E"),labels = c("Polyculture", "SCAM","SPPA")) # Set custom colors

ggplot(traits_filtered, aes(y= tot_bgb, x = factor(provenance), color = factor(species))) +
  geom_jitter(width = 0.2, height = 0.1) +  # Adjust jitter amount
  theme_bw() +
  ylim(0,26) +
  scale_x_discrete(labels = c("Blackwater", "Corn","Spartina monoculture"))  +  # Custom labels
  labs(x = "Provenance", y = "Total Belowground Biomass", color = "Competition") + # Custom labels
  scale_color_manual(values = c("#EE964B","#19647E", "#ED254E"),labels = c("Polyculture", "SCAM","SPPA")) # Set custom colors

ggplot(traits_filtered, aes(y= rsr, x = factor(provenance), color = factor(species))) +
  geom_jitter(width = 0.2, height = 0.1) +  # Adjust jitter amount
  theme_bw() +
  ylim(0,3.6) +
  scale_x_discrete(labels = c("Blackwater", "Corn", "Spartina monoculture"))  +  # Custom labels
  labs(x = "Provenance", y = "Root-to-Shoot Ratio", color = "Competition") + # Custom labels
  scale_color_manual(values = c("#EE964B","#19647E", "#ED254E"),labels = c("Polyculture", "SCAM","SPPA")) # Set custom colors


ggplot(bgb_biomass, aes(x = weight_density, y = prop_scam, color = env_treatment)) +
  geom_point() +  # Plot points for each data point
  geom_errorbar(aes(ymin = prop_scam_lower, ymax = prop_scam_upper), inherit.aes = TRUE) +  # Error bars using existing aesthetics
  facet_wrap(~ env_treatment) +  # Separate panels for each treatment
  labs(title = expression(Dry~Weight~Density~of~SCAM~(g/cm^{3})),  # Title with superscript
       x = "Total Weight Density", y = expression(Weight~Density~of~SCAM~(g/cm^{3}))) +
  theme_bw() + # Set a black and white theme (optional)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")  # Add 1:1 line


