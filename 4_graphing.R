library(tidyverse)
library(here)

(rm(list = ls()))

bgb_biomass_mod <- read_rds("output/bgb_biomass_mod.rds")
traits <- read.csv("data/Compiled_Traits_appended.csv")
traits$env_treatment <- paste0(ifelse(traits$salinity == 0, "low_sal", "high_sal"),
                                         "/",
                                         ifelse(traits$elev_bi == 0, "low_elev", "high_elev"))
traits <- traits[-which(traits$species == "empty"), ]  # Modify the original data frame

# Filtered data (assuming missing values in tot_bgb are NA)
traits_filtered <- traits %>% 
  filter(!is.na(tot_bgb) & !is.na(rsr))

ggplot(traits_filtered, aes(y = tot_bgb, x = factor(env_treatment), fill = factor(species))) +
  geom_violin(position=position_dodge(0.90),scale = 'width',trim = T) + 
  theme_bw() +
  ylim(0, 26) +
  labs(x = "Environmental Conditions", y = "Total Belowground Biomass (g)", color = "Competition") + # Custom labels
  scale_fill_manual(values = c("#EE964B","#19647E", "#ED254E"),labels = c("COMP", "SCAM","SPPA"))+ # Set custom colors
  guides(fill = guide_legend(title = "Species"))
  
# Salinity
ggplot(traits_filtered, aes(y = tot_bgb, x = factor(salinity), fill = factor(species))) +
  geom_violin(position=position_dodge(0.90),scale = 'width',trim = T) + 
  theme_bw() +
  ylim(0, 26) +
  scale_x_discrete(labels = c("Low", "High"))  +  # Custom labels
  labs(x = "Environmental Conditions", y = "Total Belowground Biomass (g)", color = "Competition") + # Custom labels
  scale_fill_manual(values = c("#EE964B","#19647E", "#ED254E"),labels = c("COMP", "SCAM","SPPA"))+ # Set custom colors
  guides(fill = guide_legend(title = "Species"))

# Inundation
ggplot(traits_filtered, aes(y = tot_bgb, x = factor(elev_bi), fill = factor(species))) +
  geom_violin(position=position_dodge(0.90),scale = 'width',trim = T) + 
  theme_bw() +
  ylim(0, 26) +
  scale_x_discrete(labels = c("Low", "High"))  +  # Custom labels
  labs(x = "Environmental Conditions", y = "Total Belowground Biomass (g)", color = "Competition") + # Custom labels
  scale_fill_manual(values = c("#EE964B","#19647E", "#ED254E"),labels = c("COMP", "SCAM","SPPA"))+ # Set custom colors
  guides(fill = guide_legend(title = "Species"))

ggplot(traits_filtered, aes(y = tot_bgb, x = factor(provenance), fill = factor(species))) +
  geom_violin(position=position_dodge(0.90),scale = 'width',trim = T) + 
  theme_bw() +
  ylim(0, 26) +
  scale_x_discrete(labels = c("Blackwater", "Corn","SPPA"))  +  # Custom labels
  labs(x = "Environmental Conditions", y = "Total Belowground Biomass (g)", color = "Competition") + # Custom labels
  scale_fill_manual(values = c("#EE964B","#19647E", "#ED254E"),labels = c("COMP", "SCAM","SPPA"))+ # Set custom colors
  guides(fill = guide_legend(title = "Species"))


