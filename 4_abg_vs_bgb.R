library(tidyverse)
library(here)

(rm(list = ls()))

bgb_biomass_layer <- readRDS("output/bgb_biomass_layer.rds")
mixture <- readRDS("output/mixture_list.rds")
traits <- read.csv("data/Compiled_Traits_appended.csv")
comp_traits <- traits[traits$species=="comp",]

comp_traits$prop_scam_abg <- NA
for(i in 1:nrow(comp_traits)){
  comp_traits$p.scam.abg[i] <- comp_traits$agb_scam[i]/comp_traits$tot_bgb[i]
}
TOP <- bgb_biomass_layer[["TOP"]]
traits_TOP <- comp_traits[comp_traits$pot_no %in% TOP$pot_no,]

merged_TOP <- merge(TOP, traits_TOP, by = "pot_no", all.x = TRUE)


ggplot(data = merged_TOP, aes(x = p.scam.abg, y = p.scam)) +
  geom_point() +
  facet_wrap(~ env_treatment)+
  theme_bw() +
  ylim(0,1) +
  geom_abline(intercept = 0,slope = 1)


MID <- bgb_biomass_layer[["MID"]]
traits_MID <- comp_traits[comp_traits$pot_no %in% TOP$pot_no,]

merged_MID <- merge(MID, traits_MID, by = "pot_no", all.x = TRUE)


ggplot(data = merged_MID, aes(x = p.scam.abg, y = p.scam)) +
  geom_point() +
  facet_wrap(~ env_treatment)+
  theme_bw() +
  ylim(0,1) +
  geom_abline(intercept = 0,slope = 1)

BTM <- bgb_biomass_layer[["BTM"]]
traits_BTM <- comp_traits[comp_traits$pot_no %in% TOP$pot_no,]

merged_BTM <- merge(BTM, traits_BTM, by = "pot_no", all.x = TRUE)


ggplot(data = merged_BTM, aes(x = p.scam.abg, y = p.scam)) +
  geom_point() +
  facet_wrap(~ env_treatment)+
  theme_bw() +
  ylim(0,1) +
  geom_abline(intercept = 0,slope = 1)
