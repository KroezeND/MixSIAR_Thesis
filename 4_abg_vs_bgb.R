(rm(list = ls()))

bgb_biomass_layer <- readRDS("output/bgb_biomass_layer.rds")
traits <- read.csv("data/Compiled_Traits_appended.csv")
comp_traits <- traits[traits$species=="comp",]

comp_traits$prop_scam_abg <- NA
for(i in 1:nrow(comp_traits)){
  comp_traits$p.scam.abg[i] <- comp_traits$agb_scam[i]/comp_traits$tot_agb[i]
}
TOP <- bgb_biomass_layer[["TOP"]]
traits_TOP <- comp_traits[comp_traits$pot_no %in% TOP$pot_no,]

merged_TOP <- merge(TOP, traits_TOP, by = "pot_no", all.x = TRUE)


ggplot2::ggplot(data = merged_TOP, ggplot2::aes(x = p.scam.abg, y = p.scam)) +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~ env_treatment)+
  ggplot2::theme_bw() +
  ggplot2::ylim(0,1) +
  ggplot2::geom_abline(intercept = 0,slope = 1)


MID <- bgb_biomass_layer[["MID"]]
traits_MID <- comp_traits[comp_traits$pot_no %in% TOP$pot_no,]

merged_MID <- merge(MID, traits_MID, by = "pot_no", all.x = TRUE)


ggplot2::ggplot(data = merged_MID, ggplot2::aes(x = p.scam.abg, y = p.scam)) +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~ env_treatment)+
  ggplot2::theme_bw() +
  ggplot2::ylim(0,1) +
  ggplot2::geom_abline(intercept = 0,slope = 1)

BTM <- bgb_biomass_layer[["BTM"]]
traits_BTM <- comp_traits[comp_traits$pot_no %in% TOP$pot_no,]

merged_BTM <- merge(BTM, traits_BTM, by = "pot_no", all.x = TRUE)


ggplot2::ggplot(data = merged_BTM, ggplot2::aes(x = p.scam.abg, y = p.scam)) +
  ggplot2::geom_point() +
  ggplot2::facet_wrap(~ env_treatment)+
  ggplot2::theme_bw() +
  ggplot2::ylim(0,1) +
  ggplot2::geom_abline(intercept = 0,slope = 1)
