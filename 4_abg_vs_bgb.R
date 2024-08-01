(rm(list = ls()))

bgb_biomass_layer <- readRDS("output/bgb_biomass_layer.rds")
traits <- read.csv("data/Compiled_Traits_appended.csv")
comp_traits <- traits[traits$species=="comp",]

comp_traits$p.scam.abg <- NA
for(i in 1:nrow(comp_traits)){
  comp_traits$p.scam.abg[i] <- comp_traits$agb_scam[i]/comp_traits$tot_agb[i] # measured abg scam proportion
}
TOP <- bgb_biomass_layer[["TOP"]] # estimated
traits_TOP <- comp_traits[comp_traits$pot_no %in% TOP$pot_no,] 
traits_TOP$env_treatment <- paste0(ifelse(traits_TOP$salinity == 0, "Rhode River", "GCREW"),
                                    "/",
                                    ifelse(traits_TOP$elev_bi == 0, "high inundation", "low inundation"))
merged_TOP <- merge(TOP, traits_TOP, by = "pot_no", all.x = TRUE)

average_p_scam_TOP <- merged_TOP |>
  group_by(env_treatment.x, gt, seed_year.x, cohort) |>
  summarize(
    mean_p_scam_abg = mean(p.scam.abg, na.rm = TRUE),
    mean_p_scam_bgb = mean(p.scam, na.rm = TRUE),
    p_scam_bgb_lower = mean(p.scam_lower, na.rm = TRUE),
    p_scam_bgb_upper = mean(p.scam_upper, na.rm = TRUE),
    .groups = "keep")
  
unique_cohort <- sort(unique(merged_TOP$cohort))
colors <- c("#2d6a4f","#74c69d","#5a189a","#c77dff")
color_map <- setNames(colors, unique_cohort)

ggTOP <- ggplot2::ggplot(data = average_p_scam_TOP, ggplot2::aes(y = mean_p_scam_bgb, x = mean_p_scam_abg, color = cohort)) +
  ggplot2::geom_point() +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = p_scam_bgb_lower, ymax = p_scam_bgb_upper),width = 3,linewidth = 0.5)+
  ggplot2::facet_wrap(~ env_treatment.x)+
  ggplot2::theme_bw() +
  ggplot2::scale_color_manual(values = color_map,
                              labels = c("blackwater\nancestral",
                                         "blackwater\nmodern",
                                         "corn\nancestral",
                                         "corn\nmodern")) +
  ggplot2::ylim(0,1) +
  ggplot2::xlim(0,1) +
  ggplot2::labs(x = "", y = "Estimated Bgb p.scam TOP") +
  ggplot2::geom_abline(intercept = 0,slope = 1, linetype = "dashed", color = "red") +
  ggplot2::theme(legend.position = "none")


MID <- bgb_biomass_layer[["MID"]] # estimated
traits_MID <- comp_traits[comp_traits$pot_no %in% MID$pot_no,] 
traits_MID$env_treatment <- paste0(ifelse(traits_MID$salinity == 0, "Rhode River", "GCREW"),
                                   "/",
                                   ifelse(traits_MID$elev_bi == 0, "high inundation", "low inundation"))
merged_MID <- merge(MID, traits_MID, by = "pot_no", all.x = TRUE)

average_p_scam_MID <- merged_MID |>
  group_by(env_treatment.x, gt, seed_year.x, cohort) |>
  summarize(
    mean_p_scam_abg = mean(p.scam.abg, na.rm = TRUE),
    mean_p_scam_bgb = mean(p.scam, na.rm = TRUE),
    p_scam_bgb_lower = mean(p.scam_lower, na.rm = TRUE),
    p_scam_bgb_upper = mean(p.scam_upper, na.rm = TRUE),
    .groups = "keep")

ggMID <- ggplot2::ggplot(data = average_p_scam_MID, ggplot2::aes(y = mean_p_scam_bgb, x = mean_p_scam_abg, color = cohort)) +
  ggplot2::geom_point() +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = p_scam_bgb_lower, ymax = p_scam_bgb_upper),width = 3,linewidth = 0.5)+
  ggplot2::facet_wrap(~ env_treatment.x)+
  ggplot2::theme_bw() +
  ggplot2::scale_color_manual(values = color_map,
                              labels = c("blackwater\nancestral",
                                         "blackwater\nmodern",
                                         "corn\nancestral",
                                         "corn\nmodern")) +
  ggplot2::ylim(0,1) +
  ggplot2::xlim(0,1) +
  ggplot2::labs(x = "", y = "Estimated Bgb p.scam MID") +
  ggplot2::geom_abline(intercept = 0,slope = 1, linetype = "dashed", color = "red") +
  ggplot2::theme(legend.position = "none")


BTM <- bgb_biomass_layer[["BTM"]] # estimated
traits_BTM <- comp_traits[comp_traits$pot_no %in% BTM$pot_no,] 
traits_BTM$env_treatment <- paste0(ifelse(traits_BTM$salinity == 0, "Rhode River", "GCREW"),
                                   "/",
                                   ifelse(traits_BTM$elev_bi == 0, "high inundation", "low inundation"))
merged_BTM <- merge(BTM, traits_BTM, by = "pot_no", all.x = TRUE)

average_p_scam_BTM <- merged_BTM |>
  group_by(env_treatment.x, gt, seed_year.x, cohort) |>
  summarize(
    mean_p_scam_abg = mean(p.scam.abg, na.rm = TRUE),
    mean_p_scam_bgb = mean(p.scam, na.rm = TRUE),
    p_scam_bgb_lower = mean(p.scam_lower, na.rm = TRUE),
    p_scam_bgb_upper = mean(p.scam_upper, na.rm = TRUE),
    .groups = "keep")

ggBTM <- ggplot2::ggplot(data = average_p_scam_BTM, ggplot2::aes(y = mean_p_scam_bgb, x = mean_p_scam_abg, color = cohort)) +
  ggplot2::geom_point() +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = p_scam_bgb_lower, ymax = p_scam_bgb_upper),width = 3,linewidth = 0.5)+
  ggplot2::facet_wrap(~ env_treatment.x)+
  ggplot2::theme_bw() +
  ggplot2::scale_color_manual(values = color_map,
                              labels = c("blackwater\nancestral",
                                         "blackwater\nmodern",
                                         "corn\nancestral",
                                         "corn\nmodern")) +
  ggplot2::ylim(0,1) +
  ggplot2::xlim(0,1) +
  ggplot2::labs(x = "Measured Agb p.scam", y = "Estimated Bgb p.scam BTM") +
  ggplot2::geom_abline(intercept = 0,slope = 1, linetype = "dashed", color = "red") +
  ggplot2::theme(legend.position = "none")

gridExtra::grid.arrange(ggTOP, ggMID, ggBTM, nrow=3)  
