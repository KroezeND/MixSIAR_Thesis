# Author: N. Kroeze
# Description: 

(rm(list = ls()))

source("src/append_functions.R")

bgb_biomass_layer <- readRDS("output/bgb_biomass_layer.rds") #contains the mixing proportions and weight density
traits <- read.csv("data/Compiled_Traits_appended.csv") # census trait data
mixture <- readRDS("output/mixture_list.rds") # list containing the three mixture distributions used in each model
bgb_biomass <- read.csv("data/bgb_biomass.csv") # all belowground biomass weights spreadsheet

# Add environmental treatment to trait census data
traits$env_treatment <- paste0(ifelse(traits$salinity == 0, "low_sal", "high_sal"),
                                         "/",
                                         ifelse(traits$elev_bi == 0, "low_elev", "high_elev"))
traits <- traits[-which(traits$species == "empty"), ]  # Modify the original data frame

# Filtered data (assuming missing values in tot_bgb are NA)
traits_filtered <- traits |>
  dplyr::filter(!is.na(tot_bgb))

bgb_biomass$env_treatment <- NA  # Create a new column for env_treatment (filled with NA initially)
for (i in 1:nrow(bgb_biomass)) {
    pot_no <- bgb_biomass[i, "pot_no"]
    matching_env <- traits[traits$pot_no == pot_no, "env_treatment",]
    if (length(matching_env) == 1) {
      bgb_biomass$env_treatment[i] <- matching_env  # Assign env_treatment if there's a unique match
    }
}

bgb_biomass |>
  dplyr::filter(segment_top == 0) |>
  dplyr::mutate(env_treatment = dplyr::if_else(env_treatment == 'high_sal/high_elev', 'GCREW/\nlow inundation', env_treatment),
                env_treatment = dplyr::if_else(env_treatment == 'high_sal/low_elev', 'GCREW/\nhigh inundation', env_treatment),
                env_treatment = dplyr::if_else(env_treatment == 'low_sal/high_elev', 'Rhode River/\nlow inundation', env_treatment),
                env_treatment = dplyr::if_else(env_treatment == 'low_sal/low_elev', 'Rhode River/\nhigh inundation', env_treatment)) |>
  ggplot2::ggplot(ggplot2::aes(y = weight_g, x = factor(env_treatment), fill = factor(species))) +
  ggplot2::geom_violin(position=ggplot2::position_dodge(0.90),scale = 'width',trim = T) + 
  ggplot2::geom_boxplot(position=ggplot2::position_dodge(0.9),width = 0.25) +
  ggplot2::theme_bw() +
  ggplot2::ylim(0, 17.5) +
  ggplot2::ggtitle("Belowground Biomass in TOP Sediment Layer by Environment") +
  ggplot2::labs(x = "Environmental Conditions", y = "Belowground Biomass (g)", color = "Competition") + # Custom labels
  ggplot2::scale_fill_manual(values = c("#EE964B","#19647E", "#ED254E"),labels = c("COMP", "SCAM","SPPA"))+ # Set custom colors
  ggplot2::guides(fill = ggplot2::guide_legend(title = "Species"))


unique_cohort <- sort(unique(mixture$TOP$cohort))
colors <- c("#2d6a4f","#74c69d","#5a189a","#c77dff")
color_map <- setNames(colors, unique_cohort)


for(i in names(bgb_biomass_layer)){
  bgb_biomass_layer[[i]] <- bgb_biomass_layer[[i]] |>
    dplyr::left_join(mixture[[i]], by = "pot_no")
}



bgb_biomass_layer[["TOP"]] |>
  dplyr::mutate(env_treatment.x = sub(pattern = '/', replacement = '/\n', x = env_treatment.x)) |>
  ggplot2::ggplot(ggplot2::aes(y=env_treatment.x, x=prop_scam, fill = cohort)) + 
  ggplot2::geom_violin(position=ggplot2::position_dodge(0.95),scale = 'width',trim = T,drop = F) +  # Dodge violin plots slightly
  #ggplot2::facet_wrap(~ env_treatment.x) +
  ggplot2::geom_boxplot(position=ggplot2::position_dodge(0.95),width = 0.25) +
  ggplot2::scale_fill_manual(values = color_map,labels = c("blackwater\nancestral","blackwater\nmodern",
                                                  "corn\nancestral","corn\nmodern")) +
  ggplot2::labs(y = "Environmental Treatment", x = expression("Dry Weight Density SCAM under Competition (g/cm"^3*")")) + 
  ggplot2::theme_bw()+
  ggplot2::xlim(0,0.022) +
  ggplot2::guides(fill = ggplot2::guide_legend(title ="Cohort",reverse=TRUE))


bgb_biomass_layer[["MID"]] |>
  dplyr::mutate(env_treatment.x = sub(pattern = '/', replacement = '/\n', x = env_treatment.x)) |>
  ggplot2::ggplot(ggplot2::aes(y=env_treatment.x, x=prop_scam, fill = cohort)) + 
  ggplot2::geom_violin(position=ggplot2::position_dodge(0.95),scale = 'width',trim = T,drop = F) +  # Dodge violin plots slightly
  ggplot2::geom_boxplot(position=ggplot2::position_dodge(0.95),width = 0.25) +
  #ggplot2::facet_wrap(~ env_treatment.x) +
  ggplot2::scale_fill_manual(values = color_map,labels = c("Blackwater Ancestral","Blackwater Modern",
                                                  "Corn Ancestral","Corn Modern")) +
  ggplot2::labs(y = "Environmental Treatment", x = expression("Dry Weight Density SCAM under Competition (g/cm"^3*")")) + 
  ggplot2::theme_bw()+
  ggplot2::xlim(0,0.01) +
  ggplot2::guides(fill = ggplot2::guide_legend(title ="Cohort",reverse=TRUE))

bgb_biomass_layer[["BTM"]] |>
  dplyr::mutate(env_treatment.x = sub(pattern = '/', replacement = '/\n', x = env_treatment.x)) |>
  ggplot2::ggplot(ggplot2::aes(y=env_treatment.x, x=prop_scam, fill = cohort)) + 
  ggplot2::geom_violin(position=ggplot2::position_dodge(0.95),scale = 'width',trim = T,drop = F) +  # Dodge violin plots slightly
  ggplot2::geom_boxplot(position=ggplot2::position_dodge(0.95),width = 0.25) +
  #ggplot2::facet_wrap(~ env_treatment.x) +
  ggplot2::scale_fill_manual(values = color_map,labels = c("Blackwater Ancestral","Blackwater Modern",
                                                  "Corn Ancestral","Corn Modern")) +  # Set color based on year
  ggplot2::labs(y = "Environmental Treatment", x = expression("Dry Weight Density SCAM under Competition (g/cm"^3*")")) + 
  ggplot2::theme_bw()+
  ggplot2::xlim(0,0.0019) +
  ggplot2::guides(fill = ggplot2::guide_legend(title ="Cohort",reverse=TRUE))

# Comparison to SCAM Monocultures
bgb <- bgb_biomass[bgb_biomass$species=="scam",]
assumed_end <- mean(as.numeric(bgb$segment_bottom[bgb$segment_top==20]),na.rm = TRUE)
# Calculate the weight density in a cylinder with radius 5cm (0.05m) and height = segment_bottom - segment_top
bgb$weight_density <- ifelse(bgb$segment_bottom == "end", 
                             bgb$weight_g / (pi * 5^2 * ((assumed_end - as.numeric(bgb$segment_top)))),
                             bgb$weight_g / (pi * 5^2 * ((as.numeric(bgb$segment_bottom) - bgb$segment_top))))

bgb <- merge(bgb, traits[, c("pot_no", "salinity", "elev_bi","gt",
                             "cohort","age_cohort","provenance","seed_year")], 
             by = "pot_no", all.x = TRUE) # add environmental and origin conditions to da
bgb$Depth <- as.character(lapply(bgb$segment_bottom, assign_depth))
bgb$env_treatment <- paste0(ifelse(bgb$salinity == 0, "Rhode River", "GCREW"),
                            "/",
                            ifelse(bgb$elev_bi == 0, "high inundation", "low inundation"))


bgb_TOP <- bgb[bgb$Depth == "TOP",]
bgb_MID <- bgb[bgb$Depth == "MID",]
bgb_BTM <- bgb[bgb$Depth == "BTM",]

bgb_TOP |>
  dplyr::mutate(env_treatment = sub(pattern = '/', replacement = '/\n', x = env_treatment)) |>
  ggplot2::ggplot(ggplot2::aes(y=env_treatment, x=weight_density, fill = cohort)) + 
  ggplot2::geom_violin(position=ggplot2::position_dodge(0.95),scale = 'width',trim = TRUE) +  # Dodge violin plots slightly
  ggplot2::geom_boxplot(position=ggplot2::position_dodge(0.95),width = 0.25) +
  #ggplot2::facet_wrap(~ env_treatment) +
  ggplot2::scale_fill_manual(values = color_map,labels = c("blackwater\nancestral","blackwater\nmodern",
                                                  "corn\nancestral","corn\nmodern")) +
  ggplot2::labs(y = "Environmental Treatment", x = expression(Weight ~ Density~SCAM~Monoculture~(g/cm^{3}))) + # Custom labels
  ggplot2::theme_bw()+
  ggplot2::xlim(0,0.022) +
  ggplot2::guides(fill = ggplot2::guide_legend(title = "Cohort",reverse=TRUE)) 

bgb_MID |>
  dplyr::mutate(env_treatment = sub(pattern = '/', replacement = '/\n', x = env_treatment)) |>
  ggplot2::ggplot(ggplot2::aes(y=env_treatment, x=weight_density, fill = cohort)) + 
  ggplot2::geom_violin(position=ggplot2::position_dodge(0.95),scale = 'width',trim = TRUE) +  # Dodge violin plots slightly
  ggplot2::geom_boxplot(position=ggplot2::position_dodge(0.95),width = 0.25) +
  #ggplot2::facet_wrap(~ env_treatment) +
  ggplot2::scale_fill_manual(values = color_map,labels = c("blackwater\nancestral","blackwater\nmodern",
                                                  "corn\nancestral","corn\nmodern")) +
  ggplot2::labs(y = "Environmental Treatment", x = expression(Weight ~ Density~SCAM~Monoculture~(g/cm^{3}))) + # Custom labels
  ggplot2::theme_bw()+
  ggplot2::xlim(0,0.01) +
  ggplot2::guides(fill = ggplot2::guide_legend(title = "Cohort",reverse=TRUE)) 

bgb_BTM |>
  dplyr::mutate(env_treatment = sub(pattern = '/', replacement = '/\n', x = env_treatment)) |>
  ggplot2::ggplot(ggplot2::aes(y=env_treatment, x=weight_density, fill = cohort)) + 
  ggplot2::geom_violin(position=ggplot2::position_dodge(0.95),scale = 'width',trim = TRUE) +  # Dodge violin plots slightly
  ggplot2::geom_boxplot(position=ggplot2::position_dodge(0.95),width = 0.25) +
  #ggplot2::facet_wrap(~ env_treatment) +
  ggplot2::scale_fill_manual(values = color_map,labels = c("blackwater\nancestral","blackwater\nmodern",
                                                  "corn\nancestral","corn\nmodern")) +
  ggplot2::labs(y = "Environmental Treatment", x = expression(Weight ~ Density~SCAM~Monoculture~(g/cm^{3}))) + # Custom labels
  ggplot2::theme_bw()+
  ggplot2::xlim(0,0.0019) +
  ggplot2::guides(fill = ggplot2::guide_legend(title = "Cohort",reverse=TRUE)) 
