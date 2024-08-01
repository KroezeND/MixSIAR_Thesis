(rm(list = ls()))

source("src/append_functions.R")


# Load in compiled weight values
bgb_biomass <- read.csv("data/bgb_biomass.csv") 

# Load in Clean and Formatted Data
mixture_TOP <- read.csv("data/full_TOP.csv") # training dataset 
mixture_MID <- read.csv("data/full_MID.csv") # training dataset 
mixture_BTM <- read.csv("data/full_BTM.csv") # training dataset 

mixture <- list(mixture_TOP,mixture_MID,mixture_BTM)

# Load jags.1 output from model run script
jags.EnvGT <- readRDS("output/Full/jags.EnvGT.rds")

bgb_biomass_layer <- list()
for(mod in 1:length(jags.EnvGT)){
  # Directly access MCMC draws for SCAM mixing proportion for each individual 
  # assumes SCAM is first mixing element
  p.scam <- jags.EnvGT[[mod]]$BUGSoutput$sims.list$p.ind[,,1]
  p.sppa <- jags.EnvGT[[mod]]$BUGSoutput$sims.list$p.ind[,,2]
  
  # Add mean mixing proportion of each individual to mixture dataset
  mixture[[mod]]$p.scam <- colMeans(p.scam)
  # Calculate lower and upper bounds (assuming 95% CrI)
  mixture[[mod]]$p.scam_lower <- apply(p.scam, 2, FUN = quantile, probs = 0.025)
  mixture[[mod]]$p.scam_upper <- apply(p.scam, 2, FUN = quantile, probs = 0.975)
  
  mixture[[mod]]$p.sppa <- colMeans(p.sppa)
  
  mixture[[mod]]$p.sppa_lower <- apply(p.sppa, 2, FUN = quantile, probs = 0.025)
  mixture[[mod]]$p.sppa_upper <- apply(p.sppa, 2, FUN = quantile, probs = 0.975)
  
  bgb_biomass_mod <- bgb_biomass[bgb_biomass$pot_no %in% mixture[[mod]]$pot_no,] # all weights that have isotopic measurements
  # Add depth column to bgb_biomass using lapply
  bgb_biomass_mod$Depth <- as.character(lapply(bgb_biomass_mod$segment_bottom, assign_depth))
  # Ensure it is only the active mixture depth layer
  bgb_biomass_mod <- bgb_biomass_mod[bgb_biomass_mod$Depth %in% mixture[[mod]]$Depth,]
  
  # Use custom function for joining mixture and bgb_biomass 
  bgb_biomass_mod <- join_and_assign_values(bgb_biomass_mod,mixture[[mod]])
  
  # Remove extraneous organizational rows
  bgb_biomass_mod <- dplyr::select(bgb_biomass_mod, -all_of(c("date_weighed","box_no","date_washed",
                                                          "notes","midpoint")))
  
  # Since some of the pots collapsed when processing the cores, it was not possible to get an accurate
  # length. In order to account for the effect of volume on root density, we are using the average length
  # of the end segment in order to include those collapsed pots, segment_bottom_adjusted
  # Calculate the weight density in a cylinder with radius 5cm (0.05m) and height = segment_bottom - segment_top
  bgb_biomass_mod$weight_density <-  bgb_biomass_mod$weight_g / 
    (pi * 5^2 * ((bgb_biomass_mod$segment_bottom_adjusted - bgb_biomass_mod$segment_top)))
  
  # Calculates prop_scam and prop_sppa
  bgb_biomass_mod <- join_and_calculate_props(bgb_biomass_mod,mixture[[mod]])
  
  bgb_biomass_layer[[mod]] <- bgb_biomass_mod
}
names(bgb_biomass_layer) <- c("TOP","MID","BTM")
names(mixture) <- c("TOP","MID","BTM")


colors <- c("#2d6a4f","#74c69d","#5a189a","#c77dff")
# Plot the mixing proportion of SCAM
mixture[["TOP"]] |>
  ggplot2::ggplot(ggplot2::aes(x = seed_year, y = p.scam, color = cohort)) +
  ggplot2::geom_point() +  # Plot points for each data point
  ggplot2::geom_errorbar(ggplot2::aes(ymin = p.scam_lower, ymax = p.scam_upper), width = 3,linewidth = 0.55) + # Add ribbon for CrI
  ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +  # Add horizontal line
  # geom_line(aes(y = p.scam), stat = "smooth", method = "loess") +  # Add smoothed line
  ggplot2::scale_color_manual(values = colors,labels = c("blackwater\nancestral","blackwater\nmodern",
                                                "corn\nancestral","corn\nmodern"))+
  ggplot2::facet_wrap(~ env_treatment) +  # Separate panels for each treatment
  ggplot2::ylim(0,1) +
  ggplot2::labs(title = "p.scam in TOP Layer by Year and Treatment (95% CrI)", x = "Year", y = "p.scam") +
  ggplot2::theme_bw() + # Set a black and white theme (optional)
  ggplot2::guides(color = ggplot2::guide_legend(title = "Age Cohort",reverse=TRUE ))

mixture[["MID"]] |>
  ggplot2::ggplot(ggplot2::aes(x = seed_year, y = p.scam, color = cohort)) +
  ggplot2::geom_point() +  # Plot points for each data point
  ggplot2::geom_errorbar(ggplot2::aes(ymin = p.scam_lower, ymax = p.scam_upper), width = 3,linewidth = 0.55) + # Add ribbon for CrI
  ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +  # Add horizontal line
  # geom_line(aes(y = p.scam), stat = "smooth", method = "loess") +  # Add smoothed line
  ggplot2::scale_color_manual(values = colors,labels = c("blackwater\nancestral","blackwater\nmodern",
                                                "corn\nancestral","corn\nmodern"))+
  ggplot2::facet_wrap(~ env_treatment) +  # Separate panels for each treatment
  ggplot2::ylim(0,1) +
  ggplot2::labs(title = "p.scam in MID Layer by Year and Treatment (95% CrI)", x = "Year", y = "p.scam") +
  ggplot2::theme_bw() + # Set a black and white theme (optional)
  ggplot2::guides(color = ggplot2::guide_legend(title = "Age Cohort" ,reverse=TRUE))

mixture[["BTM"]] |>
  ggplot2::ggplot(ggplot2::aes(x = seed_year, y = p.scam, color = cohort)) +
  ggplot2::geom_point() +  # Plot points for each data point
  ggplot2::geom_errorbar(ggplot2::aes(ymin = p.scam_lower, ymax = p.scam_upper), width = 3,linewidth = 0.55) + # Add ribbon for CrI
  ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +  # Add horizontal line
  # geom_line(aes(y = p.scam), stat = "smooth", method = "loess") +  # Add smoothed line
  ggplot2::scale_color_manual(values = colors,labels = c("blackwater\nancestral","blackwater\nmodern",
                                                "corn\nancestral","corn\nmodern"))+
  ggplot2::facet_wrap(~ env_treatment) +  # Separate panels for each treatment
  ggplot2::ylim(0,1) +
  ggplot2::labs(title = "p.scam in BTM Layer by Year and Treatment (95% CrI)", x = "Year", y = "p.scam") +
  ggplot2::theme_bw() + # Set a black and white theme (optional)
  ggplot2::guides(color = ggplot2::guide_legend(title = "Age Cohort",reverse=TRUE))


saveRDS(bgb_biomass_layer,"output/bgb_biomass_layer.rds")
saveRDS(mixture,"output/mixture_list.rds")
