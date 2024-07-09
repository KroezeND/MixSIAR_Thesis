library(tidyverse)
library(here)

(rm(list = ls()))

source("src/append_functions.R")

jags.Depth <- readRDS("output/Environment_Depth Model/jags.Depth.rds")

# Load in Clean and Formatted Data
mixture <- read.csv("data/training.csv") # training dataset 

# Directly access MCMC draws for SCAM mixing proportion for each individual 
# assumes SCAM is first mixing element
p.scam <- jags.Depth$BUGSoutput$sims.list$p.ind[,,1]
p.sppa <- jags.Depth$BUGSoutput$sims.list$p.ind[,,2]

# Add mean mixing proportion of each individual to mixture dataset
mixture$p.scam <- colMeans(p.scam)

# Calculate lower and upper bounds (assuming 95% CrI)
mixture$p.scam_lower <- apply(p.scam, 2, FUN = quantile, probs = 0.025)
mixture$p.scam_upper <- apply(p.scam, 2, FUN = quantile, probs = 0.975)

mixture$p.sppa <- colMeans(p.sppa)

mixture$p.sppa_lower <- apply(p.sppa, 2, FUN = quantile, probs = 0.025)
mixture$p.sppa_upper <- apply(p.sppa, 2, FUN = quantile, probs = 0.975)

# Load in compiled weight values
bgb_biomass <- read.csv("data/bgb_biomass.csv") 
bgb_biomass <- bgb_biomass[bgb_biomass$pot_no %in% mixture$pot_no,] # all weights that have isotopic measurements

# Add depth column to bgb_biomass using lapply
bgb_biomass$Depth <- as.character(lapply(bgb_biomass$segment_bottom, assign_depth))

# Use custom function for joining mixture and bgb_biomass 
bgb_biomass <- join_and_assign_values(bgb_biomass,mixture)

# Remove extraneous organizational rows
bgb_biomass <- dplyr::select(bgb_biomass, all_of(c("pot_no","segment_top","segment_bottom",
                                                   "env_treatment","seed_year", "weight_g","Depth",
                                                   "p.scam", "p.scam_lower", "p.scam_upper",
                                                   "p.sppa","p.sppa_lower","p.sppa_upper")))
# Since some of the pots collapsed when processing the cores, it was not possible to get an accurate
# length. In order to account for the effect of volume on root density, we are using the average length
# of the end segment in order to include those collapsed pots
assumed_end <- mean(as.numeric(bgb_biomass$segment_bottom[bgb_biomass$segment_top==20]),na.rm = TRUE)
# Calculate the weight density in a cylinder with radius 5cm (0.05m) and height = segment_bottom - segment_top
bgb_biomass$weight_density <- ifelse(bgb_biomass$segment_bottom == "end", 
                                     bgb_biomass$weight_g / (pi * 5^2 * ((assumed_end - as.numeric(bgb_biomass$segment_top)))),
                                     bgb_biomass$weight_g / (pi * 5^2 * ((as.numeric(bgb_biomass$segment_bottom) - bgb_biomass$segment_top))))

# Calculates prop_scam and prop_sppa
bgb_biomass <- join_and_calculate_props(bgb_biomass,mixture)


# Plot the mixing proportion of SCAM
ggplot(mixture, aes(x = seed_year, y = p.scam, color = env_treatment, shape = Depth)) +
  geom_point() +  # Plot points for each data point
  geom_errorbar(aes(ymin = p.scam_lower, ymax = p.scam_upper), linetype = 2, alpha = 0.15) + # Add ribbon for CrI
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +  # Add horizontal line
  # geom_line(aes(y = p.scam), stat = "smooth", method = "loess") +  # Add smoothed line
  # facet_wrap(~ env_treatment) +  # Separate panels for each treatment
  ylim(0,1) +
  labs(title = "p.scam by Year and Treatment (95% CrI)", x = "Year", y = "p.scam") +
  theme_bw()  # Set a black and white theme (optional)

unique_years <- c("pre-1950","1951-2000","post-2000")
colors <- c("#c51b7d", "#7fbc41","gray30")
color_map <- setNames(colors, unique_years)

bgb_biomass$age_group <- ifelse(bgb_biomass$seed_year < 1950, "pre-1950", 
                                ifelse(bgb_biomass$seed_year >= 1951 & bgb_biomass$seed_year <= 2000, "1951-2000", "post-2000"))


bgb_biomass %>%
  mutate(age_group = fct_relevel(age_group,"pre-1950","1951-2000","post-2000")) %>%
  ggplot(aes(y=Depth, x=prop_scam, fill = age_group)) + 
  geom_violin(position=position_dodge(0.95),scale = 'width',trim = T) +  # Dodge violin plots slightly
  geom_boxplot(position=position_dodge(0.95),width = 0.25) +
  scale_fill_manual(values = color_map) +  # Set color based on year
  labs(y = "Depth", x = "Weight Density SCAM (g/cm3)") + # Custom labels
  theme_bw()+
  geom_hline(yintercept = 1.5:2.5,color = "red",linetype = "dashed")+
  guides(fill = guide_legend(reverse=TRUE))


bgb_biomass[bgb_biomass$Depth=="BTM",] %>%
  mutate(age_group = fct_relevel(age_group,"pre-1950","1951-2000","post-2000")) %>%
  ggplot(aes(y=Depth, x=prop_scam, fill = age_group)) + 
  geom_violin(position=position_dodge(0.95),scale = 'width',trim = T) +  # Dodge violin plots slightly
  geom_boxplot(position=position_dodge(0.95),width = 0.25) +
  scale_fill_manual(values = color_map) +  # Set color based on year
  labs(y = "Depth", x = "Weight Density SCAM (g/cm3)") + # Custom labels
  theme_bw()+
  guides(fill = guide_legend(reverse=TRUE)) 
 