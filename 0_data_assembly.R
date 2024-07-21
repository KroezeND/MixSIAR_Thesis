# Author: N. Kroeze
# Description: This script cleans the raw biomass and isotope data and formats each 

# Inputs:
# "Belowground_Biomass.csv" Biomass census spreadsheet from Kleiner 2022
# "BlueGenes2023_Belowground_Biomass" - BG2023.csv Newly collected biomass spreadsheet 
# "Compiled_Traits.csv" Trait census data from Kleiner 2022
# "220729_JPM.xls" - CN sum.csv Smithsonian lab isotope data
# "220801_JPM.xls" - CN sum.csv Smithsonian lab isotope data
# "220519_JPM.xls" - CN sum.csv Smithsonian lab isotope data
# "CEST_Sample_Runs.csv" CEST lab isotope data - all runs, not including replicate testing
# "Trial_Carbon_Data_unedited.csv" CEST lab isotope data - replicate testing
# Corr_Trial_COMP.csv CEST lab isotope data - replicate testing
# Milled_replicates.csv CEST lab isotope data - replicate testing

# Outputs:
# "marsh_source_all.csv" raw source isotopes formatted to MixSIAR specifications, 151 observations of 2 variables
# "Compiled_Traits_appended.csv" traits with added weights and root-shoot-ratios, 672 observations of 35 variables
# "training_TOP.csv" mixture distribution for training MixSIAR model, 99 observations of 17 variables
# "training_MID.csv" mixture distribution for training MixSIAR model, 99 observations of 17 variables
# "training_BTM.csv" mixture distribution for training MixSIAR model, 99 observations of 17 variables
# "bgb_biomass.csv" all weights across each layer, 1692 observations of 11 variables
# "holdout_DONOTTOUCH.csv" holdout data, 90 observations of 15 variables

library(tidyverse)
library(here)

(rm(list = ls()))

source("src/data_process.R")
set.seed(1)
# Load and Append Trait data with full belowground biomass calculation from 2023 efforts
# Trait Data from HSK 2021 BlueGenes experiment, describes the environmental conditions and genotype 
traits <- read.csv("data/Biomass/Compiled_Traits.csv")

# Belowground Biomass Weights
bgb_ND <- read.csv("data/Biomass/BlueGenes2023_Belowground_Biomass - BG2023.csv",header = T) 
bgb_SI <- read.csv("data/Biomass/Belowground_Biomass.csv")

# Update bgb_ND with total weight_g per pot_no
bgb_ND_updated <- bgb_ND %>%
  dplyr::group_by(pot_no) %>%
  dplyr::summarise(total_weight_g = sum(weight_g))

# Merge bgb_ND_updated with traits dataframe, updates with new weights
traits_updated <- merge(traits, bgb_ND_updated[, c("pot_no", "total_weight_g")], by = "pot_no", all.x = TRUE) %>%
  dplyr::mutate(tot_bgb = dplyr::coalesce(tot_bgb, total_weight_g)) %>% # Combine values using coalesce
  dplyr::select(-tidyselect::one_of(c("total_weight_g")))

#Calculate rsr in traits
rsr <- traits_updated$tot_bgb/traits_updated$tot_agb
traits_updated$rsr <- rsr

bgb_biomass <- rbind(bgb_ND,bgb_SI)
# Assumed end considers the average length of each core and imposes this due to core collapse
assumed_end <- mean(as.numeric(bgb_biomass$segment_bottom[bgb_biomass$segment_top==20]),na.rm = TRUE)
bgb_biomass <- bgb_biomass %>%
  dplyr::mutate(segment_bottom_adjusted = ifelse(is.na(segment_bottom) | segment_bottom == "end", 
                                                 assumed_end,segment_bottom),
                midpoint = round((segment_top + as.numeric(segment_bottom_adjusted)) / 2,3))


# Load Raw Isotope Results
# SI Data for Low Sal/Low Elevation
ISO_SI <- do.call(rbind,lapply(c("data/Isotope/SI/220519_JPM.xls - CN sum.csv",
                                 "data/Isotope/SI/220729_JPM.xls - CN sum.csv",
                                 "data/Isotope/SI/220801_JPM.xls - CN sum.csv"),
                               function(file) process_SI_data(read.csv(file, header = TRUE))))
# process_SI_data is a custom function created to manually match the Smithsonian lab
# outputs to the CEST lab outputs

# Method Testing at CEST: SI v. ND, coarse v. fine v. mill, accuracy, etc
ISO_trial <- do.call(rbind, lapply(c("data/Isotope/CEST/Trial_Carbon_Data_unedited.csv",
                                     "data/Isotope/CEST/Corr_Trial_COMP.csv",
                                     "data/Isotope/CEST/Milled_replicates.csv"),
                                   function(file) process_CEST_data(read.csv(file, header = TRUE))))
# ISO_trial represents the test runs that assessed the repeatability of results 
# across grinding methods, across labs/instruments, and across replicates

ISO_ND <- do.call(rbind, lapply(c("data/Isotope/CEST/CEST_Sample_Runs.csv"),
                                function(file) process_CEST_data(read.csv(file, header = TRUE))))
# ISO_ND are the raw isotope from CEST, processed in the same way as the SI data

ISO_ALL <- rbind(ISO_SI, ISO_trial, ISO_ND) # bind all measurements into one dataframe
ISO_ALL <- merge(ISO_ALL, traits[, c("pot_no", "salinity", "elev_bi","gt",
                                     "cohort","age_cohort","provenance","seed_year")], 
                 by = "pot_no", all.x = TRUE) # add environmental and origin conditions to dataframe
# Create source dataframe (only SCAM and SPPA)
ISO_source <- subset(ISO_ALL, !grepl("COMP", ISO_ALL$Species))

# Filter out ecologically improbable data since the end members should be well differentiated
ISO_source <- ISO_source %>%
  dplyr::filter(!(Species == "SCAM" & d13C > -20) & !(Species == "SPPA" & d13C < -25))

# Index for just source and d13C
ISO_source <- ISO_source[c("Species","d13C")]
ISO_source <- dplyr::rename(ISO_source, Source = Species) #rename to match expected column names

ISO_source$Source <- tolower(ISO_source$Source)

# Create total competition dataframe
ISO_COMP <- subset(ISO_ALL, grepl("COMP", ISO_ALL$Species))

# Use traits to set up holdout data, choosing holdout by pot_no not individual
traits_holdout <- traits[traits$pot_no %in% ISO_COMP$pot_no,]

# Define unique scenarios of interest; salinity, elevation, and age cohort
traits_holdout$scenario <- interaction(salinity = traits_holdout$salinity,
                                       elev_bi = traits_holdout$elev_bi,
                                       age_cohort = traits_holdout$age_cohort)

# Minimum sample size for holdout (consider adjusting)
min_holdout_size <- 2# Adjust based on minimum desired holdout size

# Split data using stratified sampling
holdout_data <- traits_holdout %>%
  dplyr::group_by(scenario) %>%
  dplyr::do({
    data <- . # Reference to the group data
    holdout_sample(data, target_size = 0.23 * nrow(data))
  })


holdout_ids <- unique(holdout_data$pot_no)

holdout <- ISO_COMP[ISO_COMP$pot_no %in% holdout_ids,] # actual hold out data to be saved

# Training data is everything not in holdout
training_data <- traits_holdout[!(traits_holdout$pot_no %in% holdout_ids), ]

# Check the distribution of scenarios in both sets 
table(traits_holdout$scenario) %>% prop.table() # proportion of each scenario for all COMP
table(holdout_data$scenario) %>% prop.table()   # proportion of holdout data 
table(training_data$scenario) %>% prop.table()  # proportion of training dataset

# Get all depth layers for training data
training <- ISO_COMP[ISO_COMP$pot_no %in% training_data$pot_no,]

# Add env_treatment to training dataset
training$env_treatment <- paste0(ifelse(training$salinity == 0, "Rhode River", "GCREW"),
                                 "/",
                                 ifelse(training$elev_bi == 0, "high inundation", "low inundation"))

# Define depth conversion dictionary
depth_conversion <- c(TOP = 0, MID = 10, BTM = 20)

# Append training with midpoint value from bgb_biomass
for (i in 1:nrow(training)) {
  current_pot_no <- training$pot_no[i]
  current_depth <- training$Depth[i]
  converted_depth <- depth_conversion[[current_depth]]  # Convert depth from training
  
  # Matching rows based on pot_no and converted depth
  matching_rows <- bgb_biomass[bgb_biomass$pot_no == current_pot_no & bgb_biomass$segment_top == converted_depth,]
  
  # Check if a matching row is found
  if (nrow(matching_rows) > 0) {
    # Extract midpoint from the matching row
    matched_midpoint <- matching_rows$midpoint
    
    # Append the new row to training
    training$midpoint[training$pot_no == current_pot_no & training$Depth==current_depth] <- matched_midpoint
  }
}

# Divide by depth layer (mixture distribution)
training_TOP <- training[training$Depth=="TOP",]
training_MID <- training[training$Depth=="MID",]
training_BTM <- training[training$Depth=="BTM",]

# # Average across replicates to utilize all available information
new_training_TOP <- average_pot_no(training_TOP)
new_training_MID <- average_pot_no(training_MID)
new_training_BTM <- average_pot_no(training_BTM)


# Save modified objects
write.csv(ISO_source, "data/marsh_source_all.csv",row.names=FALSE)
write.csv(traits_updated,"data/Compiled_Traits_appended.csv")
write.csv(new_training_TOP, "data/training_TOP.csv",row.names = FALSE)
write.csv(new_training_MID, "data/training_MID.csv",row.names = FALSE)
write.csv(new_training_BTM, "data/training_BTM.csv",row.names = FALSE)
write_csv(bgb_biomass, "data/bgb_biomass.csv")
write.csv(holdout, "data/Holdout/holdout_DONOTTOUCH.csv",row.names = T)

