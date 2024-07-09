# Define a function to assign depth based on segment_bottom
assign_depth <- function(segment_bottom) {
  if (segment_bottom == 10) {
    return("TOP")
  } else if (segment_bottom == 20) {
    return("MID")
  } else {
    return("BTM")
  }
}


# Function to join bgb_biomass with mixture and assign matching values
join_and_assign_values <- function(bgb_biomass, mixture) {
  
  # Iterate through bgb_biomass rows
  for (i in 1:nrow(bgb_biomass)) {
    current_row <- bgb_biomass[i, ]  # Extract current row from bgb_biomass
    pot_no <- current_row$pot_no
    depth <- as.character(current_row$Depth)
    
    # Find matching row in mixture
    matching_row <- mixture[mixture$pot_no == pot_no & mixture$Depth == depth, ]
    
    # Assign values if a match is found
    if (nrow(matching_row) > 0) {
      bgb_biomass$p.scam[i] <- unique(matching_row$p.scam)
      bgb_biomass$p.scam_lower[i] <- unique(matching_row$p.scam_lower)
      bgb_biomass$p.scam_upper[i] <- unique(matching_row$p.scam_upper)
      bgb_biomass$env_treatment[i] <- unique(matching_row$env_treatment)
      bgb_biomass$seed_year[i] <- unique(matching_row$seed_year)
      bgb_biomass$p.sppa[i] <- unique(matching_row$p.sppa)
      bgb_biomass$p.sppa_lower[i] <- unique(matching_row$p.sppa_lower)
      bgb_biomass$p.sppa_upper[i] <- unique(matching_row$p.sppa_upper)
    }
  }
  
  # Return the updated bgb_biomass dataframe
  return(bgb_biomass)
}


# Function to join bgb_biomass with mixture and calculate prop_scam & prop_sppa
join_and_calculate_props <- function(bgb_biomass, mixture) {
  # Loop through bgb_biomass rows
  for (i in 1:nrow(bgb_biomass)) {
    # Extract current row values
    current_pot_no <- bgb_biomass$pot_no[i]
    current_depth <- bgb_biomass$Depth[i]
    
    # Find matching row with depth (if possible)
    matching_row_depth <- mixture[mixture$pot_no == current_pot_no & mixture$Depth == current_depth, ]
    
    # Check if depth match is found
    has_depth_match <- nrow(matching_row_depth) > 0
    
    # Initialize variables for prop_scam calculations
    p.scam_value <- NA
    p.scam_lower <- NA
    p.scam_upper <- NA
    
    # Initialize variables for prop_sppa calculations
    p.sppa_value <- NA
    p.sppa_lower <- NA
    p.sppa_upper <- NA
    
    # Extract values from depth match (if found)
    if (has_depth_match) {
      p.scam_value <- matching_row_depth$p.scam[1]
      p.scam_lower <- ifelse(is.na(matching_row_depth$p.scam_lower[1]), NA, matching_row_depth$p.scam_lower[1])
      p.scam_upper <- ifelse(is.na(matching_row_depth$p.scam_upper[1]), NA, matching_row_depth$p.scam_upper[1])
      
      p.sppa_value <- matching_row_depth$p.sppa[1]
      p.sppa_lower <- ifelse(is.na(matching_row_depth$p.sppa_lower[1]), NA, matching_row_depth$p.sppa_lower[1])
      p.sppa_upper <- ifelse(is.na(matching_row_depth$p.sppa_upper[1]), NA, matching_row_depth$p.sppa_upper[1])
    }
    
    # Handle cases with no depth match (use pot_no only with warning)
    if (!has_depth_match) {
      matching_row_pot_no <- mixture[mixture$pot_no == current_pot_no, ]
      if (nrow(matching_row_pot_no) > 0) {
        warning(paste0("No matching depth found for pot_no:", current_pot_no, ". Using pot_no only."))
        p.scam_value <- matching_row_pot_no$p.scam[1]
        p.scam_lower <- ifelse(is.na(matching_row_pot_no$p.scam_lower[1]), NA, matching_row_pot_no$p.scam_lower[1])
        p.scam_upper <- ifelse(is.na(matching_row_pot_no$p.scam_upper[1]), NA, matching_row_pot_no$p.scam_upper[1])
        
        p.sppa_value <- matching_row_pot_no$p.sppa[1]
        sppa_lower <- ifelse(is.na(matching_row_pot_no$p.sppa_lower[1]), NA, matching_row_pot_no$sppa_lower[1])
        sppa_upper <- ifelse(is.na(matching_row_pot_no$p.sppa_upper[1]), NA, matching_row_pot_no$sppa_upper[1])
      } else {
        warning(paste0("No matching row found for pot_no:", current_pot_no, ". Setting prop_scam to NA."))
      }
    }
    
    # Calculate prop_scam and prop_sppa
    bgb_biomass$prop_scam[i] <- bgb_biomass$weight_density[i] * p.scam_value
    bgb_biomass$prop_scam_lower[i] <- bgb_biomass$weight_density[i] * p.scam_lower
    bgb_biomass$prop_scam_upper[i] <- bgb_biomass$weight_density[i] * p.scam_upper
    
    bgb_biomass$prop_sppa[i] <- bgb_biomass$weight_density[i] * p.sppa_value
    bgb_biomass$prop_sppa_lower[i] <- bgb_biomass$weight_density[i] * p.sppa_lower
    bgb_biomass$prop_sppa_upper[i] <- bgb_biomass$weight_density[i] * p.sppa_upper
  }
  
  # Return the updated bgb_biomass dataframe
  return(bgb_biomass)
}
    