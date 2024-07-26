process_CEST_data <- function(data) {
  # Select and rename columns
  data <- data[c("Identifier.1", "Corr.d13C")] |>
    dplyr::rename(
      ID = Identifier.1,
      d13C = Corr.d13C
    )
  
  # Extract pot number
  data$pot_no <- as.numeric(gsub("2023", "", gsub("\\D", "", data$ID)))
  
  # Remove rows with missing pot numbers (standards - Sorghum flour, peach leaves, and protein std)
  data_sample <- data |> tidyr::drop_na(pot_no)
  
  # Add new columns
  data_sample$Method <- ifelse(grepl("mill", data_sample$ID), "mill",  # Add mill condition
                              ifelse(grepl("coarse", data_sample$ID), "coarse", "fine"))
  data_sample$Instrument <- "ND"
  # data_sample$Species <- gsub(pattern = "(COMP|comp|SPPA|sppa|SCAM|scam){4}", replacement = "\\U\\1", x = data_sample$ID, ignore.case = TRUE)
  data_sample$Species <- stringr::str_extract(data_sample$ID, "[SCAMPO|scampo]{4}")
  data_sample$Species <- gsub(pattern = "scam", replacement = "SCAM", x = data_sample$Species, ignore.case = TRUE)
  data_sample$Tissue <- ifelse(
    grepl("root|stem|rhiz|composite|ROOT|STEM|RHIZ", data_sample$ID, ignore.case = TRUE) &
      !grepl("COMP_", data_sample$ID),
    stringr::str_extract(data_sample$ID, "(root|stem|rhiz|composite)"),
    NA
  )
  data_sample <- data_sample |>
    dplyr::mutate(Tissue = ifelse(is.na(Tissue), "competition", Tissue))
  data_sample$Depth <- ifelse(grepl("TOP", data_sample$ID), "TOP",
                            ifelse(grepl("MID", data_sample$ID), "MID",
                                   ifelse(grepl("BTM", data_sample$ID),"BTM",NA)))
  return(data_sample)
}


process_SI_data <- function(data) {
  # Select and rename columns
  colnames(data) <- c("ID", "d15N", "wt.N", "comments", "d13C", "wt.C")
  
  data <- subset(data, select = -c(d15N,wt.N,comments,wt.C))
  
  # Extract pot number
  data$pot_no <- as.numeric(gsub("2021", "", gsub("\\D", "", data$ID)))
  
  
  # Add new columns
  data$Method <- "fine"
  data$Instrument <- "SI"
  data$Species <- stringr::str_extract(data$ID, "[SCAMPO]{4}")
  data$Tissue <- ifelse(
    grepl("root|stem|rhiz|composite|ROOT|STEM|RHIZ", data$ID) &
      !grepl("COMP_", data$ID),
    stringr::str_extract(data$ID, "(ROOT|STEM|RHIZ|COMPOSITE)"),
    NA
  )
  data <- data |>
    dplyr::mutate(Tissue = ifelse(is.na(Tissue), "competition", Tissue))
  data$Depth <- ifelse(grepl("TOP", data$ID), "TOP",
                              ifelse(grepl("MID|MIDDLE", data$ID), "MID",
                                     ifelse(grepl("BTM|BOTTOM", data$ID),"BTM","COMPOSITE")))
  return(data)
}

# Define a function to handle scenarios with limited data
holdout_sample <- function(scenario_data, target_size) {
  if (nrow(scenario_data) >= target_size) {
    return(dplyr::sample_n(scenario_data, size = target_size, replace = FALSE))
  } else {
    return(scenario_data[sample(1:nrow(scenario_data), size = min_holdout_size, replace = TRUE), ])
  }
}

average_pot_no <- function(training_layer){
  new_training <- data.frame()
  seen_pots <- vector(mode = "character", length = 0)
  
  for (i in 1:nrow(training_layer)) {
    current_pot <- training_layer[i, "pot_no"]
    if (current_pot %in% seen_pots) {
      # Pot already seen, update existing row
      # row_index <- which(training_TOP_avg$pot_no == current_pot)
      # # training_TOP_avg[row_index, "avg_d13C"] <- 
      #   mean(c(training_TOP_avg$d13C[row_index], training_TOP[i, "d13C"]))
    } else {
      # New pot, add a row with averaged d13C
      seen_pots <- c(seen_pots, current_pot)
      new_row <- training_layer[i, ]
      new_row$d13C <- mean(training_layer[training_layer$pot_no == current_pot, "d13C"])
      new_training <- rbind(new_training, new_row)
    }
  }
  return(new_training)
}
