# Function to generate predictions for a fitted model and store in list
predict_and_store <- function(data, year_data, year) {
  # Filter data for current year
  filtered_data <- data[data$seed_year == year, ]
  
  # Fit model
  model <- lm(log(weight_density) ~ midpoint, data = filtered_data)

  slope <- model$coefficients[2]
  intercept <- model$coefficients[1]
  
  x_values <- seq(min(data$midpoint), max(data$midpoint), length = 100)
  
  predicted_data <- data.frame(x = x_values, y = slope * x_values + intercept)
  predicted_data$year <- year
  
  # Store results in list (append)
  year_data[[year]] <- predicted_data
  year_data[[year]]$equation[1] <- paste("y =", round(slope,3), "x + ", round(intercept,3), sep = "")

  # Return the updated list (optional)
  return(year_data)
}
