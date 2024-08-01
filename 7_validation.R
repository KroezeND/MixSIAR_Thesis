# Author: N. Kroeze
# Description: This script is meant to set the characteristics of the mixing model, define using MixSIAR,
# run the model, and save the outputs as intermediary products

# Inputs:  "training.csv","marsh_source_all.csv","marsh_discr.csv"
# Outputs: "jags.EnvGT.rds","figs/isospace_env_plot.pdf","figs/prior_env_plot.pdf,
#          "diagnostics.pdf","diagnostics.txt","summary_statistics.txt",
#          "MixSIAR_model.txt","posterior_density_high_sal","posterior_density_low_sal"  

library(MixSIAR)
library(tidyverse)
library(here)

(rm(list = ls()))
source("src/data_process.R")

# Formatting holdout data and averaging across replicate samples, saving in data folder
holdout <- read_csv("data/Holdout/holdout_DONOTTOUCH.csv")
# Divide by depth layer (mixture distribution)
holdout_TOP <- holdout[holdout$Depth=="TOP",]
holdout_MID <- holdout[holdout$Depth=="MID",]
holdout_BTM <- holdout[holdout$Depth=="BTM",]

# # Average across replicates to utilize all available information
new_holdout_TOP <- average_pot_no(holdout_TOP)
new_holdout_MID <- average_pot_no(holdout_MID)
new_holdout_BTM <- average_pot_no(holdout_BTM)

write.csv(new_holdout_TOP,"data/holdout_avg_TOP.csv",row.names = FALSE)
write.csv(new_holdout_MID,"data/holdout_avg_MID.csv",row.names = FALSE)
write.csv(new_holdout_BTM,"data/holdout_avg_BTM.csv",row.names = FALSE)
#########################################################################
# Set holdout d13C to NA to predict d13C in R2Jags function
training_TOP <- read.csv("data/training_TOP.csv")
training_MID <- read.csv("data/training_MID.csv")
training_BTM <- read.csv("data/training_BTM.csv")

holdout$d13C <- NA

# Divide by depth layer (mixture distribution)
validation_TOP <- holdout[holdout$Depth=="TOP",]
validation_MID <- holdout[holdout$Depth=="MID",]
validation_BTM <- holdout[holdout$Depth=="BTM",]

# # Average across replicates to utilize all available information
new_validation_TOP <- average_pot_no(validation_TOP) # doesn't functionally change NA, but ensures only one replicate per pot_no
new_validation_MID <- average_pot_no(validation_MID) # doesn't functionally change NA, but ensures only one replicate per pot_no
new_validation_BTM <- average_pot_no(validation_BTM) # doesn't functionally change NA, but ensures only one replicate per pot_no

validation_TOP <- rbind(training_TOP,new_validation_TOP)
validation_MID <- rbind(training_MID,new_validation_MID)
validation_BTM <- rbind(training_BTM,new_validation_BTM)

write.csv(validation_TOP,"data/validation_TOP.csv",row.names = FALSE)
write.csv(validation_MID,"data/validation_MID.csv",row.names = FALSE)
write.csv(validation_BTM,"data/validation_BTM.csv", row.names = FALSE)


(rm(list = ls()))
######################################### TOP ##################################################
## Predicting d13C for top layer only

# In-sample + validation data
TOP <- read.csv('data/validation_TOP.csv')

# e parameter calculated in the same way as in MixSIAR source code
e_save <- matrix(rep(0, 2 * (2 - 1)), nrow = 2, ncol = (2 - 1))
e_save[,1] <- exp(c(rep(sqrt(1/(1 * (1 + 1))), 1), -sqrt(1/(1 + 1)), rep(0, 2 - 1 - 1)))
e_save[,1] <- e_save[,1] / sum(e_save[,1])

# Making mix object from which we take some of the input values for JAGS
mix <- MixSIAR::load_mix_data(filename = 'data/validation_TOP.csv', # new in-sample + validation data
                              iso_names = 'd13C', # name of column with isotope data
                              factors = c('env_treatment', 'gt'), # factor  covariates
                              fac_random = c(FALSE, TRUE), # specify which is random factor
                              fac_nested = c(FALSE, FALSE), # neither factor is nested
                              cont_effects = 'seed_year') # continuous effect covariate

# Making source object from which we take some of the input values for JAGS
source <- MixSIAR::load_source_data(filename = 'data/marsh_source_all.csv', # same source data as for original model fitting
                                    source_factors = NULL, # no source factors
                                    conc_dep = FALSE,
                                    data_type = 'raw', # raw values for source data
                                    mix = mix) # mix object from above

# Making discr object from which we take some of the input values for JAGS
discr <- MixSIAR::load_discr_data(filename = 'data/marsh_discr.csv', # same discr data as for the original model fitting
                                  mix) # mix object from above

# List of inputs for JAGS
dat <- list(X_iso = as.matrix(TOP$d13C), # isotope data (missing data are the ones we want to predict)
            N = mix$N, # from mix object
            n.sources = source$n.sources, #  from source object
            n.iso = mix$n.iso, #  from mix object
            alpha = rep(1, 2), # defined in MixSIAR source code
            frac_mu = discr$mu, # from discr object
            frac_sig2 = discr$sig2, # from discr object
            e = e_save, # defined above
            cross = array(data = NA, dim = c(nrow(TOP), 2, 2 - 1)), # defined in MixSIAR source code
            tmp.p = array(data = NA, dim = c(nrow(TOP), 2)), # defined in MixSIAR source code
            factor1_levels = length(unique(TOP$env_treatment)), # defined in MixSIAR source code
            cross.fac1 = array(data = NA, dim = c(length(unique(TOP$env_treatment)), 2, 2 - 1)), # defined in MixSIAR source code
            factor2_levels = length(unique(TOP$gt)), # defined in MixSIAR source code
            n.rep = source$n.rep, # from source object
            SOURCE_array = source$SOURCE_array, # from source object
            Factor.1 = mix$FAC[[1]]$values, # from mix object
            Factor.2 = mix$FAC[[2]]$values, # from mix object
            Cont.1 = as.vector(mix$CE[[1]])) # from mix object, using scaled & centered variable

# define initial values in the same way as in MixSIAR source code
init <- function(){list(p.global = as.vector(MCMCpack::rdirichlet(1, rep(1, 2))))}

# Same params as specified in MixSIAR source code, plus X_iso, what we want to predict
params <- c('X_iso', 'p.global', 'loglik', 'p.fac1', 'fac2.sig', 'ilr.global', 'ilr.fac1',
            'ilr.fac2', 'ilr.global', 'ilr.cont1', 'p.ind', 'resid.prop', 'deviance')

# Run JAGS model
jags_TOP <- R2jags::jags(data = dat, # Data list
                         inits = init, # inits function
                         parameters.to.save = params, # parameters to save in output
                         model.file = 'output/Depth/model_1/MixSIAR_model_1.txt', # file containing JAGS model written by MixSIAR
                         n.chains = 3, # same as for "normal" run in MixSIAR
                         n.iter = 100000, # same as for "normal" run in MixSIAR
                         n.burnin = 50000, # same as for "normal" run in MixSIAR
                         n.thin = 50) # same as for "normal" run in MixSIAR
measured_d13C <- read.csv("data/holdout_avg_TOP.csv")
estimated_d13C <- read.csv("data/validation_TOP.csv")

pred_x_iso_TOP <- jags_TOP$BUGSoutput$sims.list$X_iso # predictions for X_iso
pred_x_iso_TOP <- pred_x_iso_TOP[,,1] # make matrix instead of array
pred_x_iso_TOP <- pred_x_iso_TOP[,100:125]

# Transpose the matrix
pred_x_iso_transposed <- t(pred_x_iso_TOP)

# Convert to a data frame
pred_x_iso_df <- as.data.frame(pred_x_iso_transposed)

# Create a long format data frame
pred_x_iso_long <- pred_x_iso_df |>
  rownames_to_column(var = "pot") |>
  pivot_longer(cols = -pot, names_to = "draw", values_to = "value")
pot_no <- measured_d13C$pot_no
# Create a mapping between original pot numbers and desired pot numbers
pot_mapping <- data.frame(original_pot = seq_along(pot_no), desired_pot = pot_no)

# Join the mapping to the pred_x_iso_long data frame
pred_x_iso_long <- pred_x_iso_long |>
  mutate(pot = as.numeric(pot)) |>  # Convert pot to numeric for joining
  left_join(pot_mapping, by = c("pot" = "original_pot")) |>
  select(-pot) |>  # Remove the original pot column
  rename(pot = desired_pot)

# Extract relevant columns from estimated_d13C
pot_cohort <- estimated_d13C |>
  select(pot_no, cohort, env_treatment, seed_year)

# Join cohort information to pred_x_iso_long
pred_x_iso_long <- pred_x_iso_long |>
  left_join(pot_cohort, by = c("pot" = "pot_no"))

unique_cohort <- sort(unique(estimated_d13C$cohort))
colors <- c("#2d6a4f","#74c69d","#5a189a","#c77dff")
color_map <- setNames(colors, unique_cohort)

# Create the violin plot faceted by pot
ggplot(pred_x_iso_long, aes(x = as.factor(seed_year), y = value, fill = cohort)) +
  geom_violin(alpha = 0.7) +
  facet_wrap(~ env_treatment) +
  labs(x = "Seed Age", y = "Posterior Value of  δ¹³C (‰)") +
  scale_fill_manual(values = color_map) +
  theme_bw()

################ BIAS ###################
d13C <- jags_TOP$BUGSoutput$sims.list$X_iso

# Add mean mixing proportion of each individual to mixture dataset
pred_x_iso_mean <- as.data.frame(jags_TOP$BUGSoutput$mean$X_iso)
# Calculate lower and upper bounds (assuming 95% CrI)
pred_x_iso_mean$d13C_lower <- apply(d13C, 2, FUN = quantile, probs = 0.025)
pred_x_iso_mean$d13C_upper <- apply(d13C, 2, FUN = quantile, probs = 0.975)

estimated_d13C <- estimated_d13C[100:125,]
estimated_d13C$d13C <- pred_x_iso_mean$V1[100:125]
estimated_d13C$d13C_lower <- pred_x_iso_mean$d13C_lower[100:125]
estimated_d13C$d13C_upper <- pred_x_iso_mean$d13C_upper[100:125]

gg_TOP <- ggplot(data = estimated_d13C, aes(x = seed_year,y=d13C, color = "black")) +
  geom_point() +
  geom_errorbar(aes(ymin = d13C_lower, ymax = d13C_upper), width = 0.25) + 
  geom_point(data = measured_d13C, aes(x = seed_year,y = d13C, color = "red"))+
  facet_grid(~ env_treatment) +
  labs(x = "", y = "Posterior Value of  δ¹³C (‰)") +
  theme_bw()   +
  ggplot2::scale_color_manual(values= c("black","red"),labels = c("Estimated","Measured"))+
  guides(color = guide_legend(title = "Age Cohort"))+
  ggplot2::theme(legend.position = "none")

new_data <- data.frame("pot_no" = estimated_d13C$pot_no,
                       "est_d13C" = estimated_d13C$d13C,
                       "real_d13C" = measured_d13C$d13C, 
                       "env_treatment" = estimated_d13C$env_treatment,
                       "seed_year" = estimated_d13C$seed_year,
                       "gt" = estimated_d13C$gt,
                       "cohort" = estimated_d13C$cohort,
                       "est_d13C_lower" = estimated_d13C$d13C_lower,
                       "est_d13C_upper" = estimated_d13C$d13C_upper)


diff_d13C <- data.frame("diff_d13C" = measured_d13C$d13C - estimated_d13C$d13C)

new_data <- cbind(new_data,diff_d13C)


unique_cohort <- sort(unique(estimated_d13C$cohort))
colors <- c("#2d6a4f","#74c69d","#5a189a","#c77dff")
color_map <- setNames(colors, unique_cohort)

ggplot(data = new_data, aes(x = seed_year,y=diff_d13C, color = cohort))+
  geom_point() +
  facet_wrap(~ env_treatment) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")+
  scale_color_manual(values = color_map) +
  theme_bw()

######################################### MID ##################################################
## Predicting d13C for MID layer only

# In-sample + validation data
MID <- read.csv('data/validation_MID.csv')

# e parameter calculated in the same way as in MixSIAR source code
e_save <- matrix(rep(0, 2 * (2 - 1)), nrow = 2, ncol = (2 - 1))
e_save[,1] <- exp(c(rep(sqrt(1/(1 * (1 + 1))), 1), -sqrt(1/(1 + 1)), rep(0, 2 - 1 - 1)))
e_save[,1] <- e_save[,1] / sum(e_save[,1])

# Making mix object from which we take some of the input values for JAGS
mix <- MixSIAR::load_mix_data(filename = 'data/validation_MID.csv', # new in-sample + validation data
                              iso_names = 'd13C', # name of column with isotope data
                              factors = c('env_treatment', 'gt'), # factor  covariates
                              fac_random = c(FALSE, TRUE), # specify which is random factor
                              fac_nested = c(FALSE, FALSE), # neither factor is nested
                              cont_effects = 'seed_year') # continuous effect covariate

# Making source object from which we take some of the input values for JAGS
source <- MixSIAR::load_source_data(filename = 'data/marsh_source_all.csv', # same source data as for original model fitting
                                    source_factors = NULL, # no source factors
                                    conc_dep = FALSE,
                                    data_type = 'raw', # raw values for source data
                                    mix = mix) # mix object from above

# Making discr object from which we take some of the input values for JAGS
discr <- MixSIAR::load_discr_data(filename = 'data/marsh_discr.csv', # same discr data as for the original model fitting
                                  mix) # mix object from above

# List of inputs for JAGS
dat <- list(X_iso = as.matrix(MID$d13C), # isotope data (missing data are the ones we want to predict)
            N = mix$N, # from mix object
            n.sources = source$n.sources, #  from source object
            n.iso = mix$n.iso, #  from mix object
            alpha = rep(1, 2), # defined in MixSIAR source code
            frac_mu = discr$mu, # from discr object
            frac_sig2 = discr$sig2, # from discr object
            e = e_save, # defined above
            cross = array(data = NA, dim = c(nrow(MID), 2, 2 - 1)), # defined in MixSIAR source code
            tmp.p = array(data = NA, dim = c(nrow(MID), 2)), # defined in MixSIAR source code
            factor1_levels = length(unique(MID$env_treatment)), # defined in MixSIAR source code
            cross.fac1 = array(data = NA, dim = c(length(unique(MID$env_treatment)), 2, 2 - 1)), # defined in MixSIAR source code
            factor2_levels = length(unique(MID$gt)), # defined in MixSIAR source code
            n.rep = source$n.rep, # from source object
            SOURCE_array = source$SOURCE_array, # from source object
            Factor.1 = mix$FAC[[1]]$values, # from mix object
            Factor.2 = mix$FAC[[2]]$values, # from mix object
            Cont.1 = as.vector(mix$CE[[1]])) # from mix object, using scaled & centered variable

# define initial values in the same way as in MixSIAR source code
init <- function(){list(p.global = as.vector(MCMCpack::rdirichlet(1, rep(1, 2))))}

# Same params as specified in MixSIAR source code, plus X_iso, what we want to predict
params <- c('X_iso', 'p.global', 'loglik', 'p.fac1', 'fac2.sig', 'ilr.global', 'ilr.fac1',
            'ilr.fac2', 'ilr.global', 'ilr.cont1', 'p.ind', 'resid.prop', 'deviance')

# Run JAGS model
jags_MID <- R2jags::jags(data = dat, # Data list
                         inits = init, # inits function
                         parameters.to.save = params, # parameters to save in output
                         model.file = 'output/Depth/model_1/MixSIAR_model_1.txt', # file containing JAGS model written by MixSIAR
                         n.chains = 3, # same as for "normal" run in MixSIAR
                         n.iter = 100000, # same as for "normal" run in MixSIAR
                         n.burnin = 50000, # same as for "normal" run in MixSIAR
                         n.thin = 50) # same as for "normal" run in MixSIAR
measured_d13C <- read.csv("data/holdout_avg_MID.csv")
estimated_d13C <- read.csv("data/validation_MID.csv")

pred_x_iso_MID <- jags_MID$BUGSoutput$sims.list$X_iso # predictions for X_iso
pred_x_iso_MID <- pred_x_iso_MID[,,1] # make matrix instead of array
pred_x_iso_MID <- pred_x_iso_MID[,100:125]

# Transpose the matrix
pred_x_iso_transposed <- t(pred_x_iso_MID)

# Convert to a data frame
pred_x_iso_df <- as.data.frame(pred_x_iso_transposed)

# Create a long format data frame
pred_x_iso_long <- pred_x_iso_df |>
  rownames_to_column(var = "pot") |>
  pivot_longer(cols = -pot, names_to = "draw", values_to = "value")
pot_no <- measured_d13C$pot_no
# Create a mapping between original pot numbers and desired pot numbers
pot_mapping <- data.frame(original_pot = seq_along(pot_no), desired_pot = pot_no)

# Join the mapping to the pred_x_iso_long data frame
pred_x_iso_long <- pred_x_iso_long |>
  mutate(pot = as.numeric(pot)) |>  # Convert pot to numeric for joining
  left_join(pot_mapping, by = c("pot" = "original_pot")) |>
  select(-pot) |>  # Remove the original pot column
  rename(pot = desired_pot)

# Extract relevant columns from estimated_d13C
pot_cohort <- estimated_d13C |>
  select(pot_no, cohort, env_treatment, seed_year)

# Join cohort information to pred_x_iso_long
pred_x_iso_long <- pred_x_iso_long |>
  left_join(pot_cohort, by = c("pot" = "pot_no"))

unique_cohort <- sort(unique(estimated_d13C$cohort))
colors <- c("#2d6a4f","#74c69d","#5a189a","#c77dff")
color_map <- setNames(colors, unique_cohort)

# Create the violin plot faceted by pot
ggplot(pred_x_iso_long, aes(x = as.factor(seed_year), y = value, fill = cohort)) +
  geom_violin(alpha = 0.7) +
  facet_wrap(~ env_treatment) +
  labs(x = "Seed Age", y = "Posterior Value of  δ¹³C (‰)") +
  scale_fill_manual(values = color_map) +
  theme_bw()

################ BIAS ###################
d13C <- jags_MID$BUGSoutput$sims.list$X_iso

# Add mean mixing proportion of each individual to mixture dataset
pred_x_iso_mean <- as.data.frame(jags_MID$BUGSoutput$mean$X_iso)
# Calculate lower and upper bounds (assuming 95% CrI)
pred_x_iso_mean$d13C_lower <- apply(d13C, 2, FUN = quantile, probs = 0.025)
pred_x_iso_mean$d13C_upper <- apply(d13C, 2, FUN = quantile, probs = 0.975)

estimated_d13C <- estimated_d13C[100:125,]
estimated_d13C$d13C <- pred_x_iso_mean$V1[100:125]
estimated_d13C$d13C_lower <- pred_x_iso_mean$d13C_lower[100:125]
estimated_d13C$d13C_upper <- pred_x_iso_mean$d13C_upper[100:125]

gg_MID <- ggplot(data = estimated_d13C, aes(x = seed_year,y=d13C, color = "black")) +
  geom_point() +
  geom_errorbar(aes(ymin = d13C_lower, ymax = d13C_upper), width = 0.25) + 
  geom_point(data = measured_d13C, aes(x = seed_year,y = d13C, color = "red"))+
  facet_grid(~ env_treatment) +
  labs(x ="", y = "Posterior Value of  δ¹³C (‰)") +
  theme_bw()   +
  ggplot2::scale_color_manual(values= c("black","red"),labels = c("Estimated","Measured"))+
  guides(color = guide_legend(title = "Age Cohort"))+
  ggplot2::theme(legend.position = "none")

new_data <- data.frame("pot_no" = estimated_d13C$pot_no,
                       "est_d13C" = estimated_d13C$d13C,
                       "real_d13C" = measured_d13C$d13C, 
                       "env_treatment" = estimated_d13C$env_treatment,
                       "seed_year" = estimated_d13C$seed_year,
                       "gt" = estimated_d13C$gt,
                       "cohort" = estimated_d13C$cohort,
                       "est_d13C_lower" = estimated_d13C$d13C_lower,
                       "est_d13C_upper" = estimated_d13C$d13C_upper)


diff_d13C <- data.frame("diff_d13C" = measured_d13C$d13C - estimated_d13C$d13C)

new_data <- cbind(new_data,diff_d13C)


unique_cohort <- sort(unique(estimated_d13C$cohort))
colors <- c("#2d6a4f","#74c69d","#5a189a","#c77dff")
color_map <- setNames(colors, unique_cohort)

ggplot(data = new_data, aes(x = seed_year,y=diff_d13C, color = cohort))+
  geom_point() +
  facet_wrap(~ env_treatment) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")+
  scale_color_manual(values = color_map) +
  theme_bw()
######################################### BTM ##################################################
## Predicting d13C for BTM layer only

# In-sample + validation data
BTM <- read.csv('data/validation_BTM.csv')

# e parameter calculated in the same way as in MixSIAR source code
e_save <- matrix(rep(0, 2 * (2 - 1)), nrow = 2, ncol = (2 - 1))
e_save[,1] <- exp(c(rep(sqrt(1/(1 * (1 + 1))), 1), -sqrt(1/(1 + 1)), rep(0, 2 - 1 - 1)))
e_save[,1] <- e_save[,1] / sum(e_save[,1])

# Making mix object from which we take some of the input values for JAGS
mix <- MixSIAR::load_mix_data(filename = 'data/validation_BTM.csv', # new in-sample + validation data
                              iso_names = 'd13C', # name of column with isotope data
                              factors = c('env_treatment', 'gt'), # factor  covariates
                              fac_random = c(FALSE, TRUE), # specify which is random factor
                              fac_nested = c(FALSE, FALSE), # neither factor is nested
                              cont_effects = 'seed_year') # continuous effect covariate

# Making source object from which we take some of the input values for JAGS
source <- MixSIAR::load_source_data(filename = 'data/marsh_source_all.csv', # same source data as for original model fitting
                                    source_factors = NULL, # no source factors
                                    conc_dep = FALSE,
                                    data_type = 'raw', # raw values for source data
                                    mix = mix) # mix object from above

# Making discr object from which we take some of the input values for JAGS
discr <- MixSIAR::load_discr_data(filename = 'data/marsh_discr.csv', # same discr data as for the original model fitting
                                  mix) # mix object from above

# List of inputs for JAGS
dat <- list(X_iso = as.matrix(BTM$d13C), # isotope data (missing data are the ones we want to predict)
            N = mix$N, # from mix object
            n.sources = source$n.sources, #  from source object
            n.iso = mix$n.iso, #  from mix object
            alpha = rep(1, 2), # defined in MixSIAR source code
            frac_mu = discr$mu, # from discr object
            frac_sig2 = discr$sig2, # from discr object
            e = e_save, # defined above
            cross = array(data = NA, dim = c(nrow(BTM), 2, 2 - 1)), # defined in MixSIAR source code
            tmp.p = array(data = NA, dim = c(nrow(BTM), 2)), # defined in MixSIAR source code
            factor1_levels = length(unique(BTM$env_treatment)), # defined in MixSIAR source code
            cross.fac1 = array(data = NA, dim = c(length(unique(BTM$env_treatment)), 2, 2 - 1)), # defined in MixSIAR source code
            factor2_levels = length(unique(BTM$gt)), # defined in MixSIAR source code
            n.rep = source$n.rep, # from source object
            SOURCE_array = source$SOURCE_array, # from source object
            Factor.1 = mix$FAC[[1]]$values, # from mix object
            Factor.2 = mix$FAC[[2]]$values, # from mix object
            Cont.1 = as.vector(mix$CE[[1]])) # from mix object, using scaled & centered variable

# define initial values in the same way as in MixSIAR source code
init <- function(){list(p.global = as.vector(MCMCpack::rdirichlet(1, rep(1, 2))))}

# Same params as specified in MixSIAR source code, plus X_iso, what we want to predict
params <- c('X_iso', 'p.global', 'loglik', 'p.fac1', 'fac2.sig', 'ilr.global', 'ilr.fac1',
            'ilr.fac2', 'ilr.global', 'ilr.cont1', 'p.ind', 'resid.prop', 'deviance')

# Run JAGS model
jags_BTM <- R2jags::jags(data = dat, # Data list
                         inits = init, # inits function
                         parameters.to.save = params, # parameters to save in output
                         model.file = 'output/Depth/model_1/MixSIAR_model_1.txt', # file containing JAGS model written by MixSIAR
                         n.chains = 3, # same as for "normal" run in MixSIAR
                         n.iter = 100000, # same as for "normal" run in MixSIAR
                         n.burnin = 50000, # same as for "normal" run in MixSIAR
                         n.thin = 50) # same as for "normal" run in MixSIAR
measured_d13C <- read.csv("data/holdout_avg_BTM.csv")
estimated_d13C <- read.csv("data/validation_BTM.csv")

pred_x_iso_BTM <- jags_BTM$BUGSoutput$sims.list$X_iso # predictions for X_iso
pred_x_iso_BTM <- pred_x_iso_BTM[,,1] # make matrix instead of array
pred_x_iso_BTM <- pred_x_iso_BTM[,100:125]

# Transpose the matrix
pred_x_iso_transposed <- t(pred_x_iso_BTM)

# Convert to a data frame
pred_x_iso_df <- as.data.frame(pred_x_iso_transposed)

# Create a long format data frame
pred_x_iso_long <- pred_x_iso_df |>
  rownames_to_column(var = "pot") |>
  pivot_longer(cols = -pot, names_to = "draw", values_to = "value")
pot_no <- measured_d13C$pot_no
# Create a mapping between original pot numbers and desired pot numbers
pot_mapping <- data.frame(original_pot = seq_along(pot_no), desired_pot = pot_no)

# Join the mapping to the pred_x_iso_long data frame
pred_x_iso_long <- pred_x_iso_long |>
  mutate(pot = as.numeric(pot)) |>  # Convert pot to numeric for joining
  left_join(pot_mapping, by = c("pot" = "original_pot")) |>
  select(-pot) |>  # Remove the original pot column
  rename(pot = desired_pot)

# Extract relevant columns from estimated_d13C
pot_cohort <- estimated_d13C |>
  select(pot_no, cohort, env_treatment, seed_year)

# Join cohort information to pred_x_iso_long
pred_x_iso_long <- pred_x_iso_long |>
  left_join(pot_cohort, by = c("pot" = "pot_no"))

unique_cohort <- sort(unique(estimated_d13C$cohort))
colors <- c("#2d6a4f","#74c69d","#5a189a","#c77dff")
color_map <- setNames(colors, unique_cohort)

# Create the violin plot faceted by pot
ggplot(pred_x_iso_long, aes(x = as.factor(seed_year), y = value, fill = cohort)) +
  geom_violin(alpha = 0.7) +
  facet_wrap(~ env_treatment) +
  labs(x = "Seed Age", y = "Posterior Value of  δ¹³C (‰)") +
  scale_fill_manual(values = color_map) +
  theme_bw()

################ BIAS ###################
d13C <- jags_BTM$BUGSoutput$sims.list$X_iso

# Add mean mixing proportion of each individual to mixture dataset
pred_x_iso_mean <- as.data.frame(jags_BTM$BUGSoutput$mean$X_iso)
# Calculate lower and upper bounds (assuming 95% CrI)
pred_x_iso_mean$d13C_lower <- apply(d13C, 2, FUN = quantile, probs = 0.025)
pred_x_iso_mean$d13C_upper <- apply(d13C, 2, FUN = quantile, probs = 0.975)

estimated_d13C <- estimated_d13C[100:125,]
estimated_d13C$d13C <- pred_x_iso_mean$V1[100:125]
estimated_d13C$d13C_lower <- pred_x_iso_mean$d13C_lower[100:125]
estimated_d13C$d13C_upper <- pred_x_iso_mean$d13C_upper[100:125]

gg_BTM <- ggplot(data = estimated_d13C, aes(x = seed_year,y=d13C, color = "black")) +
  geom_point() +
  geom_errorbar(aes(ymin = d13C_lower, ymax = d13C_upper), width = 0.25) + 
  geom_point(data = measured_d13C, aes(x = seed_year,y = d13C, color = "red"))+
  facet_grid(~ env_treatment) +
  labs(x = "Seed Age", y = "Posterior Value of  δ¹³C (‰)") +
  theme_bw()   +
  ggplot2::scale_color_manual(values= c("black","red"),labels = c("Estimated","Measured"))+
  guides(color = guide_legend(title = "Age Cohort"))  +
  ggplot2::theme(legend.position = "none")

new_data <- data.frame("pot_no" = estimated_d13C$pot_no,
                       "est_d13C" = estimated_d13C$d13C,
                       "real_d13C" = measured_d13C$d13C, 
                       "env_treatment" = estimated_d13C$env_treatment,
                       "seed_year" = estimated_d13C$seed_year,
                       "gt" = estimated_d13C$gt,
                       "cohort" = estimated_d13C$cohort,
                       "est_d13C_lower" = estimated_d13C$d13C_lower,
                       "est_d13C_upper" = estimated_d13C$d13C_upper)


diff_d13C <- data.frame("diff_d13C" = measured_d13C$d13C - estimated_d13C$d13C)

new_data <- cbind(new_data,diff_d13C)


unique_cohort <- sort(unique(estimated_d13C$cohort))
colors <- c("#2d6a4f","#74c69d","#5a189a","#c77dff")
color_map <- setNames(colors, unique_cohort)

ggplot(data = new_data, aes(x = seed_year,y=diff_d13C, color = cohort))+
  geom_point() +
  facet_wrap(~ env_treatment) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed")+
  scale_color_manual(values = color_map) +
  theme_bw()

gridExtra::grid.arrange(gg_TOP, gg_MID, gg_BTM, nrow=3)  

#######################################################
saveRDS(jags.TOP,"output/Holdout/jags_TOP.rds")
saveRDS(jags.MID,"output/Holdout/jags_MID.rds")
saveRDS(jags.BTM,"output/Holdout/jags_BTM.rds")
#######################################################
