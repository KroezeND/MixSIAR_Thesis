library(MixSIAR)
library(tidyverse)
library(here)
here::here()

mix.filename <-"output/data/training.csv"
source.filename="output/data/marsh_source_all.csv"
discr.filename = "data/marsh_discr.csv"


n.mod <- 8
mix <- vector("list", n.mod) 
mix[[1]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C"),
                          factors=NULL,
                          fac_random=NULL,
                          fac_nested=NULL,
                          cont_effects=NULL)
mix[[2]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C"),
                          factors="env_treatment",
                          fac_random=FALSE,
                          fac_nested=FALSE,
                          cont_effects=NULL)
mix[[3]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C"),
                          factors=NULL,
                          fac_random=NULL,
                          fac_nested=NULL,
                          cont_effects="seed_year")
mix[[4]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C"),
                          factors=c("env_treatment"),
                          fac_random=c(FALSE),
                          fac_nested=c(FALSE),
                          cont_effects="seed_year")
mix[[5]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C"),
                          factors=c("env_treatment","Depth"),
                          fac_random=c(FALSE,FALSE),
                          fac_nested=c(FALSE,FALSE),
                          cont_effects="seed_year")
mix[[6]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C"),
                          factors=c("Depth"),
                          fac_random=c(FALSE),
                          fac_nested=c(FALSE),
                          cont_effects="seed_year")
mix[[7]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C"),
                          factors=c("env_treatment","Depth"),
                          fac_random=c(FALSE,FALSE),
                          fac_nested=c(FALSE,FALSE),
                          cont_effects=NULL)
mix[[8]] <- load_mix_data(filename=mix.filename,
                          iso_names=c("d13C"),
                          factors=c("env_treatment","pot_no"),
                          fac_random=c(FALSE,TRUE),
                          fac_nested=c(FALSE,FALSE),
                          cont_effects=NULL)


# Run the models
source <- vector("list", n.mod)
discr <- vector("list", n.mod)
jags.mod <- vector("list", n.mod)
for(mod in 1:n.mod){
  # create sub-directory and move into it
  source[[mod]] <- load_source_data(filename=source.filename,
                                    source_factors=NULL,
                                    conc_dep=FALSE,
                                    data_type="raw",
                                    mix[[mod]])
  
  # load TEF data
  discr[[mod]] <- load_discr_data(filename=discr.filename, mix[[mod]])
  
  mainDir <- getwd()
  subDir <- paste0("model_", mod)
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
  setwd(file.path(mainDir, subDir))
  
  # isospace plot
  plot_data(filename=paste0("isospace_plot_", mod),
            plot_save_pdf=TRUE,
            plot_save_png=FALSE,
            mix[[mod]], source[[mod]], discr[[mod]])
  
  # Define model structure and write JAGS model file
  model_filename <- paste0("MixSIAR_model_", mod, ".txt")
  resid_err <- TRUE
  process_err <- TRUE
  write_JAGS_model(model_filename, resid_err, process_err, mix[[mod]], source[[mod]])
  
  # Run the JAGS model
  # "short" MCMC length is plenty long for all models to converge
  jags.mod[[mod]] <- run_model(run="test", mix[[mod]], source[[mod]], discr[[mod]], model_filename, alpha.prior=1)
  
  # Process diagnostics, summary stats, and posterior plots
  output_options=list(
    summary_save = TRUE,                 # Save the summary statistics as a txt file?
    summary_name = "summary_statistics",    # If yes, specify the base file name (.txt will be appended later)
    sup_post = TRUE,                       # Suppress posterior density plot output in R?
    plot_post_save_pdf = TRUE,              # Save posterior density plots as pdfs?
    plot_post_name = "posterior_density",   # If yes, specify the base file name(s) (.pdf/.png will be appended later)
    sup_pairs = FALSE,                      # Suppress pairs plot output in R?
    plot_pairs_save_pdf = TRUE,             # Save pairs plot as pdf?
    plot_pairs_name = "pairs_plot",         # If yes, specify the base file name (.pdf/.png will be appended later)
    sup_xy = TRUE,                         # Suppress xy/trace plot output in R?
    plot_xy_save_pdf = FALSE,                # Save xy/trace plot as pdf?
    plot_xy_name = "xy_plot",               # If yes, specify the base file name (.pdf/.png will be appended later)
    gelman = TRUE,                          # Calculate Gelman-Rubin diagnostic test?
    heidel = FALSE,                          # Calculate Heidelberg-Welch diagnostic test?
    geweke = TRUE,                          # Calculate Geweke diagnostic test?
    diag_save = TRUE,                       # Save the diagnostics as a txt file?
    diag_name = "diagnostics",              # If yes, specify the base file name (.txt will be appended later)
    indiv_effect = FALSE,                   # Is Individual a random effect in the model? (already specified)
    plot_post_save_png = FALSE,             # Save posterior density plots as pngs?
    plot_pairs_save_png = FALSE,            # Save pairs plot as png?
    plot_xy_save_png = FALSE,
    diag_save_ggmcmc = FALSE,
    return_obj = FALSE)
  
  output_JAGS(jags.mod[[mod]], mix[[mod]], source[[mod]], output_options)
  graphics.off()
  # Move back up to root directory
  setwd(mainDir)
}

# Use 'compare_models' to get table with LOOic weights
names(jags.mod) <- c("Environment:Midpoint","Environment:Seed Year","Environment/Pot No",
                     "Environment/Cohort","Environment/Cohort:Midpoint",
                     "Environment/Cohort:Seed Year","Environment/Pot No:Midpoint",
                     "Pot No:Midpoint")
comparison.table <- compare_models(jags.mod)

# get multiplicative error term estimates, median(xi.C)
xi.C <- rep(NA, n.mod)
for(i in 1:n.mod){
  xi.C[i] <- round(median(jags.mod[[i]]$BUGSoutput$sims.list$resid.prop[,1]),1)
}

# add xi.C to comparison table 
y <- as.numeric(rownames(comparison.table))
# DIC <- c(647.5,610.7,634.9,647.5,649.0)
comparison.table <- cbind(comparison.table, xi.C[y])
comparison.table

saveRDS(jags.1,"jags.Depth.rds")
