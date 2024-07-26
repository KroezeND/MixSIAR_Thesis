# Author: N. Kroeze
# Description: This script is meant to set the characteristics of the mixing model, define using MixSIAR,
# run the model, and save the outputs as intermediary products

# Inputs:  "marsh_source_all.csv","marsh_discr.csv","training_TOP.csv","training_MID.csv","training_BTM.csv"
# Outputs: "Depth/jags.EnvGT.rds",
# outputs/Depth/model_X:
#               "isospace_env_plot.pdf"
#               "prior_env_plot.pdf,
#               "diagnostics.pdf"
#               "diagnostics.txt"
#               "summary_statistics.txt",
#               "MixSIAR_model.txt" 

# Required but not called by MixSIAR
require(ggplot2)

(rm(list = ls()))

source.filename = "data/marsh_source_all.csv"
discr.filename = "data/marsh_discr.csv"


n.mod <- 3
mix <- vector("list", n.mod) 
mix[[1]] <- MixSIAR::load_mix_data(filename="data/training_TOP.csv",
                          iso_names=c("d13C"),
                          factors=c("env_treatment","gt"),
                          fac_random=c(FALSE,TRUE),
                          fac_nested=c(FALSE,FALSE),
                          cont_effects="seed_year")
mix[[2]] <- MixSIAR::load_mix_data(filename="data/training_MID.csv",
                          iso_names=c("d13C"),
                          factors=c("env_treatment","gt"),
                          fac_random=c(FALSE,TRUE),
                          fac_nested=c(FALSE,FALSE),
                          cont_effects="seed_year")
mix[[3]] <- MixSIAR::load_mix_data(filename="data/training_BTM.csv",
                          iso_names=c("d13C"),
                          factors=c("env_treatment","gt"),
                          fac_random=c(FALSE,TRUE),
                          fac_nested=c(FALSE,FALSE),
                          cont_effects="seed_year")

# Run the models
source <- vector("list", n.mod)
discr <- vector("list", n.mod)
jags.EnvGT <- vector("list", n.mod)
for(mod in 1:n.mod){
  
  # create sub-directory and move into it
  source[[mod]] <- MixSIAR::load_source_data(filename=source.filename,
                                    source_factors=NULL,
                                    conc_dep=FALSE,
                                    data_type="raw",
                                    mix[[mod]])
  
  # load TEF data
  discr[[mod]] <- MixSIAR::load_discr_data(filename=discr.filename, mix[[mod]])
  
  subDir <- paste0("model_", mod)
  dir.create(file.path('output/Depth', subDir), showWarnings = FALSE)

  # isospace plot
  MixSIAR::plot_data(filename=paste0("output/Depth/", subDir, "/isospace_plot_", mod),
            plot_save_pdf=TRUE,
            plot_save_png=FALSE,
            mix[[mod]], source[[mod]], discr[[mod]])
  
  # Define model structure and write JAGS model file
  model_filename <- paste0("output/Depth/", subDir, "/MixSIAR_model_", mod, ".txt")
  resid_err <- TRUE
  process_err <- TRUE
  MixSIAR::write_JAGS_model(model_filename, resid_err, process_err, mix[[mod]], source[[mod]])
  
  # Run the described model and save as jags object
  # run == 	   Chain Length  Burn-in 	 Thin 	# Chains  # NK Run Time
  # “test” 	     1,000 	     500 	      1 	  3         Few seconds
  # “very short” 10,000 	   5,000 	    5 	  3         Around 1-3 min
  # “short” 	   50,000   	 25,000 	  25 	  3         Maybe 5 min or so
  # “normal” 	   100,000 	   50,000 	  50 	  3         Varies but likely 10-15 min, used for analysis going forward
  # “long” 	     300,000 	   200,000 	  100 	3         Untested
  # “very long”  1,000,000 	 500,000 	  500 	3         Untested
  # “extreme” 	 3,000,000 	 1,500,000  500 	3         Untested
  
  # Run the JAGS model
  jags.EnvGT[[mod]] <- MixSIAR::run_model(run="normal", mix[[mod]], source[[mod]], discr[[mod]], model_filename, alpha.prior=1)
 
  output_options <- list(
    summary_save = TRUE,                  # save summary statistics as txt file
    summary_name = "summary_statistics",  # path for summary stats
    sup_post = TRUE,                      # suppress posterior density plot, default is FALSE
    plot_post_save_pdf = TRUE,            # save posterior density plot
    plot_post_name = "posterior_density", # path for posterior density plots
    sup_pairs = TRUE,                     # suppress pairs plot
    plot_pairs_save_pdf = TRUE,           # save pairs plot
    plot_pairs_name = "pairs_plot",       # path for pairs plot
    sup_xy = TRUE,                        # suppress x/y trace plot
    plot_xy_save_pdf = TRUE,              # save x/y trace plot
    plot_xy_name = "xy_plot",             # path for x/y trace plot
    gelman = TRUE,                        # calculate Gelman-Rubin
    heidel = TRUE,                        # calculate Heidelberg-Welch
    geweke = TRUE,                        # calculate Geweke 
    diag_save = TRUE,                     # save diagnostics as txt file
    diag_name = "diagnostics",            # path for diagnostics
    indiv_effect = FALSE,                 # artifact, set to FALSE
    plot_post_save_png = FALSE,           # save posterior density plot as png
    plot_pairs_save_png = FALSE,          # save pairs plot as png
    plot_xy_save_png = FALSE,             # save x/y trace plot as png
    diag_save_ggmcmc = FALSE)             # save ggmcmc diagnostics as pdf
  
  
  MixSIAR::output_JAGS(jags.EnvGT[[mod]], mix[[mod]], source[[mod]], output_options)
  graphics.off()
}


# save jags.Env object for further analysis
saveRDS(jags.EnvGT,"output/Depth/jags.EnvGT.rds")

