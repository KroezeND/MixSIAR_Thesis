library(MixSIAR)
library(tidyverse)
library(here)
here::here()

# This script is meant to set the characteristics of the mixing model, define using MixSIAR,
# run the model, and save the outputs as intermediary products

# Load in mixture data as a list
mix <- MixSIAR::load_mix_data(filename="output/data/training.csv", 
                              iso_names="d13C", 
                              factors = c("env_treatment"), 
                              fac_random = c(FALSE), 
                              fac_nested = c(FALSE),
                              cont_effects="seed_year")
# Load in source data as a list
source <- MixSIAR::load_source_data(filename="output/data/marsh_source_all.csv", 
                                    source_factors= NULL,  
                                    conc_dep=FALSE, 
                                    data_type="raw", 
                                    mix)
# Dummy discrimination data since there will be no trophic enrichment (all zero)
discr <- MixSIAR::load_discr_data(filename = "data/marsh_discr.csv",mix)

MixSIAR::plot_data_one_iso(mix,source,discr,
                           plot_save_pdf=TRUE, plot_save_png=FALSE, "figs/isospace_plot")
# dev.off() # uncomment if pop ups are cumbersome
MixSIAR::plot_prior(alpha.prior=1,
                    source,plot_save_pdf=TRUE, plot_save_png=FALSE, "figs/prior_plot")

model_filename <- "MixSIAR_model.txt"
# When both resid_error and process_err are TRUE, the multiplicative error term
# is estimated, as described in Stock & Semmens 2016
resid_err <- TRUE # setting up error structure, residual true
process_err <- TRUE # setting up error structure, process true
MixSIAR::write_JAGS_model(model_filename, resid_err, process_err, mix, source)

# Run the described model and save as jags object
# run == 	   Chain Length  Burn-in 	 Thin 	# Chains  # NK Run Time
# “test” 	     1,000 	     500 	      1 	  3         Few seconds
# “very short” 10,000 	   5,000 	    5 	  3         Around 1-3 min
# “short” 	   50,000   	 25,000 	  25 	  3         Maybe 5 min or so
# “normal” 	   100,000 	   50,000 	  50 	  3         Varies but likely 10-15 min, used for analysis going forward
# “long” 	     300,000 	   200,000 	  100 	3         Untested
# “very long”  1,000,000 	 500,000 	  500 	3         Untested
# “extreme” 	 3,000,000 	 1,500,000  500 	3         Untested
jags.1 <- MixSIAR::run_model(run="very short", mix, source, discr, model_filename)

output_options <- list(summary_save = TRUE,                         # save summary statistics as txt file
                       summary_name = "output/summary_statistics",  # path for summary stats
                       sup_post = TRUE,                             # suppress posterior density plot, default is FALSE
                       plot_post_save_pdf = TRUE,                   # save posterior density plot
                       plot_post_name = "output/posterior_density", # path for posterior density plots
                       sup_pairs = TRUE,                            # suppress pairs plot
                       plot_pairs_save_pdf = TRUE,                  # save pairs plot
                       plot_pairs_name = "output/pairs_plot",       # path for pairs plot
                       sup_xy = TRUE,                               # suppress x/y trace plot
                       plot_xy_save_pdf = TRUE,                     # save x/y trace plot
                       plot_xy_name = "output/xy_plot",             # path for x/y trace plot
                       gelman = TRUE,                               # calculate Gelman-Rubin
                       heidel = TRUE,                               # calculate Heidelberg-Welch
                       geweke = TRUE,                               # calculate Geweke 
                       diag_save = TRUE,                            # save diagnostics as txt file
                       diag_name = "output/diagnostics",            # path for diagnostics
                       indiv_effect = FALSE,                        # artifact, set to FALSE
                       plot_post_save_png = FALSE,                  # save posterior density plot as png
                       plot_pairs_save_png = FALSE,                 # save pairs plot as png
                       plot_xy_save_png = FALSE,                    # save x/y trace plot as png
                       diag_save_ggmcmc = TRUE)                    # save ggmcmc diagnostics as pdf

# call function to save the ouputs of jags.1 following the predefined options
# input 1 Yes when prompted to safely save output in working directory
MixSIAR::output_JAGS(jags.1, mix, source, output_options) # safe to ignore warnings - NDK

# save jags.1 object for further analysis
saveRDS(jags.1,"output/jags.1.rds")
