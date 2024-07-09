# MixSIAR_Thesis
## Goals
The purpose of this analysis is to use the MixSIAR package to implement a Bayesian hierarchical mixing model to estimate the mixing proportion of two species of co-dominant marsh plants, _Schoenoplectus americanus_ and _Spartina patens_. As it is currently written, we estimate the mixing proportions using one fixed effect, environmental treatment, with four levels: high_sal/high_elev, high_sal/low_elev, low_sal/high_elev, low_sal/low_elev. We also include a continuous effect of seed year determined by <sup>210</sup>Pb dating of sediment layers. 

## Workflow
Files are organized and named in the order to be run sequentially (0-4).

- `0_data_assembly.R`: This script cleans the raw biomass and isotope data and formats each
     - Inputs:
        - `Belowground_Biomass.csv`                            Biomass spreadsheet from HSK Thesis
        - `BlueGenes2023_Belowground_Biomass - BG2023.csv`     Biomass spreadsheet from summer 2023 weighing effort
        - `Compiled_Traits.csv`                                Trait data from HSK Thesis
        - `220729_JPM.xls - CN sum.csv`                        Smithsonian lab isotope data
        - `220801_JPM.xls - CN sum.csv`                        Smithsonian lab isotope data
        - `220519_JPM.xls - CN sum.csv`                        Smithsonian lab isotope data
        - `CEST_Sample_Runs.csv`                               CEST lab isotope data - all runs, not including replicate testing
        - `Trial_Carbon_Data_unedited.csv`                     CEST lab isotope data - replicate trials
        - `Corr_Trial_COMP.csv`                                CEST lab isotope data - replicate trials
        - `Milled_replicates.csv`                              CEST lab isotope data - replicate trials
     - Outputs:
          - `marsh_source_all.csv`                             raw source isotopes formatted to MixSIAR specs, 151 observations of 2 variables
          - `Compiled_Traits_appended.csv`                     traits with added weights and root-shoot-ratios, 672 observations of 35 variables
          - `training.csv`                                     mixture distribution for training MixSIAR model, 274 observations of 17 variables
          - `bgb_biomass.csv`                                  all weights across each layer, 1692 observations of 11 variables
          - `holdout_DONOTTOUCH.csv`                           holdout data, 49 observations of 36 variables
- `1_run_model.R`: Describing model inputs and running model, adjust "run" in run_model for MCMC parameters
     - Inputs:
        - `training.csv`
        - `marsh_source_all.csv`
        - `marsh-discr`                                        null trophic enrichment factors
     - Output
        - `jags.Env.rds`                                       save full jags model run as RDS object
        - `figs/isospace_env_plot`                             d13C plot
        - `prior_env_plot`                                     uniform uninformative distribution
        - `diagnostics.pdf`                                    pdf for trace, density, and running means plot for each parameter
        - `diagnostics.txt`                                    text file containing Gelman-Rubin, Geweke, and Heidelberger-Welch diagnostics
        - `summary_statistics.txt`                             abbreviated file for main model output, group estimates, and DIC
        - `MixSIAR_model`                                      text file with model description
        - `posterior_density_high_sal`                         folder containing posterior density plots for group characteristics
        - `posterior_density_low_sal`                          folder containing posterior density plots for group characteristics
- `2_biomass_calc.R`: Calculate mixing proportion and weight densities, some plotting
     - Inputs:
          - `training.csv`
          - `jags.1.rds`
          - `bgb_biomass.csv`
     - Outputs:
          - `bgb_biomass_mod.rds`                                    modified to include mixing proportions and weight density
- `3_model_validation` temporary file holding space until model validation is pursued
     - Outputs: none
- `4_graphing.R`: Basic plotting functions for factors
     - Inputs:
          - `bgb_biomass.rds`
          - `Compiled_Traits_appended.csv`
     - Outputs: none

## File Structure
```
├── data                                  # raw biomass and isotope data
│   ├── Biomass                           # Biomass weights from two efforts
│   ├── Holdout                           # Folder to place holdout
│   └── Isotope                           # Isotope data from CEST and ND
├── figs                                  # folder for produced figures
├── output                                # intermediary files created during analysis
│   ├── data                              # modified data products
│   ├── posterior_density_high_salinity   # Empty folder to hold mixsiar outputs 
│   ├── posterior_density_high_salinity   # Empty folder to hold mixsiar outputs
│   └── jags.normal.rds                   # JAGS product used in analysis, takes a decently long time to run so provided here
├── src                                   # supporting functions
└── README.md
```

## R Version
All code was run in R version 4.4.0 (2024-04-24). Necessary packages to run script are identified at the top of each script.
- `dplyr` v1.1.4      
- `ggplot2` v3.5.1   
- `tidyverse` v2.0.0 
- `MixSIAR` v3.1.12 
