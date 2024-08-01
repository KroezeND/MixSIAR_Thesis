library(tidyverse)
library(here)
library(AICcmodavg)
library(lme4)

(rm(list = ls()))

traits <- read.csv("data/Biomass/Compiled_Traits.csv")
bgb_biomass_layer <- readRDS("output/bgb_biomass_layer.rds") #contains the mixing proportions and weight density

methane <- do.call(rbind,lapply(c("data/CH4/final_CH4_fluxes_2021-09-27.csv",
                                 "data/CH4/final_CH4_fluxes_2021-09-28.csv",
                                 "data/CH4/final_CH4_fluxes_2021-09-29.csv",
                                 "data/CH4/final_CH4_fluxes_2021-09-30.csv",
                                 "data/CH4/final_CH4_fluxes_2021-10-01.csv"),
                               function(file) read.csv(file, header = TRUE)))
methane <- methane[-c(1,34,71,104,130),]
methane$flux_CH4 <- as.numeric(methane$flux_CH4)
methane <- methane[methane$flux_CH4 > 0 & methane$flux_CH4 < 100000,]

methane$flux_CH4 <- methane$flux_CH4/24

# Use traits to set up holdout data, choosing holdout by pot_no not individual
methane_traits <- traits[traits$pot_id %in% methane$pot_id,]


merged_df <- full_join(methane_traits, methane, by = "pot_id")
merged_df <- left_join(merged_df,bgb_biomass_layer[["TOP"]],by = "pot_no")
merged_df <- merged_df[is.na(merged_df$prop_scam) ==F,]

ggplot(merged_df, aes(x = tot_bgb, y = flux_CH4))+
  geom_point()+
  geom_smooth(method = "lm",se = F)

model <- lm(flux_CH4 ~ tot_bgb, data = merged_df)
summary(model)

# mixture <- readRDS("output/mixture_list.rds")
# top <- mixture[[1]]
# 
# 
# # Model 1: Full
# model1 <- lmer(d13C ~ env_treatment + seed_year + (1|gt), data = top, REML = FALSE)
# 
# # Model 2: Random GT
# model2 <- lmer(d13C ~ (1|gt), data = top, REML = FALSE)
# 
# # Model 3: GT Seed Year
# model3 <- lmer(d13C ~ seed_year + (1|gt), data = top, REML = FALSE)
# 
# # Model 4: Env GT
# model4 <- lmer(d13C ~ env_treatment + (1|gt), data = top, REML = FALSE)
# 
# # Compare models using AIC
# model_list <- list(model1, model2, model3, model4)
# model_names <- c("Full", "Random GT", "GT:Seed Year", "Environment/GT")
# 
# aic_values <- aictab(cand.set = model_list, modnames = model_names)
# print(aic_values)
# 
# # Suppose model1 is the best model based on AIC
# summary(model1)
# 
# # Extracting variance components
# varcomp <- VarCorr(model1)
# print(varcomp)




# # Simplified model without random effect for gt
# simplified_model <- lmer(d13C ~ env_treatment + (1|gt), data = top)
# 
# # Summary of the simplified model
# summary(simplified_model)
# # Assessing the variance explained by the fixed effects
# library(car)
# 
# # ANOVA for the simplified model
# anova(simplified_model)
# 
# 
# (rm(list = ls()))