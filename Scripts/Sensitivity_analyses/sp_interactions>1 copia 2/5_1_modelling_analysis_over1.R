############################################################ #
# Script name: "5_modelling_analysis.R"
# Objective: Fit statistical models relating pollinator metrics 
#            (e.g. PDI, normalise degree, Morisita-Horn overlap)
#            to ecological and climatic predictors
#
# Input:
#   - Data/modelling_data_int_over_one.rds     (final modeling dataset)
#
# Output:
#   
# Author: L. Marini
# Reviewer: J.B. Lanuza
# Year: 2025
############################################################ #
#### Load libraries ####
library(dplyr)
library(glmmTMB)
library(performance)
library(DHARMa)
library(effects)
library(car)
#### Load data ####
modelling_data = readRDS("Data/modelling_data_int_over_one.rds")
#### Load function ####
source("Scripts/export_table_function.R")
names(modelling_data)
#### Model 1:  SES(PDI) ####
fit_SES_pdi = glmmTMB(
  SES_pdi ~  Group*scale(Dominance) + Habitat + scale(Temperature) + scale(log(Plant))*scale(Dominance)+scale(I(Dominance^2))+
    (1|study_id/network_id/date),
  dispformula = ~scale(Dominance)*scale(log(Plant))+
    (1|study_id/network_id/date),
  data = subset(modelling_data, pdi_iqr>0 & pdi_sd>1e-6 & pdi_uniq>3)
)

summary(fit_SES_pdi)
#simulateResiduals(fit_SES_pdi, plot = TRUE, n=1000)
Anova(fit_SES_pdi)
#plot(allEffects(fit_SES_pdi))
check_collinearity(fit_SES_pdi)

#### Model 2:  SES(normalized degree) ####
names(modelling_data)

fit_SES_nd <- glmmTMB(
  SES_norm_degree ~  Group*scale(Dominance) + Habitat + scale(Temperature) + scale(log(Plant))*scale(Dominance)+scale(I(Dominance^2))+
    (1|study_id/network_id/date),
  dispformula = ~scale(Dominance)+scale(log(Plant))+
    (1|study_id/network_id/date),
  data = subset(modelling_data, norm_degree_iqr>0 & norm_degree_sd>1e-6 & norm_degree_uniq>3)
)

#simulateResiduals(fit_SES_nd, plot = TRUE, n=1000)
summary(fit_SES_nd)
Anova(fit_SES_nd)
#plot(allEffects(fit_SES_nd))
check_collinearity(fit_SES_nd)

#### Model 3: SES(Morisita–Horn) ####
names(modelling_data)
fit_SES_mh <- glmmTMB(
  SES_morisita_horn ~  Group*scale(Dominance)+scale(Plant) + Habitat + scale(Temperature) + scale(Plant)+scale(Dominance)+
    (1|study_id/network_id/date),
  dispformula = ~scale(Dominance)+scale(Plant)+
    (1|study_id/network_id/date),
  data = subset(modelling_data, morisita_horn_iqr>0 & morisita_horn_sd>1e-6 & morisita_horn_uniq>3)
)

#simulateResiduals(fit_SES_mh, plot = TRUE, n=1000)
summary(fit_SES_mh)
Anova(fit_SES_mh)
#plot(allEffects(fit_SES_mh))
check_collinearity(fit_SES_mh)

#### Save data ####
saveRDS(fit_SES_pdi, "Data/fit_SES_pdi_over1.rds")
saveRDS(fit_SES_nd, "Data/fit_SES_nd_over1.rds")
saveRDS(fit_SES_mh, "Data/fit_SES_mh_over1.rds")

