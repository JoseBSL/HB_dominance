############################################################ #
# Script name: "4_4_join_predictors.R"
# Objective: Join SES metrics with ecological predictors and prepare the modeling table
#
# Input:
#   - Data/pollinator_metrics_ses.rds
#   - Data/apis_network_metrics.rds
#   - Data/temperature.rds
#   - Data/ecological_predictors.rds  (contains: degree_interactions, species_data,
#                                      plants_per_network, habitat_net)
#
# Output: Data/modelling_data_int_over_one.rds
#
# Author: L. Marini
# Reviewer: J.B. Lanuza
# Year: 2025
############################################################ #
#### Load libraries ####
library(dplyr)
library(tidyr)
library(purrr)

#### Load data ####
pollinator_metrics_ses = readRDS("Data/pollinator_metrics_ses.rds")
apis_network_metrics = readRDS("Data/apis_network_metrics.rds")
temperature = readRDS("Data/temperature.rds")   # <- ensure this matches your file name
ecological_predictors = readRDS("Data/ecological_predictors.rds")

#Recover from list ecological predictors
degree_interactions = ecological_predictors$degree_interactions
species_data = ecological_predictors$species_data
plants_per_network = ecological_predictors$plants_per_network
habitat_net = ecological_predictors$habitat_net

#### Join data ####
merged_data = pollinator_metrics_ses %>%
  left_join(degree_interactions, by = c("id", "pollinator")) %>%
  left_join(species_data, by = "pollinator") %>%
  left_join(plants_per_network, by = "id") %>%
  left_join(temperature,by = "id") %>%
  left_join(habitat_net, by = "id") %>%
  left_join(apis_network_metrics,  by = "id")
str(merged_data) ############ dataframe entering the modeling

#### Join data ####
names(merged_data)
modelling_data = subset(merged_data,
                 Pollinator_yes!="no"   &   
                   sp_interactions >1     &
                   Apis_interactions>4    &
                   Pollinator_genus!="NA" &
                   sp_degree>0            &
                   n_plants>1)##

modelling_data = modelling_data %>%
  mutate(Group=as.factor(Pollinator_yes),
         Plant=plant_n,
         Habitat=habitat_type,
         id=as.factor(id),
         pollinator.x=as.factor(pollinator.x),
         .keep = "unused")

names(modelling_data)

#Separate identifier
modelling_data = modelling_data %>%
  separate(id, into = c("study_id", "network_id", "date"), sep = "//")

#### Save data ####
saveRDS(modelling_data, "Data/modelling_data_int_over_one.rds")

